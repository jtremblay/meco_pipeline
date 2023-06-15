#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.design import *
from bio.readset import *

from bio import shotgun_metagenomics
from bio import microbial_ecology

from pipelines.illumina import illumina

# Global scope variables
log = logging.getLogger(__name__)

class MGS_augmented_assembly(illumina.Illumina):

    def trim(self):
        jobs = []
        # Merge all demultiplexed fastq files in one file. One file for reads1 and one file for reads2 of Illumina paired-end.
        # If library is Illumina single end or 454 or PacBio, only one file will be generated.
        outdir = self._root_dir
        trim_stats = []
        #sys.stderr.write('outdir: ' + outdir + '\n')
        
        for readset in self.readsets:
            trim_file_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".trim.")
            if readset.run_type == "PAIRED_END":
                #sys.stderr.write('outdir: ' + readset.fastq1 + '\n')
                if not os.path.exists(trim_file_prefix):
                    os.makedirs(trim_file_prefix)

                job = shotgun_metagenomics.trimmomatic(
                    readset.fastq1,
                    readset.fastq2,
                    trim_file_prefix + "pair1.fastq.gz",
                    trim_file_prefix + "single1.fastq.gz",
                    trim_file_prefix + "pair2.fastq.gz",
                    trim_file_prefix + "single2.fastq.gz",
                    readset.quality_offset,
                    trim_file_prefix + "out",
                    trim_file_prefix + "stats.csv"
                )

                trim_stats.append(trim_file_prefix + "stats.csv")

            elif readset.run_type == "SINGLE_END":
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (SINGLE_END not implemented yet)!")
            
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
            job.name = "trimmomatic." + readset.name
            job.subname = "trim"
            jobs.append(job)
        
        job = shotgun_metagenomics.merge_trimstats(
            trim_stats,
            os.path.join(outdir, "qced_reads", "trim_stats.tsv")
        )
        job.name = "merge_trim_stats"
        job.subname = "merge_trim_stats"
        jobs.append(job) 

        return jobs
            
    def duk(self):
        jobs=[]
        outdir = self._root_dir
        
        for readset in self.readsets:
            trim_file_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".trim.")
            outfile_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".")

            job = shotgun_metagenomics.duk_gz(
                trim_file_prefix + "pair1.fastq.gz",
                outfile_prefix + "contam_R1.fastq",
                outfile_prefix + "ncontam_R1.fastq.gz",
                outfile_prefix + "duk_contam_R1_log.txt",
                config.param('DB', 'contaminants', 1, 'filepath')
            )
            job.name = "duk_R1_" + readset.sample.name
            job.subname = "duk"
            jobs.append(job)
        
            job = shotgun_metagenomics.duk_gz(
                trim_file_prefix + "pair2.fastq.gz",
                outfile_prefix + "contam_R2.fastq",
                outfile_prefix + "ncontam_R2.fastq.gz",
                outfile_prefix + "duk_contam_R2_log.txt",
                config.param('DB', 'contaminants', 1, 'filepath')
            )
            job.name = "duk_R2_" + readset.sample.name
            job.subname = "duk"
            jobs.append(job)

            job = shotgun_metagenomics.merge_pairs(
                outfile_prefix + "ncontam_R1.fastq.gz",
                outfile_prefix + "ncontam_R2.fastq.gz",
                outfile_prefix + "ncontam_paired.fastq.gz"
            )
            job.name = "merge_pairs_" + readset.sample.name
            job.subname = "merge_pairs"
            jobs.append(job) 

        return jobs
     
    def first_assembly(self):
        root = self._root_dir
        jobs = []
        fastq_list = []
        contigs_list = []
        contigs_out_list = []
         
        job = shotgun_metagenomics.ray_big(
            fastq_list,
            os.path.join(root, "first_assembly", "big_assembly"),
            "pe"
        )
        job.name = "ray_big_assembly"
        job.subname = "ray_big"
        jobs.append(job)
        
        job = shotgun_metagenomics.compile_ray(
            os.path.join(root, "first_assembly", "big_assembly", "Contigs.fasta"),
            os.path.join(root, "first_assembly", "big_assembly", "ray_assembly_stats.txt")
        )
        job.name = "compile_ray_big_assembly"
        job.subname = "compile_ray"
        jobs.append(job)

        return jobs

    def gene_prediction(self):
        root = self._root_dir
        jobs = []
        gff_list = []
        fna_list = []
        
        infile = os.path.join(root, "first_assembly", "big_assembly", "Contigs.fasta")
        outdir = os.path.join(root, "gene_prediction", "big_assembly")
        #gff = os.path.join(root, "gene_prediction", "big_assembly", "Contigs.gff")
        #fna = os.path.join(root, "gene_prediction", "big_assembly", "Contigs.fna")
        
        if not os.path.exists(outdir):
            os.makedirs(os.path.join(outdir))
        
        job = shotgun_metagenomics.metagenemark(
            infile,
            os.path.join(outdir, "Contigs.gff"),
            os.path.join(outdir, "Contigs.fna"),
            os.path.join(outdir, "Contigs.faa"),
            os.path.join(outdir, "Contigs_renamed.fna")
        )
        job.name = "metagenemark_big_assembly"
        job.subname = "metagenemark_big"
        jobs.append(job)
         
        return jobs

#    def gene_clustering(self):
#        jobs = []
#        
#        infiles = []
#        for readset in self.readsets:
#            infiles.append(os.path.join(self._root_dir, "gene_prediction", readset.sample.name, readset.sample.name + "_renamed.fna"))
#
#        merged_fasta = os.path.join(self._root_dir, "gene_clustering", "genes.fasta")
#        barcodes = os.path.join(self._root_dir, "gene_clustering", "barcodes.fasta")
#
#        job = shotgun_metagenomics.preprocess_for_clustering(
#            infiles,
#            merged_fasta,
#            barcodes
#        )
#        job.name = "preprocess_for_clustering"
#        job.subname = "preprocess_for_clustering"
#        jobs.append(job)
#        
#        job = shotgun_metagenomics.clustering(
#            merged_fasta,
#            barcodes,
#            os.path.join(self._root_dir, "gene_clustering")
#        )   
#        job.name = "generate_clusters"
#        job.subname = "clustering"
#        jobs.append(job)
#        
#        job = shotgun_metagenomics.filter_clusters(
#            os.path.join(self._root_dir, "gene_clustering", "cdhit.fasta"),
#            os.path.join(self._root_dir, "gene_clustering", "cdhit_filtered.fasta")
#        )   
#        job.name = "filter_clusters"
#        job.subname = "filter_clusters"
#        jobs.append(job)
#        
#        return jobs

    ## Gene abundance and canopy clustering

    def abundance(self):
        jobs = []
        cov_list_contigs = []
        cov_list_genes = []
        flagstats_contigs_list = []
        trimmomatic_list = []
        
        # Will make index for bwa. and also bed file for computing reads spanning later.
        reference_contigs = os.path.join(self._root_dir, "first_assembly", "big_assembly", "Contigs.fasta")
        reference_genes = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "Contigs.fna")

        job = shotgun_metagenomics.make_index(
            reference_contigs
        )
        job.name = "make_index_contigs"
        job.subname = "make_index"
        jobs.append(job)
        
        job = shotgun_metagenomics.make_index(
            reference_genes
        )
        job.name = "make_index_genes"
        job.subname = "make_index"
        jobs.append(job)
            
        bed_contigs = os.path.join(self._root_dir, "first_assembly", "big_assembly", "Contigs.bed")
        bwt_contigs = os.path.join(self._root_dir, "first_assembly", "big_assembly", "Contigs.fasta.bwt")
        bed_genes = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "Contigs.bed")
        bwt_genes = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "Contigs.fna.bwt")

        for readset in self.readsets:
            bam_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".bam")
            flagstats_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".flagstats")
            flagstats_contigs_list.append(flagstats_contigs)
            trimmomatic = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".trim.stats.csv")
            trimmomatic_list.append(trimmomatic)
            cov_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".cov")
            cov_list_contigs.append(cov_contigs)
            bam_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".bam")
            cov_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".cov")
            cov_list_genes.append(cov_genes)
    
            outdir = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name)
            outdir_genes = os.path.join(self._root_dir, "genes_abundance", readset.sample.name)
            if not os.path.exists(outdir):
                os.makedirs(os.path.join(outdir))
            if not os.path.exists(outdir_genes):
                os.makedirs(os.path.join(outdir_genes))
       
            if( readset.run_type == "PAIRED_END"):
                infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
                out1 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
                out2 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
                
                #job = shotgun_metagenomics.split_pairs(
                #    infile,
                #    out1,
                #    out2
                #)
                #job.name = "split_pairs_" + readset.sample.name
                #job.subname = "split_pairs"
                #jobs.append(job)

                # map against contigs
                job = shotgun_metagenomics.bwa_mem_samtools(
                    reference_contigs,
                    out1,
                    out2,
                    bam_contigs
                )
                job.name = "bwa_mem-contigs" + readset.sample.name
                job.subname = "bwa"
                jobs.append(job)
                
                job = shotgun_metagenomics.flagstats(
                    bam_contigs,
                    flagstats_contigs
                )
                job.name = "flagstats" + readset.sample.name
                job.subname = "flagstats"
                jobs.append(job)
                
                # map against genes for DDA.
                job = shotgun_metagenomics.bwa_mem(
                    reference_genes,
                    out1,
                    out2,
                    bam_genes
                )
                job.name = "bwa_mem-genes" + readset.sample.name
                job.subname = "bwa"
                jobs.append(job)

            elif(readset.run_type == "SINGLE_END"):
                infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_R1.fastq.gz")
                
                job = shotgun_metagenomics.bwa_mem_se(
                    reference,
                    infile,
                    bam
                )
                job.name = "bwa_mem-" + readset.sample.name
                job.subname = "bwa_short"
                jobs.append(job)

            job = shotgun_metagenomics.coverage_bed(
                bam_contigs,
                bed_contigs,
                cov_contigs
            )
            job.name = "bedtoolsCov-contigs" + readset.sample.name
            job.subname = "bedtools"
            jobs.append(job)
            
            job = shotgun_metagenomics.coverage_bed(
                bam_genes,
                bed_genes,
                cov_genes
            )
            job.name = "bedtoolsCov-genes" + readset.sample.name
            job.subname = "bedtools"
            jobs.append(job)
            
        # Once all coverage has been computed, merge all tables.
        #sys.stderr.write('gene abundance: ' + ','.join([str(x) for x in cov_list] ) + '\n')
        job = shotgun_metagenomics.merge_counts(
            cov_list_contigs,
            os.path.join(self._root_dir, "contigs_abundance", "merged_contigs_abundance.tsv"),
            os.path.join(self._root_dir, "contigs_abundance", "merged_contigs_abundance_RPKM.tsv")
        )
        job.name = "merge_gene_abundance_contigs"
        job.subname = "merge_gene_abundance"
        jobs.append(job)
        
        job = shotgun_metagenomics.merge_counts(
            cov_list_genes,
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance_RPKM.tsv")
        )
        job.name = "merge_gene_abundance_genes"
        job.subname = "merge_gene_abundance"
        jobs.append(job)
                
        job = shotgun_metagenomics.merge_flagstats(
            flagstats_contigs_list,
            trimmomatic_list,
            os.path.join(self._root_dir, "contigs_abundance", "qc_mapping_stats.tsv")
        )
        job.name = "flagstats_merge"
        job.subname = "flagstats"
        jobs.append(job)
         
        return jobs 
    
    def gene_binning(self):
        jobs = []
        
        #job = shotgun_metagenomics.canopy(
        #    os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
        #    "canopy",
        #    os.path.join(self._root_dir, "gene_binning", "clusters.tsv"),
        #    os.path.join(self._root_dir, "gene_binning", "profiles.tsv")
        #)
        #job.name = "canopy"
        #job.subname = "canopy"
        #jobs.append(job)

        return jobs
    
    def exonerate(self):
        jobs = []

        for readset in self.readsets:
            infile_fna = os.path.join(self._root_dir, "gene_prediction", readset.sample.name, readset.sample.name + ".fna")
            infile_faa = os.path.join(self._root_dir, "gene_prediction", readset.sample.name, readset.sample.name + ".faa")
            chunks_dir = os.path.join(self._root_dir, "gene_prediction", readset.sample.name, "fasta_chunks")
            num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
            
        # Split for big assembly
        infile_fna = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "Contigs.fna")
        infile_faa = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "Contigs.faa")
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "fasta_chunks")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
        
        #FNA
        job = shotgun_metagenomics.exonerate(
            infile_fna,
            chunks_dir,
            num_chunks,
            "Contigs.fna" 
        )
        job.name = "exonerate_big_assembly_fna"
        job.subname = "exonerate"
        jobs.append(job)
        
        # FAA
        job = shotgun_metagenomics.exonerate(
            infile_faa,
            chunks_dir,
            num_chunks,
            "Contigs.faa" 
        )
        job.name = "exonerate_big_assembly_faa"
        job.subname = "exonerate"
        jobs.append(job)
 
        return jobs    

    def blastn_nt(self):
        jobs = []
        
        # Do blastn on nt for big assembly  
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastn_nt")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')

        for i in range(num_chunks):
            job = shotgun_metagenomics.blastn(
                os.path.join(chunks_dir, "Contigs.fna_chunk_{:07d}".format(i)),
                os.path.join(blast_dir, "blastn_chunk_{:07d}.tsv".format(i)),
                blast_dir 
            )
            job.name = "blastn_nt_big_assembly"
            job.subname = "blastn"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastn_nt.tsv"),
            num_chunks,
            "blastn" 
        )
        job.name = "blastn_nt_big_assembly_merge"
        job.subname = "blastn"
        jobs.append(job)
        
        job = shotgun_metagenomics.keep_blast_best_hit(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastn_nt.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastn_nt_besthit.tsv")
        )
        job.name = "blastn_nt_big_assembly_best_hit"
        job.subname = "blastn_best_hit"
        jobs.append(job)
        
        return jobs
    
    def blastp_kegg(self):
        # Blastx cog

        jobs = []

        # Do blastp on KEGG for big assembly  
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastp_kegg")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')

        for i in range(num_chunks):
            job = shotgun_metagenomics.blastp(
                os.path.join(chunks_dir, "Contigs.faa_chunk_{:07d}".format(i)),
                os.path.join(blast_dir, "blastp_chunk_{:07d}.tsv".format(i)),
                blast_dir,
                config.param('kegg', 'db') 
            )
            job.name = "blastp_kegg_big_assembly"
            job.subname = "blastp"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastp_kegg.tsv"),
            num_chunks,
            "blastp"
        )
        job.name = "blastp_kegg_big_assembly_merge"
        job.subname = "merge"
        jobs.append(job)
        
        # Generate a clean table of module/pathways.
        job = shotgun_metagenomics.parse_kegg(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastp_kegg.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastp_kegg_parsed.tsv")
        )
        job.name = "parse_kegg"
        job.subname = "merge"
        jobs.append(job)

        return jobs 
    
    def rpsblast_cog(self):
        jobs = []
        
        # Do rpsblast on COG for big assembly  
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "big_assembly", "rpsblast_cog")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')

        for i in range(num_chunks):
            job = shotgun_metagenomics.rpsblast(
                os.path.join(chunks_dir, "Contigs.faa_chunk_{:07d}".format(i)),
                os.path.join(blast_dir, "rpsblast_chunk_{:07d}.tsv".format(i)),
                blast_dir,
                config.param('cog', 'db') 
            )
            job.name = "rpsblast_cog_big_assembly"
            job.subname = "rpsblast"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "rpsblast_cog.tsv"),
            num_chunks,
            "rpsblast"
        )
        job.name = "rpsblast_cog_big_assembly_merge"
        job.subname = "merge"
        jobs.append(job)
        
        return jobs 
    
    def rpsblast_kog(self):
        jobs = []
        
        # Do rpsblast on COG for big assembly  
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "big_assembly", "rpsblast_kog")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')

        for i in range(num_chunks):
            job = shotgun_metagenomics.rpsblast(
                os.path.join(chunks_dir, "Contigs.faa_chunk_{:07d}".format(i)),
                os.path.join(blast_dir, "rpsblast_chunk_{:07d}.tsv".format(i)),
                blast_dir,
                config.param('kog', 'db') 
            )
            job.name = "rpsblast_kog_big_assembly"
            job.subname = "rpsblast"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "rpsblast_kog.tsv"),
            num_chunks,
            "rpsblast"
        )
        job.name = "rpsblast_kog_big_assembly_merge"
        job.subname = "merge"
        jobs.append(job)
        
        return jobs 
    
    def hmmscan_tigrfam(self):
        jobs = [] # Do rpsblast on COG for big assembly
 
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "big_assembly", "hmmscan_tigrfam")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint', required=True)

        for i in range(num_chunks):
            job = shotgun_metagenomics.hmmscan(
                os.path.join(chunks_dir, "Contigs.faa_chunk_{:07d}".format(i)),
                os.path.join(blast_dir, "hmmscan_chunk_{:07d}".format(i)),
                blast_dir,
                config.param('tigrfam', 'db', required=True) 
            )
            job.name = "hmmscan_tigrfam_big_assembly"
            job.subname = "hmmscan"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks_hmms(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "big_assembly"),
            num_chunks,
            "hmmscan",
            "hmmscan_tigrfam"
        )
        job.name = "hmmscan_tigrfam_big_assembly_merge"
        job.subname = "merge"      
        jobs.append(job)

        return jobs 
    
    def hmmscan_pfam(self):
        jobs = []
        
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "big_assembly", "hmmscan_pfam")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint', required=True)

        for i in range(num_chunks):
            job = shotgun_metagenomics.hmmscan(
                os.path.join(chunks_dir, "Contigs.faa_chunk_{:07d}".format(i)),
                os.path.join(blast_dir, "hmmscan_chunk_{:07d}".format(i)),
                blast_dir,
                config.param('pfam', 'db', required=True) 
            )
            job.name = "hmmscan_pfam_big_assembly"
            job.subname = "hmmscan"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks_hmms(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "big_assembly"),
            num_chunks,
            "hmmscan",
            "hmmscan_pfam"
        )
        job.name = "hmmscan_pfam_big_assembly_merge"
        job.subname = "merge"      
        jobs.append(job)
        
        return jobs
    
    def taxonomic_annotation(self):
        jobs = []
        
        job = shotgun_metagenomics.extract_taxonomy(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastn_nt_besthit.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy.tsv"),
        )
        job.name = "extract_taxonomy_big_assembly"
        job.subname = "blastn_best_hit"
        jobs.append(job)

        job = shotgun_metagenomics.generate_otu_table(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy.tsv"),
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table.tsv")
        )
        job.name = "generate_otu_table"
        job.subname = "generate_otu_table"
        jobs.append(job)
        
        job = microbial_ecology.convert_otu_to_biom(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table.biom")
        )
        job.name = "convert_otu_to_biom"
        job.subname = "convert_otu_table"
        jobs.append(job)

        job = microbial_ecology.split_otu_table(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_bacteriaArchaea.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_others.tsv"),
            "bactArch"
        )
        job.name = "split_otu_table_bactArch"
        job.name = "split_otu_table"
        jobs.append(job)

        job = microbial_ecology.convert_otu_to_biom(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_bacteriaArchaea.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_bacteriaArchaea.biom")
        )
        job.name = "convert_otu_table_to_biom_bactArch"
        job.subname = "convert_otu_table"
        jobs.append(job)
        
        job = microbial_ecology.convert_otu_to_biom(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_others.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_others.biom")
        )
        job.name = "convert_otu_table_to_biom_others"
        job.subname = "convert_otu_table"
        jobs.append(job)

        job = microbial_ecology.rarefy(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table.biom"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_rarefied.biom"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_rarefied.tsv")
        )    
        job.name = "rarefy_otu_table_all_organisms"
        job.subname = "rarefaction"
        jobs.append(job)
        
        job = microbial_ecology.rarefy(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_others.biom"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_others_rarefied.biom"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_others_rarefied.tsv")
        )    
        job.name = "rarefy_otu_table_others"
        job.subname = "rarefaction"
        jobs.append(job)
        
        job = microbial_ecology.rarefy(
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_bacteriaArchaea.biom"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_bacteriaArchaea_rarefied.biom"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_bacteriaArchaea_rarefied.tsv")
        )    
        job.name = "rarefy_otu_table_bactArchaea"
        job.subname = "rarefaction"
        jobs.append(job)
        
        # ALL
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_absolute(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "all", "absolute")
            )
            job.name = "tax_summary_absolute_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
        
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_relative(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "all", "relative")
            )
            job.name = "tax_summary_relative_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
        
        # OTHERS
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_absolute(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_others_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "others", "absolute")
            )
            job.name = "tax_summary_absolute_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
        
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_relative(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_others_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "others", "relative")
            )
            job.name = "tax_summary_relative_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
        
        # BACTERIA-ARCHAEA
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_absolute(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_bacteriaArchaea_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "bacteriaArchaea", "absolute")
            )
            job.name = "tax_summary_absolute_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
        
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_relative(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "otu_table_bacteriaArchaea_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "bacteriaArchaea", "relative")
            )
            job.name = "tax_summary_relative_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
            
        # Plot taxa - ALL
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "all", "relative", "otu_table_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "all", "relative", "plots"),
                "taxonomy_L" + str(i) + "_relative"
            )
            job.name = "plot_taxa_single_relative_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
            
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "all", "absolute", "otu_table_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "all", "absolute", "plots"),
                "taxonomy_L" + str(i) + "_absolute"
            )
            job.name = "plot_taxa_single_absolute_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
        
        # OTHERS
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "others", "relative", "otu_table_others_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "others", "relative", "plots"),
                "taxonomy_L" + str(i) + "_relative"
            )
            job.name = "plot_taxa_single_relative_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
            
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "others", "absolute", "otu_table_others_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "others", "absolute", "plots"),
                "taxonomy_L" + str(i) + "_absolute"
            )
            job.name = "plot_taxa_single_absolute_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
        
        # Bacteria-Archaea
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "bacteriaArchaea", "relative", "otu_table_bacteriaArchaea_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "bacteriaArchaea", "relative", "plots"),
                "taxonomy_L" + str(i) + "_relative"
            )
            job.name = "plot_taxa_single_relative_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
            
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "bacteriaArchaea", "absolute", "otu_table_bacteriaArchaea_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy", "bacteriaArchaea", "absolute", "plots"),
                "taxonomy_L" + str(i) + "_absolute"
            )
            job.name = "plot_taxa_single_absolute_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)

        return jobs
    
    def statistics(self):
        jobs = []
        
        job = shotgun_metagenomics.edger_single(
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "DDA", "DDA_normalized.tsv")
        )
        job.name = "DDA_single"
        job.subname = "DDA_single"
        jobs.append(job)
        
        job = shotgun_metagenomics.edger(
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            config.param('DEFAULT', 'designFile'),
            os.path.join(self._root_dir, "DDA")
        )
        job.name = "DDA"
        job.subname = "DDA"
        jobs.append(job)
        
            
        job = shotgun_metagenomics.draw_pcoa(
            os.path.join(self._root_dir, "DDA"),
            config.param('DEFAULT', 'mappingFile'),
            "_normalized.tsv",
            os.path.join(self._root_dir, "DDA", "DDA.done"),
            os.path.join(self._root_dir, "DDA", "DDA_normalized.tsv")
        )
        job.name = "drawPcoa"
        job.subname = "drawPcoa"
        jobs.append(job)
       
        types = ["pathways", "modules", "K"]
        for type in types:
            job = shotgun_metagenomics.kegg_overrep(
                os.path.join(self._root_dir, "DDA"),
                os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastp_kegg_parsed.tsv"),
                type,
                os.path.join(self._root_dir, "DDA", "DDA.done")
            )
            job.name = "getKegg_" + type
            job.subname = "getKegg"
            jobs.append(job)
        
        # COG and KOG
        job = shotgun_metagenomics.cog_overrep(
            os.path.join(self._root_dir, "DDA"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "rpsblast_cog.tsv"),
            os.path.join(self._root_dir, "DDA", "DDA.done")
        )
        job.name = "getCOG"
        job.subname = "getCOG"
        jobs.append(job)
        
        job = shotgun_metagenomics.kog_overrep(
            os.path.join(self._root_dir, "DDA"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "rpsblast_kog.tsv"),
            os.path.join(self._root_dir, "DDA", "DDA.done")
        )
        job.name = "getKOG"
        job.subname = "getKOG"
        jobs.append(job)
        
        # Pfam and tigrfam
        job = shotgun_metagenomics.pfam_overrep(
            os.path.join(self._root_dir, "DDA"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "hmmscan_pfam_tblout.tsv"),
            os.path.join(self._root_dir, "DDA", "DDA.done")
        )
        job.name = "getPfam"
        job.subname = "getPfam"
        jobs.append(job)
        
        job = shotgun_metagenomics.tigrfam_overrep(
            os.path.join(self._root_dir, "DDA"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "hmmscan_tigrfam_tblout.tsv"),
            os.path.join(self._root_dir, "DDA", "DDA.done")
        )
        job.name = "getTigrfam"
        job.subname = "getTigrfam"
        jobs.append(job)
        
        # heatmaps
        prefixes = ["DDA_normalized_significant_keggK.tsv", "DDA_normalized_significant_keggModules.tsv", "DDA_normalized_significant_keggPathways.tsv", "DDA_normalized_significant_COG.tsv", "DDA_normalized_significant_KOG.tsv", "DDA_normalized_significant_PFAM.tsv", "DDA_normalized_significant_TIGRFAM.tsv"]
        for prefix in prefixes:
            job = shotgun_metagenomics.draw_heatmap(
                os.path.join(self._root_dir, "DDA"),
                config.param('DEFAULT', 'mappingFile'),
                prefix,
                os.path.join(self._root_dir, "DDA", "DDA.done")
            )
            job.name = "draw_heatmap"
            job.subname = "draw_heatmap_" +  prefix
            jobs.append(job)
             
        return jobs
   
    # Here generate final GFF (for viewing data in a genome browser). 
    # And generate final DDA sheets. with logFC, gene_name and actual normalized values.
    def finalize(self):
        jobs = []
        
        job = shotgun_metagenomics.generate_gff(
            os.path.join(self._root_dir, "gene_prediction", "big_assembly", "Contigs.gff"),
            os.path.join(self._root_dir, "first_assembly", "big_assembly", "Contigs.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastn_nt_besthit.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastp_kegg_parsed.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "hmmscan_pfam_tblout.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "hmmscan_tigrfam_tblout.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "rpsblast_cog.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "rpsblast_kog.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "taxonomy.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "Contigs.gff"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "Contigs.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "annotations.tsv")
        )
        job.name = "generate_gff"
        job.subname = "generate_gff"
        jobs.append(job)
        
        job = shotgun_metagenomics.annotate_DDA(
            os.path.join(self._root_dir, "DDA"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "annotations.tsv")
        )
        job.name = "annotate_DDA"
        job.subname = "annotate_DDA"
        jobs.append(job)

        return jobs

    def generate_report(self):
        jobs = []
        
        i = datetime.datetime.now()
        year = "%04d" %i.year
        day = "%02d" %i.day
        month = "%02d" %i.month
        title = config.param('report', 'report.title')
        date = year + "_" + month + "_" + day + "_metagenomics_" + title
        
        job = shotgun_metagenomics.generate_report(
            self._config,
            self._root_dir,
            "Metagenomics",
            self._root_dir + "/" + date
        )
        job.name = "deliverables"
        job.subname = "report"
        jobs.append(job)
        return jobs
    
    def groopm(self):
        jobs = []
        
        bams = []
        for readset in self.readsets:
            bam_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".bam")
            bams.append(bam_contigs)
        
        job = shotgun_metagenomics.groopm_parse(
            os.path.join(self._root_dir, "first_assembly", "big_assembly", "Contigs.fasta"),
            os.path.join(self._root_dir, "binning", "db.gm"),
            os.path.join(self._root_dir, "binning", "groopm_parse.done"),
            bams
        )
        job.name = "groopm_parse"
        job.subname = "groopm_parse"
        jobs.append(job)
        
        job = shotgun_metagenomics.groopm_core(
            os.path.join(self._root_dir, "binning", "db.gm"),
            os.path.join(self._root_dir, "binning", "groopm_parse.done"),
            os.path.join(self._root_dir, "binning", "groopm_core.done")
        )
        job.name = "groopm_core"
        job.subname = "groopm_core"
        jobs.append(job)
        
        job = shotgun_metagenomics.groopm_recruit(
            os.path.join(self._root_dir, "binning", "db.gm"),
            os.path.join(self._root_dir, "binning", "groopm_core.done"),
            os.path.join(self._root_dir, "binning", "groopm_recruit.done")
        )
        job.name = "groopm_core"
        job.subname = "groopm_core"
        jobs.append(job)

        return jobs

    def cleanup(self):
        #Here, compress all .fastq files into .fastq.gz.
        jobs = []
        #job = shotgun_metagenomics.mymethod(
        #)
        #job.name = "myjobname"
        #jobs.append(job)
        sys.stderr.write('not implemented yet\n')
        return jobs
    
    # Override illumina.py readsets to make sure we are parsing a nanuq sample sheet
    # and not a readset sheet.
    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_nanuq_readset_file(self.args.readsets.name)
        return self._readsets

    @property
    def steps(self):
        # external infile?
        return [
            self.trim,
            self.duk,
            self.first_assembly,
            self.gene_prediction,
            self.abundance,
            self.exonerate,#6
            self.blastn_nt,#7
            self.blastp_kegg,#8
            self.rpsblast_cog,#9
            self.rpsblast_kog,#10
            self.hmmscan_tigrfam,#11
            self.hmmscan_pfam,#12
            self.taxonomic_annotation,#13
            #self.gene_binning,
            self.statistics, #14
            self.finalize, #15
            self.generate_report, #16
            self.groopm, #17
            self.cleanup
        ]

    def set_local_variables(self):
        self._parser_local = argparse.ArgumentParser(description='Process options.')
        self._parser_local.add_argument("-c", "--config", help="config INI-style file", nargs="+", type=file, required=True)
        self._parser_local.add_argument("-r", "--readsets", help="readset file", type=file, required=False)
        self._parser_local.add_argument("-s", "--steps", help="step range e.g. '1-5', '3,6,7', '2,4-8'", required=True)
        self._parser_local.add_argument("-o", "--output-dir", help="output directory (default: current)", default=os.getcwd())
        self._parser_local.add_argument("-j", "--job-scheduler", help="job scheduler type (default: torque)", choices=["torque", "batch", "daemon"], default="torque")
        self._parser_local.add_argument("-f", "--force", help="force creation of jobs even if up to date (default: false)", action="store_true")
        self._parser_local.add_argument("-l", "--log", help="log level (default: info)", choices=["debug", "info", "warning", "error", "critical"], default="info")

        # barcodes
        self._args_local = self._parser_local.parse_args()
        config.parse_files(self._args_local.config)
        self._config = self._args_local.config[0]
        self._config = os.path.abspath(self._config.name)
        
        self._root_dir = self._args_local.output_dir
        if not os.path.exists(self._root_dir):
            os.makedirs(self._root_dir)
        
        # Make directories
        self.make_directories(self._root_dir)
  
    # Define and make directories. Also desing initial infile.
    def make_directories(self, root_dir):
         
        if not os.path.exists(root_dir):
            os.makedirs(root_dir)
        if not os.path.exists(root_dir + "/DDA"):
            os.makedirs(os.path.join(root_dir, "DDA"))
        if not os.path.exists(root_dir + "/qced_reads"):
            os.makedirs(os.path.join(root_dir, "qced_reads"))
        if not os.path.exists(root_dir + "/first_assembly"):
            os.makedirs(os.path.join(root_dir, "first_assembly"))
        if not os.path.exists(root_dir + "/gene_prediction"):
            os.makedirs(os.path.join(root_dir, "gene_prediction"))
        if not os.path.exists(root_dir + "/gene_clustering"):
            os.makedirs(os.path.join(root_dir, "gene_clustering"))
        if not os.path.exists(root_dir + "/gene_abundance"):
            os.makedirs(os.path.join(root_dir, "gene_abundance"))
        if not os.path.exists(root_dir + "/contigs_abundance"):
            os.makedirs(os.path.join(root_dir, "contigs_abundance"))
        if not os.path.exists(root_dir + "/taxonomic_annotation"):
            os.makedirs(os.path.join(root_dir, "taxonomic_annotation"))
        if not os.path.exists(root_dir + "/gene_annotation"):
            os.makedirs(os.path.join(root_dir, "gene_annotation"))
        if not os.path.exists(root_dir + "/gene_annotation/big_assembly/taxonomy/all/absolute/plots"):
            os.makedirs(os.path.join(root_dir, "gene_annotation/big_assembly/taxonomy/all/absolute/plots"))
            os.makedirs(os.path.join(root_dir, "gene_annotation/big_assembly/taxonomy/all/relative/plots"))
        if not os.path.exists(root_dir + "/gene_annotation/big_assembly/taxonomy/others/absolute/plots"):
            os.makedirs(os.path.join(root_dir, "gene_annotation/big_assembly/taxonomy/others/absolute/plots"))
            os.makedirs(os.path.join(root_dir, "gene_annotation/big_assembly/taxonomy/others/relative/plots"))
        if not os.path.exists(root_dir + "/gene_annotation/big_assembly/taxonomy/bacteriaArchaea/absolute/plots"):
            os.makedirs(os.path.join(root_dir, "gene_annotation/big_assembly/taxonomy/bacteriaArchaea/absolute/plots"))
            os.makedirs(os.path.join(root_dir, "gene_annotation/big_assembly/taxonomy/bacteriaArchaea/relative/plots"))
        if not os.path.exists(root_dir + "/gene_abundance_clustering"):
            os.makedirs(os.path.join(root_dir, "gene_abundance_clustering"))
        if not os.path.exists(root_dir + "/binning"):
            os.makedirs(os.path.join(root_dir, "binning"))
        if not os.path.exists(root_dir + "/statistics"):
            os.makedirs(os.path.join(root_dir, "statistics"))


    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            #self._readsets = parse_nanuq_readset_file(self.args.readsets.name)
            self._readsets = parse_readset_file(self.args.readsets.name)
        return self._readsets


    def __init__(self):
        # Add pipeline specific arguments
        self.set_local_variables()
        sys.stderr.write('Running MGS (Metagenomic Genomic Specie) augmented pipeline\n')
        super(MGS_augmented_assembly, self).__init__()
                
MGS_augmented_assembly().submit_jobs()
