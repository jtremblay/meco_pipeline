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

from pipelines.illumina import illumina

# Global scope variables
log = logging.getLogger(__name__)

class Quick_assembly(illumina.Illumina):

    def trim(self):
        jobs = []
        # Merge all demultiplexed fastq files in one file. One file for reads1 and one file for reads2 of Illumina paired-end.
        # If library is Illumina single end or 454 or PacBio, only one file will be generated.
        outdir = self._root_dir
        #sys.stderr.write('outdir: ' + outdir + '\n')
        
        for readset in self.readsets:
            trim_file_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".trim.")
            if not os.path.exists(trim_file_prefix):
                os.makedirs(trim_file_prefix)
            
            if readset.run_type == "PAIRED_END":
                #sys.stderr.write('outdir: ' + readset.fastq1 + '\n')

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
            elif readset.run_type == "SINGLE_END":
                job = shotgun_metagenomics.trimmomatic_se(
                    readset.fastq1,
                    trim_file_prefix + "single1.fastq.gz",
                    readset.quality_offset,
                    trim_file_prefix + "out",
                    trim_file_prefix + "stats.csv"
                )
            
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
            job.name = "trimmomatic." + readset.name
            job.subname = "trim"
            jobs.append(job)
         
        return jobs
            
    def duk(self):
        jobs=[]
        outdir = self._root_dir
        
        for readset in self.readsets:
            trim_file_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".trim.")
            outfile_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".")

            if(readset.run_type == "SINGLE_END"):
                job = shotgun_metagenomics.duk_gz(
                    trim_file_prefix + "single1.fastq.gz",
                    outfile_prefix + "contam_R1.fastq",
                    outfile_prefix + "ncontam_R1.fastq.gz",
                    outfile_prefix + "duk_contam_R1_log.txt",
                    config.param('DB', 'contaminants', 1, 'filepath')
                )
                job.name = "duk_R1_" + readset.sample.name
                job.subname = "duk"
                jobs.append(job)
            
            elif(readset.run_type == "PAIRED_END"):
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
    
     
    def assembly(self):
        root = self._root_dir
        jobs = []

        fastq_list = []
        
        for readset in self.readsets:
            outdir = root + "assembly/" + readset.sample.name
             
            if(readset.run_type == "SINGLE_END"):
                infile = os.path.join(root, "qced_reads", readset.sample.name, readset.name + ".ncontam_R1.fastq.gz")
                type = "se"

            elif(readset.run_type == "PAIRED_END"):
                infile = os.path.join(root, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
                type = "pe"
            
            fastq_list.append(infile)

            if(self._big_assembly_only == False):
                job = shotgun_metagenomics.ray(
                    infile,
                    outdir,
                    type
                )
                job.name = "ray_" + readset.sample.name
                job.subname = "ray"
                jobs.append(job)
        

        job = shotgun_metagenomics.ray_big(
            fastq_list,
            os.path.join(root, "assembly", "big_assembly"),
            type

        )
        job.name = "ray_big_assembly"
        job.subname = "ray"
        jobs.append(job)
        
        return jobs

    def gene_prediction(self):
        root = self._root_dir
        jobs = []
        
        for readset in self.readsets:
            infile = os.path.join(root, "assembly", readset.sample.name, "Contigs.fasta")
            outdir = root + "gene_prediction/" + readset.sample.name
        
            if not os.path.exists(outdir):
                os.makedirs(os.path.join(outdir))
                    
            if(self._big_assembly_only == False):
                job = shotgun_metagenomics.metagenemark(
                    infile,
                    outdir + "/" + readset.sample.name + ".gff",
                    outdir + "/" + readset.sample.name + ".fna",
                    outdir + "/" + readset.sample.name + ".faa",
                    outdir + "/" + readset.sample.name + "_renamed.fna"
                )
                job.name = "metagenemark_" + readset.sample.name
                job.subname = "metagenemark"
                jobs.append(job)
        
        # Gene prediction for big assembly.
        outdir = os.path.join(root, "gene_prediction", "big_assembly")
        
        if not os.path.exists(outdir):
            os.makedirs(os.path.join(outdir))

        job = shotgun_metagenomics.metagenemark(
            os.path.join(root, "assembly", "big_assembly", "Contigs.fasta"),
            outdir + "/" + "big_assembly.gff",
            outdir + "/" + "big_assembly.fna",
            outdir + "/" + "big_assembly.faa",
            outdir + "/" + "big_assembly_renamed.fna"
        )
        job.name = "metagenemark_big_assembly"
        job.subname = "metagenemark"
        jobs.append(job)
        
        return jobs

     
    def exonerate(self):
        jobs = []

        for readset in self.readsets:
            infile_fna = os.path.join(self._root_dir, "gene_prediction", readset.sample.name, readset.sample.name + ".fna")
            infile_faa = os.path.join(self._root_dir, "gene_prediction", readset.sample.name, readset.sample.name + ".faa")
            chunks_dir = os.path.join(self._root_dir, "gene_prediction", readset.sample.name, "fasta_chunks")
            num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
            
            if(self._big_assembly_only == False):
                # FNA
                job = shotgun_metagenomics.exonerate(
                    infile_fna,
                    chunks_dir,
                    num_chunks,
                    readset.sample.name + "_renamed.fna" 
                )
                job.name = "exonerate_" + readset.sample.name
                job.subname = "exonerate"
                jobs.append(job)
                
                # FAA
                job = shotgun_metagenomics.exonerate(
                    infile_faa,
                    chunks_dir,
                    num_chunks,
                    readset.sample.name + "_renamed.faa" 
                )
                job.name = "exonerate_" + readset.sample.name
                job.subname = "exonerate"
                jobs.append(job)

        # Split for big assembly
        infile_fna = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "big_assembly.fna")
        infile_faa = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "big_assembly.faa")
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "fasta_chunks")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
        
        #FNA
        job = shotgun_metagenomics.exonerate(
            infile_fna,
            chunks_dir,
            num_chunks,
            "big_assembly.fna" 
        )
        job.name = "exonerate_big_assembly_fna"
        job.subname = "exonerate"
        jobs.append(job)
        
        # FAA
        job = shotgun_metagenomics.exonerate(
            infile_faa,
            chunks_dir,
            num_chunks,
            "big_assembly.faa" 
        )
        job.name = "exonerate_big_assembly_faa"
        job.subname = "exonerate"
        jobs.append(job)
 
        return jobs    

    def blastn_nt(self):
        jobs = []
        
        if(self._big_assembly_only == False):
            for readset in self.readsets:
                chunks_dir = os.path.join(self._root_dir, "gene_prediction", readset.sample.name, "fasta_chunks")
                blast_dir = os.path.join(self._root_dir, "gene_annotation", readset.sample.name, "blastn_nt")
                num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')

                for i in range(num_chunks):
                    job = shotgun_metagenomics.blastn(
                        os.path.join(chunks_dir, readset.sample.name + "_renamed.fna_chunk_{:07d}".format(i)),
                        os.path.join(blast_dir, "blastn_chunk_{:07d}.tsv".format(i)),
                        blast_dir 
                    )
                    job.name = "blastn_nt_" + readset.sample.name
                    job.subname = "blastn"
                    jobs.append(job)
      
            # Merge output chunks
          
        # Do blastn on nt for big assembly  
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastn_nt")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')

        for i in range(num_chunks):
            job = shotgun_metagenomics.blastn(
                os.path.join(chunks_dir, "big_assembly.fna_chunk_{:07d}".format(i)),
                os.path.join(blast_dir, "blastn_chunk_{:07d}.tsv".format(i)),
                blast_dir 
            )
            job.name = "blastn_nt_big_assembly_merge"
            job.subname = "blastn"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastn_nt.tsv"),
            num_chunks,
            "blastn" 
        )
        job.name = "blastn_nt_big_assembly"
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
        if(self._big_assembly_only == False):
            for readset in self.readsets:
                chunks_dir = os.path.join(self._root_dir, "gene_prediction", readset.sample.name, "fasta_chunks")
                blast_dir = os.path.join(self._root_dir, "gene_annotation", readset.sample.name, "blastx_kegg")
                num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
    
                for i in range(num_chunks):
                    job = shotgun_metagenomics.blastx(
                        os.path.join(chunks_dir, readset.sample.name + "_faa_chunk_{:07d}".format(i)),
                        os.path.join(blast_dir, "blastp_chunk_{:07d}.tsv".format(i)),
                        blast_dir,
                        config.param('kegg', 'db')
                    )
                    job.name = "blastp_kegg"
                    job.subname = "blastp"
                    jobs.append(job)

        # Do blastp on KEGG for big assembly  
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastp_kegg")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')

        for i in range(num_chunks):
            job = shotgun_metagenomics.blastp(
                os.path.join(chunks_dir, "big_assembly.faa_chunk_{:07d}".format(i)),
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
                os.path.join(chunks_dir, "big_assembly.faa_chunk_{:07d}".format(i)),
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
                os.path.join(chunks_dir, "big_assembly.faa_chunk_{:07d}".format(i)),
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
                os.path.join(chunks_dir, "big_assembly.faa_chunk_{:07d}".format(i)),
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
                os.path.join(chunks_dir, "big_assembly.faa_chunk_{:07d}".format(i)),
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

    ## Gene abundance: Mapping reads back on contigs.
    def gene_abundance(self):
        jobs = []
        cov_list = []
        
        # Will make index for bwa. and also bed file for computing reads spanning later.
        reference = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "big_assembly.fna")
        job = shotgun_metagenomics.make_index(
            reference
        )
        job.name = "make_index"
        job.subname = "make_index"
        jobs.append(job)
            
        bed = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "big_assembly.bed")

        for readset in self.readsets:
            bam = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".bam")
            cov = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".cov")
            cov_list.append(cov)
    
            outdir = os.path.join(self._root_dir, "gene_abundance/", readset.sample.name)
            if not os.path.exists(outdir):
                os.makedirs(os.path.join(outdir))
       
            
            if( readset.run_type == "PAIRED_END"):
                infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
                out1 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
                out2 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
                
                job = shotgun_metagenomics.split_pairs(
                    infile,
                    out1,
                    out2
                )
                job.name = "split_pairs_" + readset.sample.name
                job.subname = "split_pairs"
                jobs.append(job)

                job = shotgun_metagenomics.bwa_mem(
                    reference,
                    out1,
                    out2,
                    bam
                )
                job.name = "bwa_mem-" + readset.sample.name
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
                bam,
                bed,
                cov
            )
            job.name = "bedtoolsCov-" + readset.sample.name
            job.subname = "bedtools"
            jobs.append(job)
            
        # Once all coverage has been computed, merge all tables.
        #sys.stderr.write('gene abundance: ' + ','.join([str(x) for x in cov_list] ) + '\n')
        job = shotgun_metagenomics.merge_counts(
            cov_list,
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv")
        )
        job.name = "merge_gene_abundance"
        job.subname = "merge_gene_abundance"
        jobs.append(job)
         
        return jobs 
      
    def statistics(self):
        jobs = []

        # For each analysis generate heatmap + pca plot
        job = shotgun_metagenomics.reads_abundance_analysis(
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "statistics", "reads_abundance")
        )
        job.name = "gene_abundance_analysis"
        jobs.append(job)

        # For each analysis generate heatmap + pca plot + stack barplots?
        job = shotgun_metagenomics.kegg_analysis(
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "statistics", "kegg"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastp_kegg_parsed.tsv")
        )
        job.name = "kegg_analysis"
        jobs.append(job)
        
        job = shotgun_metagenomics.cog_analysis(
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "statistics", "cog"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "rpsblast_cog.tsv") 
        )
        job.name = "cog_analysis"
        jobs.append(job)
        
        job = shotgun_metagenomics.kog_analysis(
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "statistics", "kog"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "rpsblast_kog.tsv") 
        )
        job.name = "kog_analysis"
        jobs.append(job)
        
        job = shotgun_metagenomics.pfam_analysis(
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "statistics", "pfam"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "hmmscan_pfam_tblout.tsv") 
        )
        job.name = "pfam_analysis"
        jobs.append(job)
        
        job = shotgun_metagenomics.tigrfam_analysis(
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "statistics", "tigrfam"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "hmmscan_tigrfam_tblout.tsv") 
        )
        job.name = "tigrfam_analysis"
        jobs.append(job)
        
        return jobs
    
    def generate_report(self):
        jobs = []
        #job = shotgun_metagenomics.mymethod(
        #)
        #job.name = "myjobname"
        #jobs.append(job)
        sys.stderr.write('not implemented yet\n')
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
            self.trim,#1
            self.duk,#2
            self.assembly,#3
            self.gene_prediction,#4
            self.exonerate,#5
            self.blastn_nt,#6
            self.blastp_kegg,#7
            self.rpsblast_cog,#8
            self.rpsblast_kog,#9
            self.hmmscan_tigrfam,#10
            self.hmmscan_pfam,#11
            self.gene_abundance,#12
            self.generate_report,#13
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
        self._parser_local.add_argument("-b", "--big-assembly-only", help="Will generate one big assembly only (default: false)", action="store_true")

        # barcodes
        self._args_local = self._parser_local.parse_args()
        config.parse_files(self._args_local.config)
        self._config = self._args_local.config[0]
        self._config = os.path.abspath(self._config.name)
        
        self._root_dir = self._args_local.output_dir
        if not os.path.exists(self._root_dir):
            os.makedirs(self._root_dir)
        
        # Big assembly only?
        self._big_assembly_only = self._args_local.big_assembly_only
        sys.stderr.write("self._big_assembly_only: " + str(self._big_assembly_only) + "\n")

        # Make directories
        self.make_directories(self._root_dir)
  
    # Define and make directories. Also desing initial infile.
    def make_directories(self, root_dir):
         
        if not os.path.exists(root_dir):
            os.makedirs(root_dir)
        if not os.path.exists(root_dir + "/qced_reads"):
            os.makedirs(os.path.join(root_dir, "qced_reads"))
        if not os.path.exists(root_dir + "/assembly"):
            os.makedirs(os.path.join(root_dir, "assembly"))
        if not os.path.exists(root_dir + "/gene_prediction"):
            os.makedirs(os.path.join(root_dir, "gene_prediction"))
        #if not os.path.exists(root_dir + "/gene_clustering"):
        #    os.makedirs(os.path.join(root_dir, "gene_clustering"))
        if not os.path.exists(root_dir + "/gene_abundance"):
            os.makedirs(os.path.join(root_dir, "gene_abundance"))
        if not os.path.exists(root_dir + "/taxonomic_annotation"):
            os.makedirs(os.path.join(root_dir, "taxonomic_annotation"))
        if not os.path.exists(root_dir + "/gene_annotation"):
            os.makedirs(os.path.join(root_dir, "gene_annotation"))
        #if not os.path.exists(root_dir + "/gene_abundance_clustering"):
        #    os.makedirs(os.path.join(root_dir, "gene_abundance_clustering"))
        #if not os.path.exists(root_dir + "/mgs_augmented_assembly"):
        #    os.makedirs(os.path.join(root_dir, "mgs_augmented_assembly"))
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
        self.argparser.add_argument("-b", "--big-assembly-only", help="Will generate one big assembly only (default: false)", action="store_true")
        sys.stderr.write('Running quick assembly pipeline\n')
        super(Quick_assembly, self).__init__()
                
Quick_assembly().submit_jobs()
