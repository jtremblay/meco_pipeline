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
from bio import shotgun_metatranscriptomics
from bio import microbial_ecology

from pipelines.illumina import illumina

# Global scope variables
log = logging.getLogger(__name__)

class Metatranscriptomics(illumina.Illumina):

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
        
            # Merge R1 and R2 to get an interleaved file. We do this here because of memory requirements
            # for downstream duk steps. Here we assume that R1 and R2 are in the same exact order.
            job = shotgun_metagenomics.create_interleaved_fastq(
                trim_file_prefix + "pair1.fastq.gz",
                trim_file_prefix + "pair2.fastq.gz",
                trim_file_prefix + "interleaved.fastq",
                trim_file_prefix + "interleaved.fastq.gz"
            )
            job.name = "create_interleaved_fastq_" + readset.sample.name
            job.subname = "interleaved_fastq"
            jobs.append(job) 
        
        #job = shotgun_metagenomics.merge_trimstats(
        #    trim_stats,
        #    os.path.join(outdir, "qced_reads", "trim_stats.tsv")
        #)
        #job.name = "merge_trim_stats"
        #job.subname = "merge_trim_stats"
        #jobs.append(job) 

        return jobs
            
    def duk(self):
        jobs=[]
        outdir = self._root_dir
        logs_R1 = []

        for readset in self.readsets:
            trim_file_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".trim.")
            outfile_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".")
            out1up = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_unpaired_R1.fastq.gz")
            out2up = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_unpaired_R2.fastq.gz")
            out1p = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
            out2p = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
            log = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".duk_contam_interleaved_log.txt")

            job = shotgun_metagenomics.duk_gz(
                trim_file_prefix + "interleaved.fastq.gz",
                outfile_prefix + "contam.fastq.gz",
                outfile_prefix + "ncontam.fastq.gz",
                log,
                config.param('DB', 'contaminants', 1, 'filepath')
            )
            job.name = "duk_interleaved_" + readset.sample.name
            job.subname = "duk"
            jobs.append(job)
        

            job = shotgun_metagenomics.remove_unpaired_reads_and_split(
                outfile_prefix + "ncontam.fastq.gz",
                out1up,
                out2up,
                out1p,
                out2p,
                outfile_prefix + "ncontam_paired.fastq.gz"
            )
            job.name = "remove_unpaired_and_split_" + readset.sample.name
            job.subname = "remove_unpaired"
            jobs.append(job) 

            logs_R1.append(log)

        # Compile duk logs. TODO
        job = shotgun_metagenomics.merge_duk_logs_interleaved(
            logs_R1,
            os.path.join(self._root_dir, "qced_reads", "duk_merged_logs.tsv")
        )
        job.name = "merge_duk_logs"
        job.subname = "merge_duk_logs"
        jobs.append(job)

        return jobs
     
    def abundance(self):
        jobs = []
        cov_list_contigs = []
        cov_list_genes = []
        flagstats_contigs_list = []
        flagstats_genes_list = []
        trimmomatic_list = []
        
        reference_contigs = os.path.join(self._root_dir, "assembly", "Contigs.fasta")
        bed_contigs = os.path.join(self._root_dir, "assembly", "Contigs.bed")
        bwt_contigs = os.path.join(self._root_dir, "assembly", "Contigs.fasta.bwt")
        bed_genes = os.path.join(self._root_dir, "assembly", "Contigs_genes.bed")
        
        # Will make index for bwa. and also bed file for computing reads spanning later.
        #reference_genes = os.path.join(self._root_dir, "assembly", "Genes.fasta")
        #bwt_contigs = os.path.join(self._root_dir, "assembly", "Contigs.fasta.bwt")
        #bwt_genes = os.path.join(self._root_dir, "assembly", "Genes.fasta.bwt")

        #job = shotgun_metagenomics.make_index(
        #    reference_genes,
        #    bwt_genes
        #)
        #job.name = "make_index_genes"
        #job.subname = "make_index"
        #jobs.append(job)
         
        #bwt_genes = os.path.join(self._root_dir, "assembly", "Genes.fasta.bwt")

        for readset in self.readsets:
            # Trimmomatic reads.
            trimmomatic = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".trim.stats.csv")
            trimmomatic_list.append(trimmomatic)
            # BAM contigs
            bam_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".bam")
            flagstats_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".flagstats")
            flagstats_contigs_list.append(flagstats_contigs)
            # BAM genes
            #bam_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".bam")
            #flagstats_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".flagstats")
            #flagstats_genes_list.append(flagstats_genes)
            # Cov contigs
            cov_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".cov")
            cov_list_contigs.append(cov_contigs)
            # Cov genes
            cov_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".cov")
            cov_list_genes.append(cov_genes)
            # Outdir
            outdir_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name)
            outdir_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name)

            if not os.path.exists(outdir_genes):
                os.makedirs(os.path.join(outdir_genes))
            if not os.path.exists(outdir_contigs):
                os.makedirs(os.path.join(outdir_contigs))
       
            if( readset.run_type == "PAIRED_END"):
                infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
                infile_R1 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
                infile_R2 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
                
                # map against contigs
                job = shotgun_metagenomics.bwa_mem_samtools(
                    reference_contigs,
                    bwt_contigs,
                    infile_R1,
                    infile_R2,
                    bam_contigs
                )
                job.name = "bwa_mem-contigs" + readset.sample.name
                job.subname = "bwa"
                jobs.append(job)
                
                job = shotgun_metagenomics.flagstats(
                    bam_contigs,
                    flagstats_contigs
                )
                job.name = "flagstats-contigs-" + readset.sample.name
                job.subname = "flagstats"
                jobs.append(job)
               
                # map against genes
                #job = shotgun_metagenomics.bwa_mem_samtools(
                #    reference_genes,
                #    bwt_genes,
                #    infile_R1,
                #    infile_R2,
                #    bam_genes
                #)
                #job.name = "bwa_mem-genes" + readset.sample.name
                #job.subname = "bwa"
                #jobs.append(job)
                
                #job = shotgun_metagenomics.flagstats(
                #    bam_genes,
                #    flagstats_genes
                #)
                #job.name = "flagstats-genes-" + readset.sample.name
                #job.subname = "flagstats"
                #jobs.append(job)
                
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

        # merge alignment rates for contigs
        job = shotgun_metagenomics.merge_flagstats(
            flagstats_contigs_list,
            trimmomatic_list,
            os.path.join(self._root_dir, "contigs_abundance", "qc_mapping_stats_contigs.tsv")
        )
        job.name = "flagstats_merge_contigs"
        job.subname = "flagstats"
        jobs.append(job)
       
        # merge alignment rates for genes
        #job = shotgun_metagenomics.merge_flagstats(
        #    flagstats_genes_list,
        #    trimmomatic_list,
        #    os.path.join(self._root_dir, "abundance", "qc_mapping_stats_genes.tsv")
        #)
        #job.name = "flagstats_merge_genes"
        #job.subname = "flagstats"
        #jobs.append(job)
         
        return jobs 
    
    def DEG(self):
        
        jobs = []
        cov_list_contigs = []
        cov_list_genes = []
        flagstats_contigs_list = []
        trimmomatic_list = []
        #htseq_list_contigs = []
        
        # Will make index for bwa. and also bed file for computing reads spanning later.
        reference_contigs = os.path.join(self._root_dir, "assembly", "Contigs.fasta")
        reference_genes = os.path.join(self._root_dir, "assembly", "Genes.fasta")
        bed_contigs = os.path.join(self._root_dir, "assembly", "Contigs.bed")
        bwt_contigs = os.path.join(self._root_dir, "assembly", "Contigs.fasta.bwt")
        bed_genes = os.path.join(self._root_dir, "assembly", "Genes.bed")
        bwt_genes = os.path.join(self._root_dir, "assembly", "Genes.fasta.bwt")

        # Generate counts from alignments generated in step <abundance>. 
        for readset in self.readsets:
            bam_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".bam")
            bam_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".bam")
            #htseq_contigs = os.path.join(self._root_dir, "DEG", readset.sample.name, readset.name + ".htseq")
            #htseq_list_contigs.append(htseq_contigs)
            cov_genes = os.path.join(self._root_dir, "DEG", readset.sample.name, readset.name + ".genes.cov")
            cov_list_genes.append(cov_genes)
        
            outdir = os.path.join(self._root_dir, "DEG", readset.sample.name)
            if not os.path.exists(outdir):
                os.makedirs(os.path.join(outdir))
         
            # raw counts with bedtools. Bedtools used for genes abundance only
            job = shotgun_metagenomics.coverage_bed(
                bam_genes,
                bed_genes,
                cov_genes
            )
            job.name = "bedtoolsCov-genes-" + readset.sample.name
            job.subname = "bedtools"
            jobs.append(job)
            
            ## raw counts with htseq. For Htseq, use contigs abundance AND gene models (GFF).
            #job = shotgun_metatranscriptomics.htseq_count(
            #    bam_contigs,
            #    config.param('DEFAULT', 'gff', 1, 'filepath'),
            #    htseq_contigs
            #)
            #job.name = "htseq-contigs" + readset.sample.name
            #job.subname = "htseq"
            #jobs.append(job)
        
        # Count reads with bedtools and htseq
        job = shotgun_metatranscriptomics.merge_bedtools_counts(
            cov_list_genes,
            os.path.join(self._root_dir, "DEG", "merged_abundance_bedtools.tsv"),
            os.path.join(self._root_dir, "DEG", "merged_abundance_RPKM_bedtools.tsv")
        )
        job.name = "merge_gene_abundance_genes"
        job.subname = "merge_gene_abundance"
        jobs.append(job)
        
        #job = shotgun_metatranscriptomics.merge_htseq_counts(
        #    htseq_list_contigs,
        #    cov_list_genes,
        #    os.path.join(self._root_dir, "DEG", "merged_abundance_htseq.tsv"),
        #    os.path.join(self._root_dir, "DEG", "merged_abundance_RPKM_htseq.tsv")
        #)
        #job.name = "merge_gene_abundance_contigs"
        #job.subname = "merge_gene_abundance"
        #jobs.append(job)

        ## edgeR for DEG - htseq
        #job = shotgun_metatranscriptomics.edger_single(
        #    os.path.join(self._root_dir, "DEG", "merged_abundance_htseq.tsv"),
        #    os.path.join(self._root_dir, "DEG", "DEG_htseq_normalized.tsv")
        #)
        #job.name = "DEG_single_htseq"
        #job.subname = "DEG_single"
        #jobs.append(job)
        
        # edgeR for DEG - bedtools
        job = shotgun_metatranscriptomics.edger_single(
            os.path.join(self._root_dir, "DEG", "merged_abundance_bedtools.tsv"),
            os.path.join(self._root_dir, "DEG_bedtools", "DEG_bedtools_normalized.tsv")
        )
        job.name = "DEG_single_bedtools"
        job.subname = "DEG_single"
        jobs.append(job)
        
        design_files = config.param('DEFAULT', 'designFile', 1, 'string')
        designs = design_files.split(":")
        for design_file in designs:
            #job = shotgun_metatranscriptomics.edger(
            #    os.path.join(self._root_dir, "DEG", "merged_abundance_htseq.tsv"),
            #    os.path.join(self._root_dir, design_file),
            #    os.path.join(self._root_dir, "DEG")
            #)
            #job.name = "DEG"
            #job.subname = "DEG"
            #jobs.append(job)
            
            job = shotgun_metatranscriptomics.edger(
                os.path.join(self._root_dir, "DEG", "merged_abundance_bedtools.tsv"),
                os.path.join(self._root_dir, design_file),
                os.path.join(self._root_dir, "DEG_bedtools")
            )
            job.name = "DEG_bedtools"
            job.subname = "DEG"
            jobs.append(job)
        
        #job = shotgun_metagenomics.draw_pcoa(
        #    os.path.join(self._root_dir, "DEG"),
        #    config.param('DEFAULT', 'mappingFile'),
        #    "_normalized.tsv",
        #    os.path.join(self._root_dir, "DEG", "DEG.done"),
        #    os.path.join(self._root_dir, "DEG", "DEG_htseq_normalized.tsv")
        #)
        #job.name = "drawPcoa_htseq"
        #job.subname = "drawPcoa"
        #jobs.append(job)
        
        job = shotgun_metagenomics.draw_pcoa(
            os.path.join(self._root_dir, "DEG_bedtools"),
            config.param('DEFAULT', 'mappingFile'),
            "_normalized.tsv",
            os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
            os.path.join(self._root_dir, "DEG_bedtools", "DEG_bedtools_normalized.tsv")
        )
        job.name = "drawPcoa_bedtools"
        job.subname = "drawPcoa"
        jobs.append(job)
        
        job = shotgun_metagenomics.edger_glm(
            os.path.join(self._root_dir, "DEG", "merged_abundance_bedtools.tsv"),
            os.path.join(self._root_dir, "DEG_GLM")
        )
        job.name = "DEG"
        job.subname = "DEG_GLM"
        jobs.append(job)

        return jobs 
     
    def taxonomy(self):
        jobs = []
        
        # We already have taxonomy from metagenomics (assembly) pipelne.
        job = shotgun_metagenomics.generate_otu_table(
            os.path.join(self._root_dir, "assembly", "taxonomy.tsv"),
            os.path.join(self._root_dir, "DEG", "merged_abundance_bedtools.tsv"),
            os.path.join(self._root_dir, "taxonomy", "otu_table.tsv")
        )
        job.name = "generate_otu_table"
        job.subname = "generate_otu_table"
        jobs.append(job)
        
        job = microbial_ecology.convert_otu_to_biom_hdf5(
            os.path.join(self._root_dir, "taxonomy", "otu_table.tsv"),
            os.path.join(self._root_dir, "taxonomy", "otu_table.biom")
        )
        job.name = "convert_otu_to_biom"
        job.subname = "convert_otu_table"
        jobs.append(job)

        job = microbial_ecology.split_otu_table(
            os.path.join(self._root_dir, "taxonomy", "otu_table.tsv"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea.tsv"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_others.tsv"),
            "bactArch"
        )
        job.name = "split_otu_table_bactArch"
        job.name = "split_otu_table"
        jobs.append(job)

        job = microbial_ecology.convert_otu_to_biom_hdf5(
            os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea.tsv"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea.biom")
        )
        job.name = "convert_otu_table_to_biom_bactArch"
        job.subname = "convert_otu_table"
        jobs.append(job)
        
        job = microbial_ecology.convert_otu_to_biom_hdf5(
            os.path.join(self._root_dir, "taxonomy", "otu_table_others.tsv"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_others.biom")
        )
        job.name = "convert_otu_table_to_biom_others"
        job.subname = "convert_otu_table"
        jobs.append(job)

        job = microbial_ecology.rarefy_hdf5(
            os.path.join(self._root_dir, "taxonomy", "otu_table.biom"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_rarefied.biom"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_rarefied.tsv"),
            config.param('rarefaction', 'n_genes', 1, 'posint')
        )    
        job.name = "rarefy_otu_table_all_organisms"
        job.subname = "rarefaction"
        jobs.append(job)
        
        job = microbial_ecology.normalize_otu_table_hdf5(
            os.path.join(self._root_dir, "taxonomy", "otu_table.tsv"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_normalized.biom"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_normalized.tsv")
        )    
        job.name = "normalize_otu_table_all_organisms"
        job.subname = "normalization"
        jobs.append(job)
        
        job = microbial_ecology.rarefy_hdf5(
            os.path.join(self._root_dir, "taxonomy", "otu_table_others.biom"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_others_rarefied.biom"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_others_rarefied.tsv"),
            config.param('rarefaction', 'n_genes', 1, 'posint')
        )    
        job.name = "rarefy_otu_table_others"
        job.subname = "rarefaction"
        jobs.append(job)
        
        job = microbial_ecology.normalize_otu_table_hdf5(
            os.path.join(self._root_dir, "taxonomy", "otu_table_others.tsv"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_others_normalized.biom"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_others_normalized.tsv")
        )    
        job.name = "normalize_otu_table_others"
        job.subname = "normalization"
        jobs.append(job)
        
        job = microbial_ecology.rarefy_hdf5(
            os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea.biom"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea_rarefied.biom"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea_rarefied.tsv"),
            config.param('rarefaction', 'n_genes', 1, 'posint')
        )    
        job.name = "rarefy_otu_table_bactArchaea"
        job.subname = "rarefaction"
        jobs.append(job)
        
        job = microbial_ecology.normalize_otu_table_hdf5(
            os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea.tsv"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea_normalized.biom"),
            os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea_normalized.tsv")
        )    
        job.name = "normalize_otu_table_bactArchaea"
        job.subname = "normalization"
        jobs.append(job)
        
        # ALL
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_absolute(
                os.path.join(self._root_dir, "taxonomy", "otu_table_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "taxonomy", "all", "absolute")
            )
            job.name = "tax_summary_absolute_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
        
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_relative(
                os.path.join(self._root_dir, "taxonomy", "otu_table_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "taxonomy", "all", "relative")
            )
            job.name = "tax_summary_relative_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
        
        # OTHERS
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_absolute(
                os.path.join(self._root_dir, "taxonomy", "otu_table_others_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "taxonomy", "others", "absolute")
            )
            job.name = "tax_summary_absolute_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
        
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_relative(
                os.path.join(self._root_dir, "taxonomy", "otu_table_others_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "taxonomy", "others", "relative")
            )
            job.name = "tax_summary_relative_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
        
        # BACTERIA-ARCHAEA
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_absolute(
                os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "taxonomy", "bacteriaArchaea", "absolute")
            )
            job.name = "tax_summary_absolute_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
        
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_relative(
                os.path.join(self._root_dir, "taxonomy", "otu_table_bacteriaArchaea_rarefied.biom"),
                i,
                os.path.join(self._root_dir, "taxonomy", "bacteriaArchaea", "relative")
            )
            job.name = "tax_summary_relative_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
            
        # Plot taxa - ALL
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "taxonomy", "all", "relative", "otu_table_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "taxonomy", "all", "relative", "plots"),
                "taxonomy_L" + str(i) + "_relative"
            )
            job.name = "plot_taxa_single_relative_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
            
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "taxonomy", "all", "absolute", "otu_table_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "taxonomy", "all", "absolute", "plots"),
                "taxonomy_L" + str(i) + "_absolute"
            )
            job.name = "plot_taxa_single_absolute_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
        
        # OTHERS
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "taxonomy", "others", "relative", "otu_table_others_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "taxonomy", "others", "relative", "plots"),
                "taxonomy_L" + str(i) + "_relative"
            )
            job.name = "plot_taxa_single_relative_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
            
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "taxonomy", "others", "absolute", "otu_table_others_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "taxonomy", "others", "absolute", "plots"),
                "taxonomy_L" + str(i) + "_absolute"
            )
            job.name = "plot_taxa_single_absolute_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
        
        # Bacteria-Archaea
        for i in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "taxonomy", "bacteriaArchaea", "relative", "otu_table_bacteriaArchaea_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "taxonomy", "bacteriaArchaea", "relative", "plots"),
                "taxonomy_L" + str(i) + "_relative"
            )
            job.name = "plot_taxa_single_relative_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
            
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._root_dir, "taxonomy", "bacteriaArchaea", "absolute", "otu_table_bacteriaArchaea_rarefied" + "_L" + str(i) + ".txt"),
                os.path.join(self._root_dir, "taxonomy", "bacteriaArchaea", "absolute", "plots"),
                "taxonomy_L" + str(i) + "_absolute"
            )
            job.name = "plot_taxa_single_absolute_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)

        return jobs
    
    def overrepresentation (self):
        jobs = []

        types = ["pathways", "modules", "K"]
        for type in types:
            job = shotgun_metagenomics.kegg_overrep(
                os.path.join(self._root_dir, "DEG_bedtools"),
                os.path.join(self._root_dir, "gene_annotation", "blastp_kegg_parsed.tsv"),
                type,
                os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
                "DEG"
            )
            job.name = "getKegg_" + type
            job.subname = "getKegg"
            jobs.append(job)
         
        # COG and KOG
        job = shotgun_metagenomics.cog_overrep(
            os.path.join(self._root_dir, "DEG_bedtools"),
            os.path.join(self._root_dir, "gene_annotation", "rpsblast_cog.tsv"),
            os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
            "DEG"
        )
        job.name = "getCOG"
        job.subname = "getCOG"
        jobs.append(job)
        
        job = shotgun_metagenomics.kog_overrep(
            os.path.join(self._root_dir, "DEG_bedtools"),
            os.path.join(self._root_dir, "gene_annotation", "rpsblast_kog.tsv"),
            os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
            "DEG"
        )
        job.name = "getKOG"
        job.subname = "getKOG"
        jobs.append(job)
        
        # Pfam and tigrfam
#        job = shotgun_metagenomics.pfam_overrep(
#            os.path.join(self._root_dir, "DEG_bedtools"),
#            os.path.join(self._root_dir, "gene_annotation", "hmmscan_pfam_tblout.tsv"),
#            os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
#            "DEG"
#        )
#        job.name = "getPfam"
#        job.subname = "getPfam"
#        jobs.append(job)
#         
#        job = shotgun_metagenomics.tigrfam_overrep(
#            os.path.join(self._root_dir, "DEG_bedtools"),
#            os.path.join(self._root_dir, "gene_annotation", "hmmscan_tigrfam_tblout.tsv"),
#            os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
#            "DEG"
#        )
#        job.name = "getTigrfam"
#        job.subname = "getTigrfam"
#        jobs.append(job)
#        
#        # heatmaps
#        prefixes = ["DEG_normalized_significant_keggK.tsv", "DEG_normalized_significant_keggModules.tsv", "DEG_normalized_significant_keggPathways.tsv", "DEG_normalized_significant_COG.tsv", "DEG_normalized_significant_KOG.tsv", "DEG_normalized_significant_PFAM.tsv", "DEG_normalized_significant_TIGRFAM.tsv"]
#        for prefix in prefixes:
#            job = shotgun_metagenomics.draw_heatmap(
#                os.path.join(self._root_dir, "DEG_bedtools"),
#                config.param('DEFAULT', 'mappingFile'),
#                prefix,
#                os.path.join(self._root_dir, "DEG_bedtools", "DEG.done")
#            )
#            job.name = "draw_heatmap"
#            job.subname = "draw_heatmap_" +  prefix
#            jobs.append(job)
         
        return jobs

    # Here generate final GFF (for viewing data in a genome browser). 
    # And generate final DEG sheets. with logFC, gene_name and actual normalized values.
    def finalize(self):
        jobs = []
        
        DEGs = ["DEG","DEG_bedtools"]

        for DEG in DEGs:
            job = shotgun_metagenomics.annotate_DDA(
                os.path.join(self._root_dir, DEG),
                config.param('DEFAULT', 'annotations', 1, 'filepath')
            )
            job.name = "annotate_" + DEG
            job.subname = "annotate_DEG"
            jobs.append(job)

        return jobs

    def generate_report(self):
        jobs = []
        
        i = datetime.datetime.now()
        year = "%04d" %i.year
        day = "%02d" %i.day
        month = "%02d" %i.month
        title = config.param('report', 'report.title')
        date = year + "_" + month + "_" + day + "_metatranscriptomics_" + title
        
        job = shotgun_metagenomics.concatenate_reports(
            self._config,
            self._root_dir,
            "Metatranscriptomics",
            self._root_dir + "/" + date
        )
        job.name = "deliverables"
        job.subname = "report"
        jobs.append(job)
        return jobs
    
    def exploratory(self):
        jobs = []
        
        i = datetime.datetime.now()
        year = "%04d" %i.year
        day = "%02d" %i.day
        month = "%02d" %i.month
        title = config.param('report', 'report.title')
        date = year + "_" + month + "_" + day

        log_fc_string = config.param('exploratory', 'cutoff_logFC', 1, 'string')
        log_fcs = log_fc_string.split(":")
                
        gene_list_files = config.param('exploratory', 'gene_list', 1, 'string')
        gene_lists = gene_list_files.split(":")
        
        for gene_list in gene_lists:
            gene_list_prefix = os.path.splitext(gene_list)[0]
            gene_list_prefix = gene_list_prefix.replace("./","")
            for log_fc in log_fcs:
                if not os.path.exists( os.path.join(self._root_dir, "DEG_bedtools", "exploratory", "log_fc-"+log_fc) ):
                    os.makedirs( os.path.join(self._root_dir, "DEG_bedtools", "exploratory", "log_fc-"+log_fc) )
                
                job = shotgun_metatranscriptomics.explore_heatmaps(
                    log_fc,
                    os.path.join(self._root_dir, "DEG_bedtools", "DEG_normalized_significant.tsv"),
                    os.path.abspath(os.path.join(self._root_dir, "DEG_bedtools")),
                    os.path.abspath(os.path.join(self._root_dir, "DEG_bedtools", "exploratory", "log_fc-" + log_fc)),
                    date,
                    gene_list
                )
                job.name = "exploratory_" + gene_list_prefix + "_" + log_fc
                job.subname = "exploratory"
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
            self.trim,#1
            self.duk,#2
            self.abundance,#3
            self.DEG,#4
            self.overrepresentation,#5
            self.finalize,#6
            self.taxonomy,#7
            self.generate_report,#8
            self.exploratory,#9
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
        if not os.path.exists(root_dir + "/DEG"):
            os.makedirs(os.path.join(root_dir, "DEG"))
        if not os.path.exists(root_dir + "/qced_reads"):
            os.makedirs(os.path.join(root_dir, "qced_reads"))
        if not os.path.exists(root_dir + "/abundance"):
            os.makedirs(os.path.join(root_dir, "abundance"))
        if not os.path.exists(root_dir + "/taxonomy"):
            os.makedirs(os.path.join(root_dir, "taxonomy"))
        if not os.path.exists(root_dir + "/taxonomy/all/absolute/plots"):
            os.makedirs(os.path.join(root_dir, "taxonomy/all/absolute/plots"))
            os.makedirs(os.path.join(root_dir, "taxonomy/all/relative/plots"))
        if not os.path.exists(root_dir + "/taxonomy/others/absolute/plots"):
            os.makedirs(os.path.join(root_dir, "taxonomy/others/absolute/plots"))
            os.makedirs(os.path.join(root_dir, "taxonomy/others/relative/plots"))
        if not os.path.exists(root_dir + "/taxonomy/bacteriaArchaea/absolute/plots"):
            os.makedirs(os.path.join(root_dir, "taxonomy/bacteriaArchaea/absolute/plots"))
            os.makedirs(os.path.join(root_dir, "taxonomy/bacteriaArchaea/relative/plots"))
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
        sys.stderr.write('Metatranscriptomics pipeline\n')
        super(Metatranscriptomics, self).__init__()
                
Metatranscriptomics().submit_jobs()
