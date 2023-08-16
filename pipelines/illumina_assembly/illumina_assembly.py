#!/usr/bin/env python

#######################################################################################################
# Copyright (C) 2023 Julien Tremblay
#
# This file is part of MECO Pipelines.
#
# MECO Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MECO Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MECO Pipelines.  If not, see <http://www.gnu.org/licenses/>.
#######################################################################################################

# Python Standard Modules
import logging
import os
import sys
import time
import argparse
# Append meco_pipelines directory to Python library path
#sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MECO Modules
from core.config import *
from core.job import *
from bio.readset import *

from bio import genome_assembly
from bio import shotgun_metagenomics
from bio import microbial_ecology
from pipelines import common

log = logging.getLogger(__name__)

class GenomeAssembly(common.MECOPipeline):
    """
    Illumina Genome Assembly Pipeline
    Written by Julien Tremblay
    ========================

    Simple pipeline that will perform QC (Trimmomatic) and de novo assembly using either metaSPAdes or Megahit - and 
    probably more choices to come in the future. Written for bacterial genomes assembly.

    """
    def trim(self):
        """
        Step trim(): Raw fastqs will be trimmed using Trimmomatic. Interleaved fastqs will be generated after trimming. 
        """

        jobs = []

        for sample in self.samples:   
           #log.info("Sample: " + sample.name)
           if not os.path.exists(os.path.join("genomes", sample.name, "qced_reads")):
               os.makedirs(os.path.join("genomes", sample.name, "qced_reads"))
           
           for readset in sample.readsets:
               #log.info("Readset: " + readset.sample.name) 
               #log.info("Readset: " + readset.fastq1) 
               #log.info("Readset: " + readset.fastq2) 
               #if not os.path.exists(os.path.join("genomes", sample.name, "qced_reads", readset.name)):
               #    os.makedirs(os.path.join("genomes", sample.name, "qced_reads", readset.name))
               
               trim_file_prefix = os.path.join("genomes", sample.name, "qced_reads", readset.name + "_trim_")
           
               job = genome_assembly.trimmomatic(
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
               job.name = "trimmomatic_" + readset.sample.name
               job.subname = "trim"
               jobs.append(job)

        return jobs

    def remove_contam(self):
        """
        Step remove_contam(): Trimmed fastqs will be filtered for contaminant sequences (e.g. Illumina adapters,
                              known primer sequences, etc). A second round of contaminant filtering will be done 
                              to filter out PhiX sequences which are usually spiked-in in Illumina sequencing runs.
        """
        
        jobs = []
        
        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "qced_reads")):
                os.makedirs(os.path.join("genomes", sample.name, "qced_reads"))
            
            for readset in sample.readsets:
                
                #trim_file_prefix = os.path.join("genomes", sample.name, "qced_reads", readset.name + "_trim_interleaved")
                
                log = os.path.join("genomes", sample.name, "qced_reads", readset.name + ".bbduk_contam_paired_log.txt")

                job = shotgun_metagenomics.bbduk_paired(
                    os.path.join("genomes", sample.name, "qced_reads", readset.name + "_trim_" + "pair1.fastq.gz"),
                    os.path.join("genomes", sample.name, "qced_reads", readset.name + "_trim_" + "pair2.fastq.gz"),
                    os.path.join("genomes", sample.name, "qced_reads", readset.name + "_" + "contam_R1.fastq"),
                    os.path.join("genomes", sample.name, "qced_reads", readset.name + "_" + "contam_R2.fastq"),
                    os.path.join("genomes", sample.name, "qced_reads", readset.name + "_ncontam_paired_R1.fastq.gz"),
                    os.path.join("genomes", sample.name, "qced_reads", readset.name + "_ncontam_paired_R2.fastq.gz"),
                    log,
                    config.param('DB', 'contaminants', 1, 'filepath')
                )
                job.name = "bbduk_paired_" + readset.sample.name
                job.subname = "duk"
                jobs.append(job)

                #job = genome_assembly.flash(
                #    trim_file_prefix + "_ncontam_nphix_paired.fastq.gz",
                #    os.path.join(sample.name, "qced_reads", readset.name),
                #    readset.name
                #)
                #job.name = "FLASH_" + sample.name + "_" + readset.name
                #job.subname = "flash"
                #jobs.append(job)

        return jobs

    def assembly(self):
        """
        Step assembly(): Perform assembly using Megahit and/or Spades
        """

        jobs = []
        
        for sample in self.samples:

            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            assemblers = assemblers_string
            reads_type = "pe"
           
            if not os.path.exists(os.path.join("genomes", sample.name, "assembly")):
                os.makedirs(os.path.join("genomes", sample.name, "assembly"))
            
            
            for assembler in assemblers:
                infiles = []
                infiles_merged = []
                infiles_R1 = []
                infiles_R2 = []

                if assembler == "megahit_merged":
                    for readset in sample.readsets:
                        trim_file_prefix = os.path.join("genomes", sample.name, "qced_reads", readset.name + "_trim")
                        infiles.append(trim_file_prefix + "_ncontam_nphix_paired.fastq.gz")
                        infiles_merged.append(os.path.join("genomes", sample.name, "qced_reads", readset.name + ".extendedFrags.fastq"))
                        infiles_R1.append(os.path.join("genomes", sample.name, "qced_reads", readset.name + ".notCombined_1.fastq"))
                    
                    job = genome_assembly.megahit(
                        infiles,
                        os.path.join("genomes", sample.name, "assembly", "megahit"),
                        reads_type
                    )
                    job.name = "megahit_" + sample.name
                    job.subname = "megahit"
                    jobs.append(job)
                    
                    job = genome_assembly.compile_assembly_stats(
                        os.path.join("genomes", sample.name, "assembly", "final.contigs.fa"),
                        os.path.join("genomes", sample.name, "assembly", "assembly_stats.txt")
                    )    
                    job.name = "compile_assembly_stats_megahit_" + sample.name
                    job.subname = "compile_assembly_stats"
                    jobs.append(job)

                elif assembler == "megahit":
                    for readset in sample.readsets:
                        infiles_R1.append(os.path.join("genomes", sample.name, "qced_reads", readset.name + "_ncontam_paired_R1.fastq.gz"))
                        infiles_R2.append(os.path.join("genomes", sample.name, "qced_reads", readset.name + "_ncontam_paired_R2.fastq.gz"))

                    job = shotgun_metagenomics.megahit_R1_R2(
                        infiles_R1,
                        infiles_R2,
                        os.path.join("genomes", sample.name, "assembly"),
                        reads_type
                    )
                    job.name = "megahit"
                    job.subname = "megahit"
                    jobs.append(job)

                    job = genome_assembly.compile_assembly_stats(
                        os.path.join("genomes", sample.name, "assembly", "final.contigs.fa"),
                        os.path.join("genomes", sample.name, "assembly", "assembly_stats.txt")
                    )    
                    job.name = "compile_assembly_stats_megahit_" + sample.name
                    job.subname = "compile_assembly_stats"
                    jobs.append(job)

                
                if assembler == "spades":
                    # Do spades --12
                    job = genome_assembly.spades(
                        infiles,
                        os.path.join("genomes", sample.name, "assembly"),
                        reads_type
                    )
                    job.name = "spades_" + sample.name
                    job.subname = "spades"
                    jobs.append(job)
                    
                    job = genome_assembly.compile_assembly_stats(
                        os.path.join("genomes", sample.name, "assembly", "contigs.fasta"),
                        os.path.join("genomes", sample.name, "assembly", "assembly_stats.txt")
                    )    
                    job.name = "compile_assembly_stats_spades_" + sample.name
                    job.subname = "compile_assembly_stats"
                    jobs.append(job)
                    
                if assembler == "spades_merged":
                    # Do spades --merged
                    job = genome_assembly.spades_merged(
                        infiles_merged,
                        infiles_R1,
                        infiles_R2,
                        os.path.join("genomes", sample.name, "assembly")
                    )
                    job.name = "spades_merged_" + sample.name
                    job.subname = "spades"
                    jobs.append(job)
                    
                    job = genome_assembly.compile_assembly_stats(
                        os.path.join("genomes", sample.name, "assembly", "contigs.fasta"),
                        os.path.join("genomes", sample.name, "assembly", "assembly_stats.txt")
                    )
                    job.name = "compile_assembly_stats_spades_merged" + sample.name
                    job.subname = "compile_assembly_stats"
                    jobs.append(job)


        return jobs


    def gene_prediction(self):
        """
        Step gene_prediction(): Use Prodigal to predict ORFs of bacterial genomes.
        """
        jobs = []
        
        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "gene_prediction")):
                os.makedirs(os.path.join("genomes", sample.name, "gene_prediction"))
            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            #assemblers = assemblers_string.split(":")
            #assemblers = ["megahit"]
            assemblers = assemblers_string
            
            for assembler in assemblers:
                if not os.path.exists(os.path.join("genomes", sample.name, "gene_prediction")):
                    os.makedirs(os.path.join("genomes", sample.name, "gene_prediction"))
            
            
                infile = os.path.join("genomes", sample.name, "assembly", "Contigs.fasta")
                outdir = os.path.join("genomes", sample.name, "gene_prediction")
                
                if not os.path.exists(outdir):
                    os.makedirs(os.path.join(outdir))
                
                job = shotgun_metagenomics.prodigal(
                    infile,
                    os.path.join(outdir, "Contigs.gff"),
                    os.path.join(outdir, "Contigs.fna"),
                    os.path.join(outdir, "Contigs.faa"),
                    os.path.join(outdir, "Contigs_renamed.gff"),
                    os.path.join(outdir, "Contigs_renamed.fna"),
                    os.path.join(outdir, "Contigs_renamed.faa")
                )
                job.name = "prodigal" + "_" + sample.name + "_" + assembler 
                job.subname = "prodigal"
                jobs.append(job)
 
        return jobs


    def coverage(self):
        """
        Step coverage(): Perform assembly using Megahit and/or Spades
        """

        jobs = []
        
        for sample in self.samples:
            #if not os.path.exists(os.path.join(sample.name, "qced_reads")):
            #    os.makedirs(os.path.join(sample.name, "qced_reads"))
            if not os.path.exists(os.path.join("genomes", sample.name, "contigs_abundance")):
                os.makedirs(os.path.join("genomes", sample.name, "contigs_abundance"))
            if not os.path.exists(os.path.join("genomes", sample.name, "gene_abundance")):
                os.makedirs(os.path.join("genomes", sample.name, "gene_abundance"))
            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            #assemblers = assemblers_string.split(":")
            #assemblers = ["megahit"]
            assemblers = assemblers_string

            for assembler in assemblers:
                if not os.path.exists(os.path.join("genomes", sample.name, "contigs_abundance")):
                    os.makedirs(os.path.join("genomes", sample.name, "contigs_abundance"))
                if not os.path.exists(os.path.join("genomes", sample.name, "gene_abundance")):
                    os.makedirs(os.path.join("genomes", sample.name, "gene_abundance"))

                # Will make index for bwa. and also bed file for computing reads spanning later.
                reference_contigs = os.path.join("genomes", sample.name, "assembly", "Contigs.fasta")
                bed_contigs = os.path.join("genomes", sample.name, "assembly", "Contigs.bed")
                bwt_contigs = os.path.join("genomes", sample.name, "assembly", "Contigs.fasta.bwt")
                bed_genes = os.path.join("genomes", sample.name, "gene_prediction", "Contigs_genes.bed")

                job = genome_assembly.make_index(
                    reference_contigs,
                    bwt_contigs
                )
                job.name = "make_index_contigs_" + assembler + "_" + sample.name
                job.subname = "make_index"
                jobs.append(job)
        
                job = genome_assembly.fasta_to_bed(
                    reference_contigs,
                    bed_contigs
                )
                job.name = "fasta_to_bed_" + assembler + "_" + sample.name
                job.subname = "fasta_to_bed"
                jobs.append(job)
        
                job = shotgun_metagenomics.gff_to_bed(
                    os.path.join("genomes", sample.name, "gene_prediction", "Contigs_renamed.gff"),
                    bed_genes
                )
                job.name = "gff_to_bed"
                job.subname = "gff_to_bed"
                jobs.append(job)

                flagstats_contigs_list = []
                trimmomatic_list = []
                cov_list_contigs = []
                cov_list_genes = []
                
                for readset in sample.readsets:
                    out1 = os.path.join("genomes", sample.name, "qced_reads", readset.name + "_ncontam_paired_R1.fastq.gz")
                    out2 = os.path.join("genomes", sample.name, "qced_reads", readset.name + "_ncontam_paired_R2.fastq.gz")
                    flag = "0x2"
                    bam_contigs = os.path.join("genomes", sample.name, "contigs_abundance", readset.name + ".bam")
                    flagstats_contigs = os.path.join("genomes", sample.name, "contigs_abundance", readset.name + ".flagstats")
                    flagstats_contigs_list.append(flagstats_contigs)
                    trimmomatic = os.path.join("genomes", sample.name, "qced_reads", readset.name + "_trim_stats.csv")
                    trimmomatic_list.append(trimmomatic)
                    cov_contigs = os.path.join("genomes", sample.name, "contigs_abundance", readset.name + ".cov")
                    cov_list_contigs.append(cov_contigs)
                    cov_genes = os.path.join("genomes", sample.name, "gene_abundance", readset.name + ".cov")
                    cov_list_genes.append(cov_genes)

                    ## map against contigs
                    job = genome_assembly.bwa_mem_samtools(
                        reference_contigs,
                        bwt_contigs,
                        out1,
                        out2,
                        bam_contigs
                    )
                    job.name = "bwa_mem-contigs" + readset.name + "_" + assembler
                    job.subname = "bwa"
                    jobs.append(job)

                    job = genome_assembly.flagstats(
                        bam_contigs,
                        flagstats_contigs
                    )
                    job.name = "flagstats" + readset.name + "_" + assembler
                    job.subname = "flagstats"
                    jobs.append(job)
            
                    job = genome_assembly.coverage_bed(
                        bam_contigs,
                        bed_contigs,
                        cov_contigs,
                        flag
                    )
                    job.name = "bedtoolsCov-contigs-" + readset.name + "_" + assembler
                    job.subname = "bedtools"
                    jobs.append(job)
 
                    job = shotgun_metagenomics.coverage_bed(
                        bam_contigs,
                        bed_genes,
                        cov_genes,
                        flag
                    )
                    job.name = "bedtoolsCov-genes-" + readset.name + "_" + assembler
                    job.subname = "bedtools"
                    jobs.append(job)
        
                job = shotgun_metagenomics.merge_counts(
                    cov_list_contigs,
                    os.path.join("genomes", sample.name, "contigs_abundance", "merged_contigs_abundance.tsv"),
                    os.path.join("genomes", sample.name, "contigs_abundance", "merged_contigs_abundance_cpm.tsv"),
                    "contigs"
                )
                job.name = "merge_contigs_abundance_contigs" + "_" + readset.name + "_" + assembler
                job.subname = "merge_contigs_abundance"
                jobs.append(job)
                
                job = shotgun_metagenomics.merge_counts(
                    cov_list_genes,
                    os.path.join("genomes", sample.name, "gene_abundance", "merged_gene_abundance.tsv"),
                    os.path.join("genomes", sample.name, "gene_abundance", "merged_gene_abundance_cpm.tsv"),
                    "genes"
                )
                job.name = "merge_gene_abundance_contigs" + "_" + readset.name + "_" + assembler
                job.subname = "merge_gene_abundance"
                jobs.append(job)
                
                job = genome_assembly.merge_flagstats(
                    flagstats_contigs_list,
                    trimmomatic_list,
                    os.path.join("genomes", sample.name, "contigs_abundance", "qc_mapping_stats.tsv")
                )
                job.name = "flagstats_merge_" + sample.name + "_" + assembler 
                job.subname = "flagstats"
                jobs.append(job)

                #Generate files for potential gam-ngs merge
                job = genome_assembly.get_insert_size(
                    os.path.join("genomes", sample.name, "contigs_abundance"),
                    os.path.join("genomes", sample.name, "contigs_abundance", "lib_stats.tsv"),
                    os.path.join("genomes", sample.name, "contigs_abundance", "qc_mapping_stats.tsv")
                )
                job.name = "get_insert_size_" + sample.name + "_" + assembler
                job.subname = "get_insert_size"
                jobs.append(job)

        return jobs

    def report_multiple(self):
        """
        Step report_muiltiple(): Generate report of all assemblies using Quast. This method is
                                 for the case when multiple assemblers are used and comparisons
                                 between assembly methods are needed.
        """
        
        jobs = []
        
        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "reports")):
                os.makedirs(os.path.join("genomes", sample.name, "reports"))
            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            assemblers = assemblers_string
          
            contigs_list = []   
            bam_list = []

            for assembler in assemblers:
                if not os.path.exists(os.path.join("genomes", sample.name, "contigs_abundance", assembler)):
                    os.makedirs(os.path.join("genomes", sample.name, "contigs_abundance", assembler))

                contigs_list.append(os.path.join("genomes", sample.name, "assembly", assembler, "Contigs.fasta"))

                if len(sample.readsets) == 1:
                    for readset in sample.readsets:
                        bam_list.append(os.path.join("genomes", sample.name, "contigs_abundance", assembler, readset.name + ".bam"))

                else:
                    bam_list_to_merge = []
                    
                    for readset in sample.readsets:
                        bam_contigs = os.path.join("genomes", sample.name, "contigs_abundance", assembler, readset.name + ".bam")
                        bam_list_to_merge.append(bam_contigs)
                    
                    ## map against contigs
                    job = genome_assembly.merge_bams(
                        bam_list_to_merge,
                        os.path.join("genomes", sample.name, "contigs_abundance", assembler, sample.name + ".bam")
                    )
                    job.name = "merge_bams_" + sample.name + "_" + assembler
                    job.subname = "bwa"
                    jobs.append(job)
                    
                    bam_list.append(os.path.join("genomes", sample.name, "contigs_abundance", assembler, sample.name + ".bam"))

        
            job = genome_assembly.quast(
                contigs_list,
                bam_list,
                os.path.join("genomes", sample.name, "reports")
            )
            job.name = "quast_" + sample.name
            job.subname = "flagstats"
            jobs.append(job)

        return jobs

    
    """
    Step exonerate(): Split fasta files in smaller chunks for annotation.
    """

    def exonerate(self):
        jobs = []
        
        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "gene_prediction")):
                os.makedirs(os.path.join("genomes", sample.name, "gene_prediction"))
            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            #assemblers = assemblers_string.split(":")
            #assemblers = ["megahit"]
            assemblers = assemblers_string
          
            contigs_list = []   
            bam_list = []

            for assembler in assemblers:
                if not os.path.exists(os.path.join("genomes", sample.name, "gene_prediction")):
                    os.makedirs(os.path.join("genomes", sample.name, "gene_prediction"))

                # Split for big assembly
                infile_contigs = os.path.join("genomes", sample.name, "assembly", "Contigs.fasta")
                chunks_dir_contigs = os.path.join("genomes", sample.name, "assembly", "fasta_chunks")
                number_of_chunks_file_contigs = os.path.join("genomes", sample.name, "assembly", "estimated_number_of_chunks_contigs.txt")
                infile_faa = os.path.join("genomes", sample.name, "gene_prediction", "Contigs_renamed.faa")
                chunks_dir = os.path.join("genomes", sample.name, "gene_prediction", "fasta_chunks")
                number_of_chunks_file_genes = os.path.join("genomes", sample.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
                
                #FNA contigs assembly (i.e. contigs)
                job = shotgun_metagenomics.estimate_number_of_chunks(
                    infile_contigs,
                    number_of_chunks_file_contigs,
                    str(config.param('exonerate', 'targeted_chunk_file_size_contigs', 1, 'int'))
                )
                job.name = "estimate_chunks_file_size" + "_" + sample.name + "_" + assembler
                jobs.append(job)

                job = shotgun_metagenomics.exonerate(
                    infile_contigs,
                    chunks_dir_contigs,
                    number_of_chunks_file_contigs,
                    "Contigs.fasta" 
                )
                job.name = "exonerate_contigs_fna" + "_" + sample.name + "_" + assembler
                job.subname = "exonerate"
                jobs.append(job)
                
                #FAA genes (i.e. predicted genes)
                job = shotgun_metagenomics.estimate_number_of_chunks(
                    infile_faa,
                    number_of_chunks_file_genes,
                    str(config.param('exonerate', 'targeted_chunk_file_size_genes', 1, 'int'))
                )
                job.name = "estimate_chunks_file_size" + "_" + sample.name + "_" + assembler
                jobs.append(job)
                
                job = shotgun_metagenomics.exonerate(
                    infile_faa,
                    chunks_dir,
                    number_of_chunks_file_genes,
                    "Contigs_renamed.faa" 
                )
                job.name = "exonerate_genes_fna" + "_" + sample.name + "_" + assembler
                job.subname = "exonerate"
                jobs.append(job)
        
        return jobs    
    

    def blastn_ncbi_genomes(self):
        """
        Step blastn_nt_contigs(): BLASTn contigs against NCBI nt
        """
        
        jobs = []

        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations"))
            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            #assemblers = assemblers_string.split(":")
            #assemblers = ["megahit"]
            assemblers = assemblers_string
          
            contigs_list = []   

            for assembler in assemblers:
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations"))
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "taxonomy_contigs")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations", "taxonomy_contigs"))

                # Do blastn on nt for big assembly  
                chunks_dir = os.path.join("genomes", sample.name, "assembly", "fasta_chunks")
                blast_dir = os.path.join("genomes", sample.name, "annotations", "blastn_nt_contigs")
                number_chunks_file = os.path.join("genomes", sample.name, "assembly", "estimated_number_of_chunks_contigs.txt")
                
                infiles = []
                outfiles = []
                dones = []
                
                num_chunks = 0
                if os.path.exists(number_chunks_file) and os.path.getsize(number_chunks_file) > 0:
                    with open(number_chunks_file, 'r') as file:
                        num_chunks = file.read().replace('\n', '')
                        num_chunks = int(num_chunks)
                else:
                    raise Exception(str(number_chunks_file) + " file does not exist\nPlease run exonerate step before running array jobs (blastn, diamond-blastp, hmmscan, rpsblast etc.\n")
                
                for i in range(num_chunks):
                    infiles.append(os.path.join(chunks_dir, "Contigs.fasta_chunk_{:07d}".format(i)))
                    outfiles.append(os.path.join(blast_dir, "blastn_chunk_{:07d}.tsv".format(i)))
                    dones.append(os.path.join(blast_dir, "blastn_chunk_{:07d}.done".format(i)))
             
                job = shotgun_metagenomics.blastn_array_job(
                    os.path.join(chunks_dir),
                    "Contigs.fasta_chunk_",
                    os.path.join(blast_dir),
                    "blastn_chunk_",
                    "blastn",
                    infiles,
                    outfiles,
                    dones,
                    self._curr_scheduler
                )
                job.name = "blastn_nt_contigs" + "_" + sample.name + "_" + assembler
                job.subname = "blastn"
                job.job_array_num_task = num_chunks
                jobs.append(job)
              
                # Merge output chunks
                job = shotgun_metagenomics.merge_chunks(
                    blast_dir,
                    os.path.join("genomes", sample.name, "annotations", "blastn_nt_contigs.tsv"),
                    num_chunks,
                    "blastn" 
                )
                job.name = "blastn_nt_big_contigs_merge" + "_" + sample.name + "_" + assembler
                job.subname = "blastn"
                jobs.append(job)
                
                job = shotgun_metagenomics.keep_blast_best_hit(
                    os.path.join("genomes", sample.name, "annotations", "blastn_nt_contigs.tsv"),
                    os.path.join("genomes", sample.name, "annotations", "blastn_nt_contigs_besthit.tsv")
                )
                job.name = "blastn_nt_big_contigs_best_hit" + "_" + sample.name + "_" + assembler
                job.subname = "keep_best_hit"
                jobs.append(job)
                    
                job = shotgun_metagenomics.extract_taxonomy(
                    os.path.join("genomes", sample.name, "annotations", "blastn_nt_contigs_besthit.tsv"),
                    os.path.join("genomes", sample.name, "annotations", "taxonomy_contigs", "taxonomy.tsv"),
                    config.param('ncbi_tax', 'accession_to_tax', required=True)
                )
                job.name = "extract_taxonomy_contigs" + "_" + sample.name + "_" + assembler
                job.subname = "ncbi_tax"
                jobs.append(job)
           
                #job = shotgun_metagenomics.generate_otu_table(
                #    os.path.join("genomes", sample.name, "annotations", assembler, "taxonomy_contigs", "taxonomy.tsv"),
                #    os.path.join("genomes", sample.name, "contigs_abundance", assembler, "merged_contigs_abundance.tsv"),
                #    os.path.join("genomes", sample.name, "annotations", assembler, "taxonomy_contigs", "otu_table.tsv")
                #)
                #job.name = "generate_otu_table_contigs" + "_" + sample.name + "_" + assembler
                #job.subname = "generate_otu_table"
                #jobs.append(job)
        
        return jobs
    
    def hmmsearch_kofam(self):
        jobs = []
        
        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "gene_prediction")):
                os.makedirs(os.path.join("genomes", sample.name, "gene_prediction"))
        
            chunks_dir = os.path.join("genomes", sample.name, "gene_prediction", "fasta_chunks")
            hmmsearch_outdir = os.path.join("genomes", sample.name, "annotations", "hmmsearch_kofam")
            number_chunks_file = os.path.join("genomes", sample.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "hmmsearch_kofam")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations", "hmmsearch_kofam"))
            
            infiles = []
            tblouts = []
            domtblouts = []
            pfamtblouts = []
            #outfiles = []
            dones = []
            
            num_chunks = 0
            if os.path.exists(number_chunks_file) and os.path.getsize(number_chunks_file) > 0:
                with open(number_chunks_file, 'r') as file:
                    num_chunks = file.read().replace('\n', '')
                    num_chunks = int(num_chunks)
            else:
                raise Exception(str(number_chunks_file) + " file does not exist\nPlease run exonerate step before running array jobs (blastn, diamond-blastp, hmmscan, rpsblast etc.\n")
            
            num_chunks = num_chunks 
            for i in range(0, num_chunks):
                infiles.append(os.path.join(chunks_dir, "Contigs_renamed.faa_chunk_{:07d}".format(i)))
                tblouts.append(os.path.join(hmmsearch_outdir, "hmmsearch_chunk_{:07d}.tblout".format(i)))
                domtblouts.append(os.path.join(hmmsearch_outdir, "hmmsearch_chunk_{:07d}.domtblout".format(i)))
                pfamtblouts.append(os.path.join(hmmsearch_outdir, "hmmsearch_chunk_{:07d}.pfamtblout".format(i)))
                dones.append(os.path.join(hmmsearch_outdir, "hmmsearch_chunk_{:07d}.done".format(i)))
            
            job = shotgun_metagenomics.hmmscan_array_job(
                os.path.join(chunks_dir),
                "Contigs_renamed.faa_chunk_",
                os.path.join(hmmsearch_outdir),
                "hmmsearch_chunk_",
                infiles,
                tblouts, domtblouts, pfamtblouts,
                dones,
                config.param('kegg', 'KO_profiles', 1, 'filepath'),
                self._curr_scheduler
            )
            job.name = "hmmsearch_kofam_array_job"
            job.subname = "hmmsearch"
            job.job_array_num_task = num_chunks
            jobs.append(job)

            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks_hmms(
                hmmsearch_outdir,
                os.path.join("genomes", sample.name, "annotations"),
                num_chunks,
                "hmmsearch",
                "hmmsearch_kofam"
            )
            job.name = "hmmsearch_kofam_merge"
            job.subname = "merge"      
            jobs.append(job)
            
            job = shotgun_metagenomics.parse_kofam(
                os.path.join("genomes", sample.name, "annotations", "hmmsearch_kofam_tblout.tsv"),
                os.path.join("genomes", sample.name, "annotations", "KOs_parsed.tsv"),
                hmmsearch=True
            )
            job.name = "parse_kegg"
            job.subname = "merge"
            jobs.append(job)
        
        return jobs 

    
    def diamond_blastp_kegg(self):
        jobs = []
        
        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations"))
            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            #assemblers = assemblers_string.split(":")
            #assemblers = ["megahit"]
            assemblers = assemblers_string
          
            contigs_list = []   

            for assembler in assemblers:
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations"))

                # TO FIX eventually
                fname = os.path.join("genomes", sample.name, "annotations", "blastp_nr_annotated.tsv")
                open(fname, 'a').close()
                os.utime(fname, None)
                
                # Do blastp on KEGG for co-assembly 
                chunks_dir = os.path.join("genomes", sample.name, "gene_prediction", "fasta_chunks")
                blast_dir = os.path.join("genomes", sample.name, "annotations", "blastp_kegg")
                number_chunks_file = os.path.join("genomes", sample.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
                
                infiles = []
                outfiles = []
                dones = []
                
                num_chunks = 0
                if os.path.exists(number_chunks_file) and os.path.getsize(number_chunks_file) > 0:
                    with open(number_chunks_file, 'r') as file:
                        num_chunks = file.read().replace('\n', '')
                        num_chunks = int(num_chunks)
                else:
                    raise Exception(str(number_chunks_file) + " file does not exist\nPlease run exonerate step before running array jobs (blastn, diamond-blastp, hmmscan, rpsblast etc.\n")

                for i in range(num_chunks):
                    infiles.append(os.path.join(chunks_dir, "Contigs_renamed.faa_chunk_{:07d}".format(i)))
                    outfiles.append(os.path.join(blast_dir, "blastp_chunk_{:07d}.tsv".format(i)))
                    dones.append(os.path.join(blast_dir, "blastp_chunk_{:07d}.done".format(i)))
              
                job = shotgun_metagenomics.diamond_blastp_kegg_array_job(
                    os.path.join(chunks_dir),
                    "Contigs_renamed.faa_chunk_",
                    os.path.join(blast_dir),
                    "blastp_chunk_",
                    "blastp",
                    infiles,
                    outfiles,
                    config.param('blastp_kegg', 'db_kegg', 1, 'string'),
                    dones,
                    self._curr_scheduler
                )
                job.name = "blastp_array_kegg" + "_" + sample.name + "_" + assembler
                job.subname = "blastp_kegg"
                job.job_array_num_task = num_chunks
                jobs.append(job)
                
                # Merge output chunks
                job = shotgun_metagenomics.merge_chunks(
                    blast_dir,
                    os.path.join("genomes", sample.name, "annotations", "blastp_kegg.tsv"),
                    num_chunks,
                    "blastp"
                )
                job.name = "blastp_kegg" + "_" + sample.name + "_" + assembler
                job.subname = "merge"
                jobs.append(job)
                
                # Generate a clean table of module/pathways.
                job = shotgun_metagenomics.parse_kegg(
                    os.path.join("genomes", sample.name, "annotations", "blastp_kegg.tsv"),
                    os.path.join("genomes", sample.name, "annotations", "blastp_kegg_parsed.tsv")
                )
                job.name = "parse_kegg" + "_" + sample.name + "_" + assembler
                job.subname = "merge"
                jobs.append(job)

        return jobs 

    """
    Step rpsblast_cog(): Perform RPSBLAST of genes (faa) against NCBI COG database.
    """

    def rpsblast_cog(self):
        jobs = []
        
        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations"))
            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            #assemblers = assemblers_string.split(":")
            #assemblers = ["megahit"]
            assemblers = assemblers_string
          
            contigs_list = []   

            for assembler in assemblers:
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations"))
        
                number_chunks_file = os.path.join("genomes", sample.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
                num_chunks = 0
                if os.path.exists(number_chunks_file) and os.path.getsize(number_chunks_file) > 0:
                    with open(number_chunks_file, 'r') as file:
                        num_chunks = file.read().replace('\n', '')
                        num_chunks = int(num_chunks)
                else:
                    raise Exception(str(number_chunks_file) + " file does not exist\nPlease run exonerate step before running array jobs (blastn, diamond-blastp, hmmscan, rpsblast etc.\n")
                
                # Do rpsblast on COG for big assembly  
                if(config.param('blastp', 'skip_cog', 1, 'string') == 'yes'):
                    fname = os.path.join("genomes", sample.name, "annotations", "rpsblast_cog.tsv")
                    open(fname, 'a').close()
                    os.utime(fname, None)
                else:
                    chunks_dir = os.path.join("genomes", sample.name, "gene_prediction", "fasta_chunks")
                    blast_dir = os.path.join("genomes", sample.name, "annotations", "rpsblast_cog")
                    #num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
                    infiles = []
                    outfiles = []
                    dones = []
                    
                    for i in range(num_chunks):
                        infiles.append(os.path.join(chunks_dir, "Contigs_renamed.faa_chunk_{:07d}".format(i)))
                        outfiles.append(os.path.join(blast_dir, "rpsblast_chunk_{:07d}.tsv".format(i)))
                        dones.append(os.path.join(blast_dir, "rpsblast_chunk_{:07d}.done".format(i)))
            
                    job = shotgun_metagenomics.rpsblast_array_job(
                        os.path.join(chunks_dir),
                        "Contigs_renamed.faa_chunk_",
                        os.path.join(blast_dir),
                        "rpsblast_chunk_",
                        infiles,
                        outfiles,
                        dones,
                        config.param('cog', 'db'),
                        self._curr_scheduler
                    )
                    job.name = "rpsblast_cog" + "_" + sample.name + "_" + assembler
                    job.subname = "rpsblast"
                    job.job_array_num_task = num_chunks
                    jobs.append(job)
            
                    # Merge output chunks
                    job = shotgun_metagenomics.merge_chunks(
                        blast_dir,
                        os.path.join("genomes", sample.name, "annotations", "rpsblast_cog.tsv"),
                        num_chunks,
                        "rpsblast"
                    )
                    job.name = "rpsblast_cog" + "_" + sample.name + "_" + assembler
                    job.subname = "merge"
                    jobs.append(job)
        
        return jobs 
   
    """
    Step hmmacan_pfam(): HMMSCAN of all predicted genes against PFAM-A database.
    """

    def hmmscan_pfam(self):
        
        jobs = []

        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations"))
            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            #assemblers = assemblers_string.split(":")
            #assemblers = ["megahit"]
            assemblers = assemblers_string
          
            for assembler in assemblers:
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations"))
        
                chunks_dir = os.path.join("genomes", sample.name, "gene_prediction", "fasta_chunks")
                blast_dir = os.path.join("genomes", sample.name, "annotations", "hmmscan_pfam")
                number_chunks_file = os.path.join("genomes", sample.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
                
                infiles = []
                tblouts = []
                domtblouts = []
                pfamtblouts = []
                dones = []
                
                num_chunks = 0
                if os.path.exists(number_chunks_file) and os.path.getsize(number_chunks_file) > 0:
                    with open(number_chunks_file, 'r') as file:
                        num_chunks = file.read().replace('\n', '')
                        num_chunks = int(num_chunks)
                else:
                    raise Exception(str(number_chunks_file) + " file does not exist\nPlease run exonerate step before running array jobs (blastn, diamond-blastp, hmmscan, rpsblast etc.\n")
                
                for i in range(num_chunks):
                    infiles.append(os.path.join(chunks_dir, "Contigs_renamed.faa_chunk_{:07d}".format(i)))
                    tblouts.append(os.path.join(blast_dir, "hmmscan_chunk_{:07d}.tblout".format(i)))
                    domtblouts.append(os.path.join(blast_dir, "hmmscan_chunk_{:07d}.domtblout".format(i)))
                    pfamtblouts.append(os.path.join(blast_dir, "hmmscan_chunk_{:07d}.pfamtblout".format(i)))
                    dones.append(os.path.join(blast_dir, "hmmscan_chunk_{:07d}.done".format(i)))

                job = shotgun_metagenomics.hmmscan_array_job(
                    os.path.join(chunks_dir),
                    "Contigs_renamed.faa_chunk_",
                    os.path.join(blast_dir),
                    "hmmscan_chunk_",
                    infiles,
                    tblouts, domtblouts, pfamtblouts,
                    dones,
                    config.param('pfam', 'db', required=True),
                    self._curr_scheduler
                )
                job.name = "hmmscan_pfam" + "_" + sample.name + "_" + assembler
                job.subname = "hmmscan"
                job.job_array_num_task = num_chunks
                jobs.append(job)

                # Merge output chunks
                job = shotgun_metagenomics.merge_chunks_hmms(
                    blast_dir,
                    os.path.join("genomes", sample.name, "annotations"),
                    num_chunks,
                    "hmmscan",
                    "hmmscan_pfam"
                )
                job.name = "hmmscan_pfam_merge" + "_" + sample.name + "_" + assembler
                job.subname = "merge"      
                jobs.append(job)
            
        return jobs

    def diamond_blastp_nr(self):
        jobs = []
       
        # For now, create/touch dummy file. To improve/fix eventually.
        #if(config.param('blastp', 'skip_blastp_nr', 1, 'string') == 'yes'):
        #sys.stderr.write("Skipping diamond_blastp_nr step...\n")
        #fname = os.path.join(self._root_dir, "annotations", "blastp_nr_annotated.tsv")
        #open(fname, 'a').close()
        #os.utime(fname, None)
        
        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations"))
            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            #assemblers = assemblers_string.split(":")
            #assemblers = ["megahit"]
            assemblers = assemblers_string
          
            for assembler in assemblers:
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations"))
                #if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "taxonomy_genes")):
                #    os.makedirs(os.path.join("genomes", sample.name, "annotations", "taxonomy_genes"))
        
                # Do blastp on KEGG for big assembly  
                chunks_dir = os.path.join("genomes", sample.name, "gene_prediction", "fasta_chunks")
                blast_dir = os.path.join("genomes", sample.name, "annotations", "blastp_nr")
                number_chunks_file = os.path.join("genomes", sample.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
                
                infiles = []
                outfiles = []
                dones = []

                num_chunks = 0
                if os.path.exists(number_chunks_file) and os.path.getsize(number_chunks_file) > 0:
                    with open(number_chunks_file, 'r') as file:
                        num_chunks = file.read().replace('\n', '')
                        num_chunks = int(num_chunks)
                else:
                    raise Exception(str(number_chunks_file) + " file does not exist\nPlease run exonerate step before running array jobs (blastn, diamond-blastp, hmmscan, rpsblast etc.\n")

                for i in range(num_chunks):
                    infiles.append(os.path.join(chunks_dir, "Contigs_renamed.faa_chunk_{:07d}".format(i)))
                    outfiles.append(os.path.join(blast_dir, "blastp_chunk_{:07d}.tsv".format(i)))
                    dones.append(os.path.join(blast_dir, "blastp_chunk_{:07d}.done".format(i)))
                
                job = shotgun_metagenomics.diamond_blastp_nr_array_job(
                    os.path.join(chunks_dir),
                    "Contigs_renamed.faa_chunk_",
                    os.path.join(blast_dir),
                    "blastp_chunk_",
                    "blastp",
                    infiles,
                    outfiles,
                    config.param('blastp_nr', 'db_nr', 1, 'string'),
                    dones,
                    self._curr_scheduler
                )
                job.name = "blastp_array_nr" + "_" + sample.name + "_" + assembler
                job.subname = "blastp_nr"
                job.job_array_num_task = num_chunks
                jobs.append(job)
                
                # Merge output chunks
                job = shotgun_metagenomics.merge_chunks(
                    blast_dir,
                    os.path.join("genomes", sample.name, "annotations", "blastp_nr.tsv"),
                    num_chunks,
                    "blastp"
                )
                job.name = "merge_diamond_blastp_nr" + "_" + sample.name + "_" + assembler
                job.subname = "merge"
                jobs.append(job)
                
        
    
        return jobs 
   
    def taxonomic_annotation(self):
        """
        Step taxonomic_annotation(): Implements the CAT package to retrieve consensus taxonomy of each contig.
        """

        jobs = []
        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations"))
            
            assemblers_string = config.param('DEFAULT', 'assembler', 1, 'string').split(":")
            #assemblers = assemblers_string.split(":")
            #assemblers = ["megahit"]
            assemblers = assemblers_string
          
            for assembler in assemblers:
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations"))
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus"))
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "all")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "all"))
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "all", "absolute")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "all", "absolute"))
                if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "all", "relative")):
                    os.makedirs(os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "all", "relative"))
        
                job = genome_assembly.convert_orf_ids_for_cat(
                   os.path.join("genomes", sample.name, "gene_prediction", "Contigs_renamed.gff"),
                   os.path.join("genomes", sample.name, "annotations",  "blastp_nr.tsv"),
                   os.path.join("genomes", sample.name, "annotations",  "blastp_nr_alt_orf_ids.tsv")
                )
                job.name = "convert_ord_ids_for_CAT"
                job.subname = "convert_ord_ids_for_CAT"
                jobs.append(job)
                
                job = genome_assembly.CAT(
                   os.path.join("genomes", sample.name, "assembly", "Contigs.fasta"),
                   os.path.join("genomes", sample.name, "gene_prediction","Contigs.faa"),
                   os.path.join("genomes", sample.name, "annotations", "blastp_nr_alt_orf_ids.tsv"),
                   os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "out")
                )
                job.name = "CAT"
                job.subname = "CAT"
                jobs.append(job)
                
                #job = genome_assembly.generate_otu_table_consensus(
                job = shotgun_metagenomics.generate_feature_table_consensus(
                   os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "out.contig2classification_with_names.tsv"),
                   os.path.join("genomes", sample.name, "contigs_abundance", "merged_contigs_abundance.tsv"),
                   os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "taxonomy.tsv"),
                   os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "feature_table.tsv")

                )
                job.name = "generate_feature_table_consensus"
                job.subname = "generate_feature_table"
                jobs.append(job)
           
                prefixes = [
                    "consensus"
                ]
                infiles = [
                    os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "feature_table.tsv")
                ]
                directories = [
                    os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus")
                ]
                abundances = [
                    os.path.join("genomes", sample.name, "contigs_abundance", "merged_contigs_abundance.tsv") # this file is only used to generate the normalized otu table downstream 
                ] 
          
                # Then proces all otu tables including the extended ones if -e option
                for i in range(0, len(infiles)):
                    # This job doesnt really do any normalization at all. It takes cpm values and creates an otu table.
                    normalizations = [""]
                    types = ["absolute", "relative"]
                    #organisms = ["", "_others", "_bacteriaArchaea"]
                    #organisms2 = ["all", "others", "bacteriaArchaea"]
                    organisms = [""]
                    organisms2 = ["all"]

                    for n in range(0, len(normalizations)):
                        for j in range(0, len(types)):
                            for k in range(0, len(organisms)):
                                for m in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
                                    if types[j] == "absolute":
                                        job = microbial_ecology.summarize_taxonomy_absolute(
                                            os.path.join(directories[i], "feature_table" + normalizations[n] + organisms[k] + ".tsv"),
                                            m,
                                            os.path.join(directories[i], organisms2[k], types[j]) 
                                        )
                                        job.name = "tax_summary" + normalizations[n] + "_" + types[j] + "_" + organisms2[k]  + "_L" + str(m) + "_" + prefixes[i]
                                        job.subname = "summarize_taxonomy"
                                        jobs.append(job)
                                    elif types[j] == "relative":
                                        job = microbial_ecology.summarize_taxonomy_relative(
                                            os.path.join(directories[i], "feature_table" + normalizations[n] + organisms[k] + ".tsv"),
                                            m,
                                            os.path.join(directories[i], organisms2[k], types[j]) 
                                        )
                                        job.name = "tax_summary" + normalizations[n] + "_" + types[j] + "_" + organisms2[k]  + "_L" + str(m) + "_" + prefixes[i]
                                        job.subname = "summarize_taxonomy"
                                        jobs.append(job)
                        
                    
        return jobs
   

    def ncrna(self):
        """
        Step ncrna(): scan for non-coding RNA in contigs: rRNAs only
        """

        jobs = []
        
        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "rrna")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations", "rrna"))
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap"))
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "rrna", "abundance")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations", "rrna", "abundance"))
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "rrna", "rdp")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations", "rrna", "rdp"))
            # Do rRNA prediction from co-assembly with barrnap
            prefixes = ["bac", "arc", "euk", "mito"]
            for prefix in prefixes:
                job = shotgun_metagenomics.barrnap(
                    os.path.join("genomes", sample.name, "assembly", "Contigs.fasta"),
                    prefix,
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + ".fna"),
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + ".tsv")
                )
                job.name = "barrnap_" + prefix
                job.subname = "barrnap"
                jobs.append(job)
           
                job = shotgun_metagenomics.split_barrnap_fna(
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + ".fna"),
                    prefix,
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + "_5S.fna"), 
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + "_SSU.fna"), 
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + "_LSU.fna"), 
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + "_12S.fna"), 
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + "_5-8S.fna") 
                )
                job.name = "split_barrnap_fna_" + prefix + "_" + sample.name
                jobs.append(job)
        
            prefixes = ["arc_SSU", "arc_LSU", "bac_SSU", "bac_LSU", "euk_SSU", "euk_LSU"]
            for prefix in prefixes:
                # create bed files for abundance - we have to only get the rea
                # Just take taxonomy file which contains coord for each rRNA genes, then generate bed file.
                job = shotgun_metagenomics.barrnap_fna_to_bed(
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + ".fna"),
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + ".bed")
                )
                job.name = "barrnap_to_bed_" + prefix + "_" + sample.name
                job.subname = "barrnap_to_bed"
                jobs.append(job)

                # Then from this new bed file, process all bams to just get reads that falls into coords of new bed file
                # which corresponds to rRNA genes. Abundance will be used for both blastn and rdp taxonomy.
                cov_list = []
                for readset in sample.readsets:
                
                    if(readset.run_type == "SINGLE_END"):
                        flag = "0x0"
                    elif(readset.run_type == "PAIRED_END"):
                        flag = "0x2"

                    bam_contigs = os.path.join("genomes", sample.name, "contigs_abundance", readset.name + ".bam")
                
                    job = shotgun_metagenomics.coverage_bed(
                        bam_contigs,
                        os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + ".bed"),
                        os.path.join("genomes", sample.name, "annotations", "rrna", "abundance", prefix + "_" + readset.name + ".cov"),
                        flag,
                        prefix
                    )
                    job.name = "bedtoolsCov-contigs-rrna_" + readset.sample.name + "_" + prefix
                    job.subname = "bedtools_rrna"
                    jobs.append(job)
                    
                    job = shotgun_metagenomics.rename_rrna_bed(
                        os.path.join("genomes", sample.name, "annotations", "rrna", "abundance", prefix + "_" + readset.name + ".cov"),
                        os.path.join("genomes", sample.name, "annotations", "rrna", "abundance", prefix + "_" + readset.name + "_renamed.cov"),
                    )
                    job.name = "rename_bed_" + readset.sample.name + "_" + prefix
                    job.subname = "rename_bed"
                    jobs.append(job)

                    cov_list.append(os.path.join("genomes", sample.name, "annotations", "rrna", "abundance", prefix + "_" + readset.name + "_renamed.cov"))
            
                job = shotgun_metagenomics.merge_counts(
                    cov_list,
                    os.path.join("genomes", sample.name, "annotations", "rrna", "abundance", prefix + "_merged_abundance.tsv"),
                    os.path.join("genomes", sample.name, "annotations", "rrna", "abundance", prefix + "_merged_abundance_cpm.tsv"),
                    "contigs",
                    True
                )
                job.name = "merge_contig_abundance_" + prefix + "_" + sample.name
                job.subname = "merge_gene_abundance"
                jobs.append(job)
                
                # RDP classifier with rrna results.
                job = microbial_ecology.rdp_wrapper(
                    os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + ".fna"),
                    os.path.join("genomes", sample.name, "annotations", "rrna", "rdp", prefix + "_rdp.tsv"),
                    0
                )
                job.name = "classify_" + prefix + "_" + sample.name
                job.subname = "RDP"
                jobs.append(job)
               
                # Convert RDP table to in-house taxonomy format.
                job = shotgun_metagenomics.rdp_to_taxonomy(
                    os.path.join("genomes", sample.name, "annotations", "rrna", "rdp", prefix + "_rdp.tsv"),
                    os.path.join("genomes", sample.name, "annotations", "rrna", "rdp", prefix + "_rdp_taxonomy.tsv")
                )
                job.name = "rdp_to_taxonomy_rrna_" + prefix + "_" + sample.name
                job.subname = "rdp_to_taxonomy"
                jobs.append(job)

                # Then with rdp output, generate otu table
                job = shotgun_metagenomics.generate_feature_table_rrna(
                    os.path.join("genomes", sample.name, "annotations", "rrna", "rdp", prefix + "_rdp_taxonomy.tsv"),
                    os.path.join("genomes", sample.name, "annotations", "rrna", "abundance", prefix + "_merged_abundance.tsv"),
                    os.path.join("genomes", sample.name, "annotations", "rrna", "rdp", prefix +  "_feature_table.tsv")
                )
                job.name = "generate_feature_table_" + prefix + "_" + sample.name
                job.subname = "generate_feature_table"
                jobs.append(job)

        return jobs


    def finalize(self):
        """
        Step finalize(): Generate summary tsv and gff files compiling all annotation data.
                         Also generate global feature table and file containing all 16S sequences.
        """
        
        jobs = []
           
        feature_tables = []
        fasta_16S = []
        prefixes = []
        
        if not os.path.exists(os.path.join("annotations")):
            os.makedirs(os.path.join("annotations"))
        if not os.path.exists(os.path.join("annotations", "absolute")):
            os.makedirs(os.path.join("annotations", "absolute"))
        if not os.path.exists(os.path.join("annotations", "relative")):
            os.makedirs(os.path.join("annotations", "relative"))

        for sample in self.samples:
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations")):
                os.makedirs(os.path.join("genomes", sample.name, "annotations"))
                
            if not (os.path.join("genomes", sample.name, "annotations", "rpsblast_kog.tsv")):
                fname = os.path.join("genomes", sample.name, "annotations", "rpsblast_kog.tsv")
                open(fname, 'a').close()
                os.utime(fname, None)
            
            #if not (os.path.join("genomes", sample.name, "annotations", "blastp_nr_annotated.tsv")):
            #    fname = os.path.join("genomes", sample.name, "annotations", "blastp_nr_annotated.tsv")
            #    sys.stderr.write(fname + "\n")
            #    open(fname, 'a').close()
            #    os.utime(fname, None)
            if not os.path.exists(os.path.join("genomes", sample.name, "annotations", "blastp_nr_annotated.tsv")):
                os.mknod(os.path.join("genomes", sample.name, "annotations", "blastp_nr_annotated.tsv"))


            job = genome_assembly.generate_gff(
                # infiles
                os.path.join("genomes", sample.name, "gene_prediction", "Contigs_renamed.gff"),
                os.path.join("genomes", sample.name, "assembly", "Contigs.fasta"),
                os.path.join("genomes", sample.name, "annotations", "KOs_parsed.tsv"),
                os.path.join("genomes", sample.name, "annotations", "hmmscan_pfam_tblout.tsv"),
                os.path.join("genomes", sample.name, "annotations", "rpsblast_cog.tsv"),
                os.path.join("genomes", sample.name, "annotations", "rpsblast_kog.tsv"),
                os.path.join("genomes", sample.name, "annotations", "taxonomy_consensus", "taxonomy.tsv"),
                os.path.join("genomes", sample.name, "annotations", "blastp_nr_annotated.tsv"),
                # outfiles
                os.path.join("genomes", sample.name, "annotations", "Contigs.gff"),
                os.path.join("genomes", sample.name, "annotations", "Contigs.fasta"),
                os.path.join("genomes", sample.name, "annotations", "annotations.tsv")
            )
            job.name = "generate_gff" + "_" + sample.name
            job.subname = "generate_gff"
            jobs.append(job)
        
            prefix = "bac_SSU"
            feature_tables.append(os.path.join("genomes", sample.name, "annotations", "rrna", "rdp", prefix +  "_feature_table.tsv"))
            fasta_16S.append(os.path.join("genomes", sample.name, "annotations", "rrna", "barrnap", prefix + ".fna"))
            prefixes.append(sample.name)

        # generate 16S centric results.
        job = genome_assembly.merge_feature_tables(
            prefixes,
            feature_tables,
            os.path.join("annotations", "feature_table_16S.tsv")
        )
        job.name = "merge_feature_table_16S"
        job.subname = "merge"
        jobs.append(job)
        
        job = genome_assembly.merge_fasta(
            prefixes,
            fasta_16S,
            os.path.join("annotations", "barrnap_16S.fna")
        )
        job.name = "merge_fasta_16S"
        job.subname = "merge"
        jobs.append(job)
                                
        types = ["absolute", "relative"]
        for j in range(0, len(types)):
            for m in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
                if types[j] == "absolute":
                    job = microbial_ecology.summarize_taxonomy_absolute(
                        os.path.join("annotations", "feature_table_16S.tsv"),
                        m,
                        os.path.join("annotations", types[j])#, "feature_table_16S_" +  types[j] + "L" + str(m) + ".tsv") 
                    )
                    job.name = "tax_summary" + "_" + types[j] + "_L" + str(m)
                    job.subname = "summarize_taxonomy"
                    jobs.append(job)
                elif types[j] == "relative":
                    job = microbial_ecology.summarize_taxonomy_relative(
                        os.path.join("annotations", "feature_table_16S.tsv"),
                        m,
                        os.path.join("annotations", types[j])#, "feature_table_16S_" +  types[j] + "L" + str(m) + ".tsv") 
                    )
                    job.name = "tax_summary" + "_" + types[j] + "_L" + str(m)
                    job.subname = "summarize_taxonomy"
                    jobs.append(job)

        return jobs
    
    #def set_local_variables(self):
        #self._parser_local = argparse.ArgumentParser(description='Process options.')
        #self._args_local = self._parser_local.parse_args()
        #self._root_dir = self._args_local.output_dir

    def __init__(self):
        version = open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "..", "VERSION"), 'r').read().split('\n')[0]
        ga_string = """
######################################################################################
   ____                                    _                           _     _       
  / ___| ___ _ __   ___  _ __ ___   ___   / \   ___ ___  ___ _ __ ___ | |__ | |_   _ 
 | |  _ / _ \ '_ \ / _ \| '_ ` _ \ / _ \ / _ \ / __/ __|/ _ \ '_ ` _ \| '_ \| | | | |
 | |_| |  __/ | | | (_) | | | | | |  __// ___ \\\\__ \__ \  __/ | | | | | |_) | | |_| |
  \____|\___|_| |_|\___/|_| |_| |_|\___/_/   \_\___/___/\___|_| |_| |_|_.__/|_|\__, |
                                                                               |___/ 
                 Support: jtremblay514@gmail.com
               Home page: jtremblay.github.io/pipelines.html
                 Version: """ + version + """

######################################################################################"""
        sys.stderr.write(ga_string + '\n')
        time.sleep(1)
        #self.set_local_variables()
        self.argparser.add_argument("-r", "--readsets", help="readset file",  type=argparse.FileType('r'))
        self._parser_local = self.argparser
        self._args_local = self._parser_local.parse_args() 
        self._curr_scheduler = self._args_local.job_scheduler
        super(GenomeAssembly, self).__init__()

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            if self.args.readsets:
                self._readsets = parse_readset_file(self.args.readsets.name)
            else:
                self.argparser.error("argument -r/--readsets is required!")

        return self._readsets

    @property
    def steps(self):
        return [
            self.trim,
            self.remove_contam,
            self.assembly,
            self.gene_prediction,
            self.coverage,
            #elf.report_multiple
            self.exonerate,
            self.blastn_ncbi_genomes,
            self.hmmsearch_kofam,#self.diamond_blastp_kegg,
            self.rpsblast_cog,
            self.hmmscan_pfam,
            self.diamond_blastp_nr,
            self.taxonomic_annotation,
            self.ncrna,
            self.finalize
        ]

#if __name__ == '__main__':
#    GenomeAssembly()

GenomeAssembly().submit_jobs()
