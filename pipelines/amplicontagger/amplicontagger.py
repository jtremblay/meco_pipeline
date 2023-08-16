#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys
import time
from pipelines import common

# Append meco_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# MECO Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.readset import *

from bio import rrna_amplicons
from bio import microbial_ecology
from pipelines import common

# Global scope variables
log = logging.getLogger(__name__)

class AmpliconTagger(common.MECOPipeline):

    def validate(self):
        jobs = []

        job = rrna_amplicons.validate_map_and_barcodes(
            os.path.relpath(self._barcodes.name),
            os.path.join(self._root_dir, "validateMapAndBarcodes.done")
        )
        job.name = "validate_map_and_barcodes"
        jobs.append(job)

        return jobs

    def remove_contam(self):
        jobs = []
        
        validated_barcodes = os.path.join(self._root_dir, "validateMapAndBarcodes.done")
        
        if 'self._external_infile' in globals():
            infile = self._external_infile
        else:
            infile = self._infile

        sys.stderr.write('external infile: ' + infile + '\n')

        job = rrna_amplicons.bbduk(
            infile,
            os.path.join(self._duk_dir, "contam.fastq"),
            os.path.join(self._duk_dir, "ncontam.fastq"),
            os.path.join(self._duk_dir, "ncontam_singletons.fastq"),
            os.path.join(self._duk_dir, "duk_contam_log.txt"),
            config.param('DB', 'contaminants', 1, 'filepath'),
            validated_barcodes
        );
        job.name = "bbduk_contam"
        job.subname = "bbduk"
        jobs.append(job)
        
        job = rrna_amplicons.bbduk(
            os.path.join(self._duk_dir, "ncontam.fastq"),
            os.path.join(self._duk_dir, "ncontam_phix.fastq"),
            os.path.join(self._duk_dir, "ncontam_nphix.fastq"),
            os.path.join(self._duk_dir, "ncontam_nphix_singletons.fastq"),
            os.path.join(self._duk_dir, "duk_phix_log.txt"),
            config.param('DB', 'phix', 'filepath'),
            False
        );
        job.name = "bbduk_phix"
        job.subname = "bbduk"
        jobs.append(job)

        return jobs

    def split_barcodes(self):
        jobs = []
        
        if(self._lib_type == "nc1"):
            barcodes_output = self._root_dir + "/reads_1/fastqs/ncontam_nphix_1.fastq"
        elif(self._lib_type == "nc2"):
            barcodes_output = self._root_dir + "/reads_2/fastqs/ncontam_nphix_2.fastq"
        elif(self._lib_type == "ex" or self._lib_type == "nc1nc2"):
            barcodes_output = self._root_dir + "/reads_12/fastqs/ncontam_nphix.fastq"

        # split barcodes
        job = rrna_amplicons.split_barcodes(
          self._duk_dir + "/ncontam_nphix.fastq",
          os.path.relpath(self._barcodes.name),
          barcodes_output,
          self._lib_type_dir + "/fastqs/ncontam_nphix_barcodes.txt"
        )
        job.name = "split_barcodes"
        job.subname = "barcodes"
        jobs.append(job)
        
        # Only performs remove_unpaired and split pairs if library is ex (i.e. overlapping paired end reads).
        if(self._lib_type == "ex" or self._lib_type == "nc1nc2"):
            # remove unpaired reads.
            job = rrna_amplicons.remove_unpaired_reads(
                self._lib_type_dir + "/fastqs/ncontam_nphix.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_paired.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_unpaired_1.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_unpaired_2.fastq"
            )
            job.name = "remove_unpaired"
            job.subname = "remove_unpaired"
            jobs.append(job)
        
            # split by pairs.
            job = rrna_amplicons.split_pairs(
                self._lib_type_dir + "/fastqs/ncontam_nphix_paired.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_1.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_2.fastq"
            )
            job.name = "split_pairs"
            job.subname = "split_pairs"
            jobs.append(job)

        return jobs    

    def qscores(self):
        jobs = []
        
        # Do reads 1
        if(self._lib_type == "ex" or self._lib_type == "nc1"):
            job = rrna_amplicons.generate_qscore_sheet(
                self._lib_type_dir + "/fastqs/ncontam_nphix_1.fastq",
                "reads_1",
                self._lib_type_dir + "/qscores/qual_stats_1_log.txt",
                self._lib_type_dir + "/qscores/qscores_1.tsv",
                os.path.relpath(self._barcodes.name)
            )
            job.name = "generate_qscore_reads1"
            job.subname = "qscore_sheet"
            jobs.append(job)
        
        # Do reads 2
        if(self._lib_type == "ex" or self._lib_type == "nc2"):
            job = rrna_amplicons.generate_qscore_sheet(
                self._lib_type_dir + "/fastqs/ncontam_nphix_2.fastq",
                "reads_2",
                self._lib_type_dir + "/qscores/qual_stats_12_log.txt",
                self._lib_type_dir + "/qscores/qscores_2.tsv",
                os.path.relpath(self._barcodes.name)
            )
            job.name = "generate_qscore_reads2"
            job.subname = "qscore_sheet"
            jobs.append(job)
        
        # Generate qscores graph for paired end reads
        if(self._lib_type == "ex"):
            job = rrna_amplicons.generate_qscore_graph_paired(
                self._lib_type_dir + "/qscores/qscores_1.tsv",
                self._lib_type_dir + "/qscores/qscores_2.tsv",
                self._lib_type_dir + "/qscores/qual_stats_unassembled.pdf"
            )
            job.name = "generate_qscore_graph_paired"
            job.subname = "qscore_plot"
            jobs.append(job)
            
        # Generate qscore graph for single end reads
        if(self._lib_type == "nc1" or self._lib_type == "nc2"):

            if(self._lib_type == "nc1"):
                qscore_infile = self._lib_type_dir + "/qscores/qscores_1.tsv"
            if(self._lib_type == "nc2"):
                qscore_infile = self._lib_type_dir + "/qscores/qscores_2.tsv"
            
            job = rrna_amplicons.generate_qscore_graph_single(
                qscore_infile,
                "qual_stats_unassembled",
                self._lib_type_dir + "/qual_stats_unassembled.pdf"
            )
            job.name = "generate_qscore_graph_single"
            job.subname = "qscore_plot"
            jobs.append(job)

        return jobs

    def qc(self):
        jobs = []

        clustering_method = config.param('clustering', 'clustering_method', 1, 'string')
        
        ##################################
        # For overlapping paired end reads
        # 
        if(self._lib_type == "ex" and clustering_method != "dada2"):
            
            # cut reads 1
            job = rrna_amplicons.cut_reads(
                self._lib_type_dir + "/fastqs/ncontam_nphix_1.fastq",
                config.param('cut_reads', 'R1_start', 1, 'int'),
                config.param('cut_reads', 'R1_end', 1, 'int'),
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed_1.fastq"
            )
            job.name = "cut_reads1"
            jobs.append(job)
    
            # cut reads 2
            job = rrna_amplicons.cut_reads(
                self._lib_type_dir + "/fastqs/ncontam_nphix_2.fastq",
                config.param('cut_reads', 'R2_start', 1, 'int'),
                config.param('cut_reads', 'R2_end', 1, 'int'),
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed_2.fastq"
            )
            job.name = "cut_reads2"
            jobs.append(job)
        
            job = rrna_amplicons.tags_qc_ex(
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed_1.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed_2.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed_1_Nfilt.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed_2_Nfilt.fastq",
                config.param('DEFAULT', 'primer_file', 1, 'filepath'),
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed_1_Nfilt_noprimer.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed_2_Nfilt_noprimer.fastq",
                self._lib_type_dir + "/fastqs/primer_removal_log.txt"
            )
            job.name = "tags_qc_ex"
            job.subname = "separate_readss_qc"
            jobs.append(job)
            
            job = rrna_amplicons.flash(
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed_1_Nfilt_noprimer.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed_2_Nfilt_noprimer.fastq",
                "ncontam_nphix_trimmed",
                self._lib_type_dir + "/fastqs"
            )
            job.name = "flash"
            job.subname = "flash"
            jobs.append(job)
            
            job = rrna_amplicons.tags_qc_ex_merged(
                self._lib_type_dir + "/fastqs/assembly_complete/ncontam_nphix_trimmed.extendedFrags.fastq",
                self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed.extendedFrags_QCpassed.fastq",
                self._lib_type_dir + "/fastqs/merged_qc_log.txt"
            )
            job.name = "tags_qc_ex_merged"
            job.subname = "merged_reads_qc"
            jobs.append(job)
            
            qc_passed_reads = self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed.extendedFrags_QCpassed.fastq";
            
            # generate qscores
            job = rrna_amplicons.generate_qscore_sheet(
                qc_passed_reads,
                "assembled_filtered",
                self._lib_type_dir + "/qscores/barcodes_assembled_QCpassed_log.txt",
                self._lib_type_dir + "/qscores/qscores_assembled_QCpassed.tsv",
                os.path.relpath(self._barcodes.name)
            )
            job.name = "qscore_sheet_assembled_qced"
            job.subname = "qscore_sheet"
            jobs.append(job)
            
            # generate qscore graphs
            job = rrna_amplicons.generate_qscore_graph_single(
                self._lib_type_dir + "/qscores/qscores_assembled_QCpassed.tsv",
                "qual_stats_assembled_filtered",
                self._lib_type_dir + "/qscores/qscores_assembled_QCpassed.pdf"
            )
            job.name = "qscore_graph_assembled_qced"
            job.subname = "qscore_plot"
            jobs.append(job)
        
        ##################################
        # For single end reads
        #     
        elif(self._lib_type == "nc1" or self._lib_type == "nc2" or self._lib_type == "nc1nc2"):

            if(self._lib_type == "nc1"):
                curr_infile =  self._lib_type_dir + "/fastqs/ncontam_nphix_1.fastq"
            elif(self._lib_type == "nc2"):
                curr_infile =  self._lib_type_dir + "/fastqs/ncontam_nphix_2.fastq"
            elif(self._lib_type == "nc1nc2"):
                curr_infile_R1 =  self._lib_type_dir + "/fastqs/ncontam_nphix_1.fastq"
                curr_infile_R2 =  self._lib_type_dir + "/fastqs/ncontam_nphix_2.fastq"

            # cut reads (reads1 or reads2)
            param_start = config.param('cut_reads', 'R1_start', 1, 'int') if self._lib_type == "nc1" else config.param('cut_reads', 'R2_start', 1, 'int')
            param_end = config.param('cut_reads', 'R1_end', 1, 'int') if self._lib_type == "nc1" else config.param('cut_reads', 'R2_end', 1, 'int')
            
            if(self._lib_type == "nc1"):
                curr_reads_output = self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed.fastq"
            elif(self._lib_type == "nc2"):
                curr_reads_output = self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed.fastq"
           
            if(self._lib_type == "nc1" or self._lib_type == "nc2"):
                job = rrna_amplicons.cut_reads(
                    curr_infile,
                    param_start,
                    param_end,
                    curr_reads_output
                )
                job.name = "cut_reads_single_end_reads"
                jobs.append(job)
            
            elif(self._lib_type == "nc1nc2"):
                job = rrna_amplicons.cut_reads(
                    self._lib_type_dir + "/fastqs/ncontam_nphix_1.fastq",
                    config.param('cut_reads', 'R1_start', 'int'),
                    config.param('cut_reads', 'R1_end', 'int'),
                    self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed.fastq"
                )
                job.name = "cut_reads_single_end_fwd_reads"
                jobs.append(job)
                
                job = rrna_amplicons.cut_reads(
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2.fastq",
                    config.param('cut_reads', 'R2_start', 'int'),
                    config.param('cut_reads', 'R2_end', 'int'),
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed.fastq"
                )
                job.name = "cut_reads_single_end_rev_reads"
                jobs.append(job)
     
            # tagsQC (including removing primers)
            if(self._lib_type == "nc1"):
               
                if(config.param('clustering', 'reads_type', 1, 'string') == "long_reads"):
                    job = rrna_amplicons.tags_qc_se(
                        self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed.fastq",
                        self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_Nfilt.fastq",
                        config.param('DEFAULT', 'primer_file', 1, 'filepath'),
                        self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq",
                        self._lib_type_dir + "/fastqs/qc_log.txt",
                    )
                    job.name = "tags_qc_single_end_reads"
                    job.name = "tags_QC"
                    jobs.append(job)
    
                else:# short reads
                    job = rrna_amplicons.tags_qc_se(
                        self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed.fastq",
                        self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_Nfilt.fastq",
                        config.param('DEFAULT', 'primer_file', 1, 'filepath'),
                        self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq",
                        self._lib_type_dir + "/fastqs/qc_log.txt",
                    )
                    job.name = "tags_qc_single_end_reads"
                    job.name = "tags_QC"
                    jobs.append(job)

                qc_passed_reads = self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq";
                
                # generate qscores
                job = rrna_amplicons.generate_qscore_sheet(
                    qc_passed_reads,
                    "reads_1_QCed",
                    self._lib_type_dir + "/fastqs/barcodes_log_reads_1_QCpassed.txt",
                    self._lib_type_dir + "/fastqs/qscores_reads_1_QCpassed.tsv",
                    os.path.relpath(self._barcodes.name)
                )
                job.name = "qscore_sheet"
                job.subname = "qscore_sheet"
                jobs.append(job)
       
                # generate qscore graphs
                job = rrna_amplicons.generate_qscore_graph_single(
                    self._lib_type_dir + "/fastqs/qscores_reads_1_QCpassed.tsv",
                    "qual_stats_reads_1_QCed",
                    self._lib_type_dir + "/qscores_reads_1_QCpassed.pdf"
                )
                job.name = "qscore_plot"
                job.subname = "qscore_plot"
                jobs.append(job)

            elif(self._lib_type == "nc2"):
                job = rrna_amplicons.tags_qc_se(
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed.fastq",
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_Nfilt.fastq",
                    config.param('DEFAULT', 'primer_file', 1, 'filepath'),
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_QCpassed.fastq",
                    self._lib_type_dir + "/fastqs/qc_log.txt",
                )
                job.name = "tags_qc_single_end_reads_R2"
                job.name = "tags_QC"
                jobs.append(job)
        
                qc_passed_reads = self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_QCpassed.fastq";
       
                # generate qscores
                job = rrna_amplicons.generate_qscore_sheet(
                    qc_passed_reads,
                    "reads_2_QCed",
                    self._lib_type_dir + "/fastqs/barcodes_log_reads_2_QCpassed.txt",
                    self._lib_type_dir + "/fastqs/qscores_reads_2_QCpassed.tsv",
                    os.path.relpath(self._barcodes.name)
                )
                job.name = "qscore_sheet_R2"
                job.subname = "qscore_sheet"
                jobs.append(job)
       
                # generate qscore graphs
                job = rrna_amplicons.generate_qscore_graph_single(
                    self._lib_type_dir + "/fastqs/qscores_reads_2_QCpassed.tsv",
                    "qual_stats_reads_2_QCed",
                    self._lib_type_dir + "/qscores_reads_2_QCpassed.pdf"
                )
                job.name = "qscore_plot_R2"
                job.subname = "qscore_plot"
                jobs.append(job)
                  
            # Essentially for dada2.
            elif(self._lib_type == "nc1nc2" and clustering_method == "dada2"):
                
                job = rrna_amplicons.tags_qc_R1R2(
                    self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed.fastq",
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed.fastq",
                    self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmednoN.fastq",
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmednoN.fastq",
                    config.param('DEFAULT', 'primer_file', 1, 'filepath'),
                    self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq",
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_QCpassed.fastq",
                    self._lib_type_dir + "/fastqs/qc_log.txt"
                )
                job.name = "tags_qc_paired_end_nc1_nc2"
                job.subname = "tags_QC"
                jobs.append(job)
            
                # generate qscores
                job = rrna_amplicons.generate_qscore_sheet(
                    self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq",
                    "reads_1_QCed",
                    self._lib_type_dir + "/fastqs/barcodes_log_reads_1_QCpassed.txt",
                    self._lib_type_dir + "/fastqs/qscores_reads_1_QCpassed.tsv",
                    os.path.relpath(self._barcodes.name)
                )
                job.name = "qscore_sheet"
                job.subname = "qscore_sheet"
                jobs.append(job)
           
                # generate qscore graphs
                job = rrna_amplicons.generate_qscore_graph_single(
                    self._lib_type_dir + "/fastqs/qscores_reads_1_QCpassed.tsv",
                    "qual_stats_reads_1_QCed",
                    self._lib_type_dir + "/qscores_reads_1_QCpassed.pdf"
                )
                job.name = "qscore_plot"
                job.subname = "qscore_plot"
                jobs.append(job)
            
                # generate qscores
                job = rrna_amplicons.generate_qscore_sheet(
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_QCpassed.fastq",
                    "reads_2_QCed",
                    self._lib_type_dir + "/fastqs/barcodes_log_reads_2_QCpassed.txt",
                    self._lib_type_dir + "/fastqs/qscores_reads_2_QCpassed.tsv",
                    os.path.relpath(self._barcodes.name)
                )
                job.name = "qscore_sheet"
                job.subname = "qscore_sheet"
                jobs.append(job)
           
                # generate qscore graphs
                job = rrna_amplicons.generate_qscore_graph_single(
                    self._lib_type_dir + "/fastqs/qscores_reads_2_QCpassed.tsv",
                    "qual_stats_reads_2_QCed",
                    self._lib_type_dir + "/qscores_reads_2_QCpassed.pdf"
                )
                job.name = "qscore_plot"
                job.subname = "qscore_plot"
                jobs.append(job)
    
            else:
                sys.stderr.write("No acceptable option combination for library type and clustering method...\n")
                sys.stderr.write("Valid parameters : lib_type=ex;     clustering_method=dnaclust\n")
                sys.stderr.write("Valid parameters : lib_type=ex;     clustering_method=vsearch\n")
                sys.stderr.write("Valid parameters : lib_type=ex;     clustering_method=deblur\n")
                sys.stderr.write("Valid parameters : lib_type=nc1nc2; clustering_method=dada2\n")
                sys.stderr.write("Valid parameters : lib_type=nc1;    clustering_method=dnaclust\n")
                sys.stderr.write("Valid parameters : lib_type=nc1;    clustering_method=vsearch\n")
                sys.stderr.write("Valid parameters : lib_type=nc1;    clustering_method=deblur\n")
                sys.stderr.write("Valid parameters : lib_type=nc1;    clustering_method=dada2\n")
                exit(1);

        return jobs

    def generate_clusters(self):
        jobs = []
        
        clustering_method = config.param('clustering', 'clustering_method', 1, 'string')
        
        if(clustering_method == "dnaclust"):
            job = rrna_amplicons.clustering_dnaclust(
                self._qc_passed_reads,
                os.path.relpath(self._barcodes.name),
                self._lib_type_dir + "/obs",
                config.param('clustering', 'reads_type', 1, 'string')
            )
            job.name = "generate_clusters_dnaclust"
            job.subname = "clustering"
            jobs.append(job)

        elif(clustering_method == "vsearch"):
            job = rrna_amplicons.clustering_vsearch(
                self._qc_passed_reads,
                os.path.relpath(self._barcodes.name),
                self._lib_type_dir + "/obs",
                config.param('clustering', 'reads_type', 1, 'string')
            )
            job.name = "generate_clusters_vsearch"
            job.subname = "clustering"
            jobs.append(job)
        
        elif(clustering_method == "deblur"):

            # First demultiplex assembled (FLASH) fastqs:
            job = rrna_amplicons.split_barcodes_dir(
                self._qc_passed_reads,
                os.path.relpath(self._barcodes.name),
                os.path.join(self._lib_type_dir, "obs", "split_barcodes.done"),
                os.path.join(self._lib_type_dir, "obs", "split_barcodes.txt"),
                os.path.join(self._lib_type_dir, "obs", "demultiplexed")
            )
            job.name = "split_barcodes_deblur"
            job.subname = "barcodes"
            jobs.append(job)
    
            # Do the clustering (well variants in that case).
            job = rrna_amplicons.clustering_deblur(
                os.path.join(self._lib_type_dir, "obs", "split_barcodes.done"),
                os.path.join(self._lib_type_dir, "obs", "demultiplexed"),
                os.path.join(self._lib_type_dir, "obs", "deblur"),
                os.path.join(self._lib_type_dir, "obs", "deblur", "all.biom"),
                os.path.join(self._lib_type_dir, "obs", "deblur", "all.seqs.fa"),
                os.path.join(self._lib_type_dir, "obs", "deblur", "deblur.done")
            )
            job.name = "generate_asv_deblur"
            job.subname = "clustering"
            jobs.append(job)
            
            job = rrna_amplicons.convert_biom_to_tsv(
                os.path.join(self._lib_type_dir, "obs", "deblur", "all.biom"),
                os.path.join(self._lib_type_dir, "obs", "deblur", "all.tsv")
            )
            job.name = "convert_biom_to_tsv_deblur"
            jobs.append(job)
            
            job = rrna_amplicons.postprocess_asvs(
                os.path.join(self._lib_type_dir, "obs", "deblur", "deblur.done"),
                os.path.join(self._lib_type_dir, "obs", "deblur", "all.tsv"),
                os.path.join(self._lib_type_dir, "obs", "all.tsv"),
                os.path.join(self._lib_type_dir, "obs", "all.fasta")
            )
            job.name = "postprocess_deblur"
            jobs.append(job)
            
            job = rrna_amplicons.remove_chimeras(
                os.path.join(self._lib_type_dir, "obs", "all.tsv"),
                os.path.join(self._lib_type_dir, "obs", "all.fasta"),
                os.path.join(self._lib_type_dir, "obs"),
                os.path.join(self._lib_type_dir, "obs", "obs.tsv"),
                os.path.join(self._lib_type_dir, "obs", "obs.fasta")
                
            )
            job.name = "remove_chimeras"
            job.subname = "clustering"
            jobs.append(job)

        elif(clustering_method == "dada2"):
            # dada2 can't process dada2 as merged (assembled) sequences.
            barcodes_done_files = []
            if(self._lib_type == "ex"):
                sys.stderr.write("Can't use lib_type=ex with dada2. Please change lib_type to lib_type=nc1nc2\nStrickly speaking, DADA2 takes in input fwd and rev reads separately...\n")
                exit(1)
            if(self._lib_type == "nc1nc2"):

                job = rrna_amplicons.split_barcodes_dir(
                    self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq",
                    os.path.relpath(self._barcodes.name),
                    os.path.join(self._lib_type_dir, "obs", "split_barcodes_R1.done"),
                    os.path.join(self._lib_type_dir, "obs", "split_barcodes_R1.txt"),
                    os.path.join(self._lib_type_dir, "obs", "demultiplexed"),
                    "_R1",
                    True
                )
                job.name = "split_barcodes_for_dada2_R1"
                job.subname = "barcodes"
                jobs.append(job)

                barcodes_done_files.append(os.path.join(self._lib_type_dir, "obs", "split_barcodes_R1.done"))
            
                job = rrna_amplicons.split_barcodes_dir(
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_QCpassed.fastq",
                    os.path.relpath(self._barcodes.name),
                    os.path.join(self._lib_type_dir, "obs", "split_barcodes_R2.done"),
                    os.path.join(self._lib_type_dir, "obs", "split_barcodes_R2.txt"),
                    os.path.join(self._lib_type_dir, "obs", "demultiplexed"),
                    "_R2",
                    True
                )
                job.name = "split_barcodes_for_dada2_R2"
                job.subname = "barcodes"
                jobs.append(job)
                
                barcodes_done_files.append(os.path.join(self._lib_type_dir, "obs", "split_barcodes_R2.done"))
            
            elif(self._lib_type == "nc1"):
                
                job = rrna_amplicons.split_barcodes_dir(
                    self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq",
                    os.path.relpath(self._barcodes.name),
                    os.path.join(self._lib_type_dir, "obs", "split_barcodes_R1.done"),
                    os.path.join(self._lib_type_dir, "obs", "split_barcodes_R1.txt"),
                    os.path.join(self._lib_type_dir, "obs", "demultiplexed"),
                    "_R1"
                )
                job.name = "split_barcodes_for_dada2_nc1"
                job.subname = "barcodes"
                jobs.append(job)
                
                barcodes_done_files.append(os.path.join(self._lib_type_dir, "obs", "split_barcodes_R1.done"))
            
            elif(self._lib_type == "nc2"):
                
                job = rrna_amplicons.split_barcodes_dir(
                    self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_QCpassed.fastq",
                    os.path.relpath(self._barcodes.name),
                    os.path.join(self._lib_type_dir, "obs", "split_barcodes_R2.done"),
                    os.path.join(self._lib_type_dir, "obs", "split_barcodes_R2.txt"),
                    os.path.join(self._lib_type_dir, "obs", "demultiplexed"),
                    "_R2"
                )
                job.name = "split_barcodes_for_dada2_nc2"
                job.subname = "barcodes"
                jobs.append(job)
                
                barcodes_done_files.append(os.path.join(self._lib_type_dir, "obs", "split_barcodes_R2.done"))
            
            # Do the clustering (well variants in that case.
            # Determine dada_type
            data_type = ""
            if(config.param('clustering', 'reads_type', 1, 'string') == "long_reads" and self._lib_type == "nc1"):
                data_type = "lse"
            elif(config.param('clustering', 'reads_type', 1, 'string') == "short_reads" and self._lib_type == "nc1nc2"):
                data_type = "spe"
            elif(config.param('clustering', 'reads_type', 1, 'string') == "short_reads" and self._lib_type == "nc1"):
                data_type = "sse"
            elif(config.param('clustering', 'reads_type', 1, 'string') == "short_reads" and self._lib_type == "nc2"):
                data_type = "sse_nc2"
            else:
                sys.stderr.write("Supported methods: 'long and single-end reads', 'short and paired-end reads' or 'short and single-end reads'\n")
                
            job = rrna_amplicons.clustering_dada2(
                barcodes_done_files,
                os.path.join(self._lib_type_dir, "obs", "demultiplexed"),
                os.path.join(self._lib_type_dir, "obs", "dada2"),
                os.path.join(self._lib_type_dir, "obs", "dada2", "all_nochimera.tsv"),
                os.path.join(self._lib_type_dir, "obs", "dada2", "dada2.done"),
                data_type
            )
            job.name = "generate_asv_dada2"
            job.subname = "clustering"
            jobs.append(job)
            
            job = rrna_amplicons.postprocess_asvs(
                os.path.join(self._lib_type_dir, "obs", "dada2", "dada2.done"),
                os.path.join(self._lib_type_dir, "obs", "dada2", "all_nochimera.tsv"),
                os.path.join(self._lib_type_dir, "obs", "all_nochimera.tsv"),
                os.path.join(self._lib_type_dir, "obs", "all_nochimera.fasta")
                
            )
            job.name = "postprocess_dada2"
            #job.subname = "clustering"
            jobs.append(job)
            
            job = rrna_amplicons.remove_chimeras(
                os.path.join(self._lib_type_dir, "obs", "all_nochimera.tsv"),
                os.path.join(self._lib_type_dir, "obs", "all_nochimera.fasta"),
                os.path.join(self._lib_type_dir, "obs"),
                os.path.join(self._lib_type_dir, "obs", "obs.tsv"),
                os.path.join(self._lib_type_dir, "obs", "obs.fasta"),
                "ref_only"
                
            )
            job.name = "remove_chimeras_dada2"
            job.subname = "clustering"
            jobs.append(job)
        else:
            sys.stderr.write("Supported methods: 'dnaclust', 'vsearch', 'dada2' or 'deblur'\n")

        return jobs

    def classify(self): # Add PacBio's blast here!
        # classify with RDP
        jobs = []
        job = microbial_ecology.rdp_wrapper(
            self._lib_type_dir + "/obs/obs.fasta",
            self._lib_type_dir + "/rdp/rdp.tsv"
        )
        job.name = "classify"
        job.subname = "RDP"
        jobs.append(job)
 
        return jobs
    
    def generate_feature_table(self):
        jobs = []
            
        #sys.stderr.write(os.path.join(self._lib_type_dir, "feature_tables", "feature_table.tsv") + "\n")

        job = microbial_ecology.generate_feature_table(
            os.path.join(self._lib_type_dir, "obs", "obs.tsv"),
            os.path.join(self._lib_type_dir, "rdp", "rdp.tsv"),
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table.tsv"),
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table_failed.tsv")
        )
        job.name = "generate_feature_table_open"
        job.subname = "add_taxonomy"
        jobs.append(job)
         
        return jobs

    ##################################
    # remove empty columns from otu table + split otu tables by organism types + convert otu tables from tsv to biom format.
    #
    def filter_feature_table(self):
        jobs = []
        job = microbial_ecology.filter_feature_table(
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table.tsv"),
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered.tsv")
        )
        job.name = "filter_feature_table"
        job.subname = "filter_feature_table"
        jobs.append(job)
        
        job = microbial_ecology.split_feature_table(
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered.tsv"),
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + ".tsv"),
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_others.tsv"),
            self._organism_type
        )
        job.name = "split_feature_table_" + self._organism_type
        job.name = "split_feature_table"
        jobs.append(job)
            
        # Rarefy tables to the specified arbitrary thresholds from the .ini file.
        rarefaction_values_string = config.param('rarefaction', 'n_normalization', 1, 'string')
        rarefaction_values = rarefaction_values_string.split(':')
        
        for rarefaction_value in rarefaction_values:
            if not os.path.exists(os.path.join(self._lib_type_dir, "feature_tables", "rarefactions", str(rarefaction_value))):
                os.makedirs(os.path.join(self._lib_type_dir, "feature_tables", "rarefactions", str(rarefaction_value)))
             
            job = microbial_ecology.rtk_multi_even(
                os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + ".tsv"),
                os.path.join(self._lib_type_dir, "feature_tables", "rarefactions", str(rarefaction_value), "rarefaction"),
                "#FEATURE_ID",
                remove_last_col = True
            )
            job.name = "rarefy_multi_norm_" + rarefaction_value
            job.subname = "rarefaction"
            jobs.append(job)

            job = microbial_ecology.merge_rarefied_tables(
                os.path.join(self._lib_type_dir, "feature_tables", "rarefactions", str(rarefaction_value), "rarefaction"),
                os.path.join(self._lib_type_dir, "feature_tables", "feature_table.tsv"),
                os.path.join(self._lib_type_dir, "feature_tables", "rarefactions", str(rarefaction_value), "rarefaction", "rtk.done"),
                os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + "_rarefied_" + str(rarefaction_value) + ".tsv")
            )
            job.name = "merge_rarefied_tables_norm_" + rarefaction_value
            job.subname = "merge_rarefied_tables"
            jobs.append(job)
        
        # Then, normalize using edgeR with CPMs.
        job = microbial_ecology.normalize_feature_table(
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + ".tsv"),
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + "_normalized.tsv")
        )
        job.name = "normalize_feature_table_target_organisms"
        job.subname = "normalization"
        jobs.append(job)
        
        # Convert tsv to biom
        #job = microbial_ecology.convert_tsv_to_biom_hdf5(
        #    os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + ".tsv"),
        #    os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + ".biom")
        #)
        #job.name = "convert_feature_table_to_biom_target_organisms"
        #job.subname = "convert_feature_table"
        #jobs.append(job)
            
        return jobs

    # AMPLICONTAGGER_DIVERSITY
    def prepare_tables_for_diversity(self):
       
        # Here set a common name for OTU table that will be used in downstream steps - makes life easier.
        feature_table_rarefied   = os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + "_rarefied_" + self._rarefaction_value + ".tsv")
        feature_table_normalized = os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + "_normalized.tsv")
        
        jobs = []
        
        # Remove empty columns in otu table. It rarely happens, but just to be sure.
        job = microbial_ecology.filter_and_sort_feature_table(
            feature_table_normalized,
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table_final_normalized.tsv")
        )
        job.name = "filter_feature_table_normalized"
        job.subname = "filter_feature_table"
        jobs.append(job)
        
        job = microbial_ecology.filter_and_sort_feature_table(
            feature_table_rarefied,
            os.path.join(self._lib_type_dir, "feature_tables", "feature_table_final_rarefied.tsv")
        )
        job.name = "filter_feature_table_rarefied"
        job.subname = "filter_feature_table"
        jobs.append(job)
           
        return jobs

    def align(self):
        
        curr_feature_table = os.path.join(self._lib_type_dir, "feature_tables", "feature_table_final_rarefied.tsv")
        
        jobs = []

        job = microbial_ecology.filter_for_mafft(
            curr_feature_table,
            os.path.join(self._lib_type_dir, "obs", "obs.fasta"),
            os.path.join(self._lib_type_dir, "feature_tables", "features.fasta")
        )
        job.name = "filter_for_mafft"
        job.subname = "filter_for_mafft"
        jobs.append(job)
        
        job = microbial_ecology.mafft(
            os.path.join(self._lib_type_dir, "feature_tables", "features.fasta"),
            os.path.join(self._lib_type_dir, "mafft", "mafft.fasta")
        )
        job.name = "mafft"
        job.subname = "mafft"
        jobs.append(job)
       
        job = microbial_ecology.fasttree(
          os.path.join(self._lib_type_dir, "mafft", "mafft.fasta"),
          os.path.join(self._lib_type_dir, "fasttree", "tree.fasttree")
        );
        job.name = "fasttree"
        job.subname = "fasttree"
        jobs.append(job)
          
        return jobs
        
    def summarize_taxonomy(self):
       
        jobs = []

        curr_feature_table = os.path.join(self._lib_type_dir, "feature_tables", "feature_table_final_rarefied.tsv")
        
        plot_taxa_string_abs = ""
        plot_taxa_string_rel = ""
        feature_table_prefix = "feature_table_final_rarefied"
        
        for i in range(1, config.param('summarize_taxonomy', 'taxonomy_depth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_absolute(
                curr_feature_table,
                i,
                os.path.join(self._lib_type_dir, "tax_summary", "absolute")
            )
            job.name = "tax_summary_absolute_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
            plot_taxa_string_abs += self._lib_type_dir + "/tax_summary/absolute/" + feature_table_prefix + "_L" + str(i) + ".tsv,"
        
        # Taxonomy relative raw
        for i in range(1, config.param('summarize_taxonomy', 'taxonomy_depth', 1, 'int')):
            job = microbial_ecology.summarize_taxonomy_relative(
                curr_feature_table,
                i,
                os.path.join(self._lib_type_dir, "tax_summary", "relative")
            )
            job.name = "tax_summary_relative_L" + str(i)
            job.subname = "summarize_taxonomy"
            jobs.append(job)
            plot_taxa_string_rel += self._lib_type_dir + "/tax_summary/relative/" + feature_table_prefix + "_L" + str(i) + ".tsv,"
        
        plot_taxa_string_abs = plot_taxa_string_abs[:-1]
        plot_taxa_string_rel = plot_taxa_string_rel[:-1]
        
        # draw taxonomy plots with custom R script. 
        for i in range(1, config.param('summarize_taxonomy', 'taxonomy_depth', 1, 'int')):
            job = microbial_ecology.plot_taxa_single(
                os.path.join(self._lib_type_dir, "tax_summary", "relative", feature_table_prefix + "_L" + str(i) + ".tsv"),
                os.path.join(self._lib_type_dir, "tax_summary", "relative", "plots"),    
                "taxonomy_L" + str(i) + "_relative"
            )
            job.name = "plot_taxa_single_relative_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
            
            job = microbial_ecology.plot_taxa_single(
                os.path.join(self._lib_type_dir, "tax_summary", "absolute", feature_table_prefix + "_L" + str(i) + ".tsv"),
                os.path.join(self._lib_type_dir, "tax_summary", "absolute", "plots"),    
                "taxonomy_L" + str(i) + "_absolute"
            )
            job.name = "plot_taxa_single_absolute_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
            
            # plot with mapping file.
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._lib_type_dir, "tax_summary", "relative", feature_table_prefix + "_L" + str(i) + ".tsv"),
                os.path.join(self._lib_type_dir, "tax_summary", "relative", "plots"),    
                "taxonomy_L" + str(i) + "_relative"
            )
            job.name = "plot_taxa_multi_relative_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
            
            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                os.path.join(self._lib_type_dir, "tax_summary", "absolute", feature_table_prefix + "_L" + str(i) + ".tsv"),
                os.path.join(self._lib_type_dir, "tax_summary", "absolute", "plots"),    
                "taxonomy_L" + str(i) + "_absolute"
            )
            job.name = "plot_taxa_multi_absolute_L" + str(i)
            job.subname = "plot_taxa"
            jobs.append(job)
        
        return jobs

    def beta_diversity(self):
        jobs = []
        
        curr_feature_table = os.path.join(self._lib_type_dir, "feature_tables", "feature_table_final_rarefied.tsv")
        feature_table_prefix = "feature_table_final_rarefied"
        
        job = microbial_ecology.beta_diversity(
            curr_feature_table,
            "weighted-unifrac",
            os.path.join(self._lib_type_dir, "beta_div"),
            os.path.join(self._lib_type_dir, "fasttree", "tree.fasttree"),
            "weighted_unifrac_" + feature_table_prefix
        )
        job.name = "beta_diversity_unifrac_weighted"
        job.subname = "beta_diversity"
        jobs.append(job)
        
        job = microbial_ecology.beta_diversity(
            curr_feature_table,
            "bray-curtis",
            os.path.join(self._lib_type_dir, "beta_div"),
            os.path.join(self._lib_type_dir, "fasttree", "tree.fasttree"),
            "bray_curtis_" + feature_table_prefix
        )
        job.name = "beta_diversity_bray_curtis"
        job.subname = "beta_diversity"
        jobs.append(job)

        job = microbial_ecology.principal_coordinates(
            os.path.join(self._lib_type_dir, "beta_div", "weighted_unifrac_" + feature_table_prefix + ".tsv"),
            os.path.join(self._lib_type_dir, "beta_div", "weighted_unifrac_" + feature_table_prefix + "_coords.tsv")
        )
        job.name = "beta_diversity_unifrac_weighted_PCoA"
        job.subname = "pca"
        jobs.append(job)
        
        job = microbial_ecology.principal_coordinates(
            os.path.join(self._lib_type_dir, "beta_div", "bray_curtis_" + feature_table_prefix + ".tsv"),
            os.path.join(self._lib_type_dir, "beta_div", "bray_curtis_" + feature_table_prefix + "_coords.tsv")
        )
        job.name = "beta_diversity_bray_curtis_PCoA"
        job.subname = "pca"
        jobs.append(job)

        job = microbial_ecology.pca_3d_plot(
            os.path.join(self._lib_type_dir, "beta_div", "weighted_unifrac_" + feature_table_prefix + "_coords.tsv"),
            os.path.join(self._lib_type_dir, "beta_div", "plots", "3d_weighted_unifrac")
        )
        job.name = "beta_diversity_weighted_unifrac_3d_plot"
        job.subname = "pca_plot"
        jobs.append(job)

        job = microbial_ecology.pca_3d_plot(
            os.path.join(self._lib_type_dir, "beta_div", "bray_curtis_" + feature_table_prefix + "_coords.tsv"),
            os.path.join(self._lib_type_dir, "beta_div", "plots", "3d_bray_curtis")
        )
        job.name = "beta_diversity_bray_curtis_3d_plot"
        job.subname = "pca_plot"
        jobs.append(job)
        
        job = microbial_ecology.plot_pca(
            os.path.join(self._lib_type_dir, "beta_div", "weighted_unifrac_" + feature_table_prefix + "_coords.tsv"),
            os.path.join(self._lib_type_dir, "beta_div", "plots"),
            "pca_weighted_unifrac"
        )
        job.name = "weighted_unifrac_pcoa_plots"
        job.subname = "pca_plot"
        jobs.append(job)
        
        job = microbial_ecology.plot_pca(
            os.path.join(self._lib_type_dir, "beta_div", "bray_curtis_" + feature_table_prefix + "_coords.tsv"),
            os.path.join(self._lib_type_dir, "beta_div", "plots"),
            "pca_bray_curtis"
        )
        job.name = "bray_curtis_pcoa_plots"
        job.subname = "pca_plot"
        jobs.append(job) 

        return jobs

    def alpha_diversity(self):
        jobs = []
        
        curr_feature_table_tsv = os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + ".tsv")
     
        # Rarefaction curves using all data - get an idea if sequencing depth is appropriate for each sample.
        job = microbial_ecology.rtk(
            curr_feature_table_tsv,
            os.path.join(self._lib_type_dir, "alpha_div"),
            "#FEATURE_ID",
            remove_last_col = True
        )
        job.name = "alpha_diversity_saturation"
        job.subname = "alpha_diversity"
        jobs.append(job)
        
        job = microbial_ecology.rarefaction_plots(
            os.path.join(self._lib_type_dir, "alpha_div"),
            os.path.join(self._lib_type_dir, "alpha_div", "plots")
        )
        job.name = "rarefaction_plots"
        job.subname = "rarefaction_plots"
        jobs.append(job)
        
        return jobs

    def heatmap(self):
        jobs = []

        range_string = config.param('otu_heatmap', 'n', 1, 'string')
        range = range_string.split(':')
        for n in range:

            job = microbial_ecology.otu_heatmap(
                os.path.join(self._lib_type_dir, "feature_tables", "feature_table_final_rarefied.tsv"),
                os.path.join(self._lib_type_dir, "heatmap"),
                "OTU_heatmap",
                n
            )
            job.name = "otu_heatmap_" + n
            job.subname = "otu_heatmap"
            jobs.append(job)

        return jobs
    
    def blast(self):
        jobs = []
        job = microbial_ecology.blast(
            os.path.join(self._lib_type_dir, "feature_tables", "features.fasta"),
            os.path.join(self._lib_type_dir, "blast", "blast.out"),
            config.param('blast', 'db')
        )
        job.name = "blast"
        job.subname = "blast"
        jobs.append(job)

        return jobs

    def statistics(self):
        jobs = []
        
        os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + ".tsv")
        
        # Compute GLM diff. expr. of OTUs following provided mapping file.
        if( config.param('DOA_edger', 'do_glm', 1, 'string').lower() == "yes".lower() ):
            mapping_files = config.param('DOA_edger', 'mapping_files', 1, 'string')
            mapping_files = mapping_files.split(":")
            for mapping_file in mapping_files:
                job = microbial_ecology.edger_feature_glm(
                    os.path.join(self._lib_type_dir, "feature_tables", "feature_table_filtered_" + self._organism_type + ".tsv"),
                    os.path.join(self._lib_type_dir, "DOA","GLM"),
                    mapping_file
                )
                job.name = "DOA_edger_" + mapping_file
                job.subname = "DOA_edger"
                jobs.append(job)
        
        else:
            sys.stderr.write("Skipping GLM differential abundance of OTUs.\n")
        
        if( config.param('DOA_ancom', 'do_ancom', 1, 'string').lower() == "yes".lower() ):
            mapping_files = config.param('DOA_ancom', 'mapping_files', 1, 'string')
            mapping_files = mapping_files.split(":")
            for mapping_file in mapping_files:
                job = microbial_ecology.ancom_otus(
                    os.path.join(self._lib_type_dir, "feature_tables", "feature_table_final_rarefied.tsv"),
                    os.path.join(self._lib_type_dir, "DOA", "ancom"),
                    mapping_file
                )
                job.name = "DOA_ancom_" + mapping_file
                job.subname = "DOA_ancom"
                jobs.append(job)
        
        else:
            sys.stderr.write("Skipping ANCOM differential abundance of OTUs.\n")


        # Populate lists of files to be analyzed for reads count.
        report_files = []
        report_names = []
        
        if 'self._external_infile' in globals():
            infile = self._external_infile
        else:
            infile = self._infile

        report_files.append(infile)
        report_files.append(self._duk_dir + "/contam.fastq")
        report_files.append(self._duk_dir + "/ncontam_phix.fastq")
        report_files.append(self._duk_dir + "/ncontam_nphix.fastq")
        report_names.append("total_reads")
        report_names.append("contaminants_reads")
        report_names.append("phix_reads")
        report_names.append("non_contam_non_phix_reads")
       
        if(re.match("ex", self._lib_type)):
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_1.fastq")
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_2.fastq")
            report_files.append(self._lib_type_dir + "/fastqs/assembly_complete/ncontam_nphix_trimmed.extendedFrags.fastq")  
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed.extendedFrags_QCpassed.fastq")
            report_names.append("non_contam_non_phix_reads_1")
            report_names.append("non_contam_non_phix_reads_2")
            report_names.append("assembled_reads")  
            report_names.append("assembled_reads_QC_passed")  
        elif(re.match("nc1nc2", self._lib_type)):
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_1.fastq")
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_2.fastq")
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq")
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_QCpassed.fastq")
            report_names.append("non_contam_non_phix_reads_1")
            report_names.append("non_contam_non_phix_reads_2")
            report_names.append("reads_1_QC_passed")
            report_names.append("reads_2_QC_passed")
        elif(re.match("nc1", self._lib_type)):
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_1.fastq")
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq") 
            report_names.append("non_contam_non_phix_reads_1")
            report_names.append("reads_1_QC_passed")
        elif(re.match("nc2", self._lib_type)):
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_2.fastq")
            report_files.append(self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_QCpassed.fastq")
            report_names.append("non_contam_non_phix_reads_2")
            report_names.append("reads_2_QC_passed")
        
        ## Analysis type
        if(re.match("ex", self._lib_type)):
            analysis_type = "assembled_clustered"
            barcodes_dist = self._lib_type_dir + "/fastqs/ncontam_nphix_barcodes.txt"
            feature_table = self._lib_type_dir + "/feature_tables/feature_table.tsv"
            obs_table = self._lib_type_dir + "/obs/obs.tsv"
        elif(re.match("nc1nc2", self._lib_type)):
            analysis_type = "assembled_variants"
            barcodes_dist = self._lib_type_dir + "/fastqs/ncontam_nphix_barcodes.txt"
            feature_table = self._lib_type_dir + "/feature_tables/feature_table.tsv"
            obs_table = self._lib_type_dir + "/obs/obs.tsv"
        elif(re.match("nc1", self._lib_type)):
            analysis_type = "reads_1_clustered"
            barcodes_dist = self._lib_type_dir + "/fastqs/barcodes_log_reads_1_QCpassed.txt"
            feature_table = self._lib_type_dir + "/feature_tables/feature_table.tsv"
            obs_table = self._lib_type_dir + "/obs/obs.tsv"
        elif(re.match("nc2", self._lib_type)):
            analysis_type = "reads_2_clustered"
            barcodes_dist = self._lib_type_dir + "/fastqs/barcodes_log_reads_2_QCpassed.txt"
            feature_table = self._lib_type_dir + "/feature_tables/feature_table.tsv"
            obs_table = self._lib_type_dir + "/obs/obs.tsv"
 
        job = rrna_amplicons.count_report(
            report_files,
            report_names,
            analysis_type,
            barcodes_dist,
            feature_table,
            obs_table,
            self._lib_type_dir + "/countReport.tsv"
        )
        job.name = "count_report"
        job.subname = "count_report"
        jobs.append(job)

        return jobs

    def deliverables(self):
        jobs = []
         
        # Compute date to generate report folder.
        i = datetime.datetime.now()
        year = "%04d" %i.year
        day = "%02d" %i.day
        month = "%02d" %i.month
        prefix = year + "_" + month + "_" + day + "_" + self._project_name + "_amplicontagger_complete"

        job = microbial_ecology.generate_tarball(
            os.path.join(self._lib_type_dir, "blast", "blast.out"),
            self._config,
            self._lib_type_dir,
            self._root_dir + "/" + prefix,
            self.get_output_files
        )
        job.name = "deliverables_tarball"
        job.subname = "report"
        jobs.append(job)
        
        return jobs

    def picrust(self):
        jobs = []
        
        # assign taxonomy with blast for picrust's closed reference method.
        job = microbial_ecology.blast(
            os.path.join(self._lib_type_dir, "obs", "obs.fasta"),
            os.path.join(self._lib_type_dir, "blast", "blast_greengenes.out"),
            config.param('blast', 'db_greengenes', 1, 'filepath')
        )
        job.name = "picrust_blast_greengenes"
        job.subname = "blast"
        jobs.append(job)
  
        job = microbial_ecology.generate_feature_table_closed(
            os.path.join(self._lib_type_dir, "obs", "obs.tsv"),
            os.path.join(self._lib_type_dir, "blast", "blast_greengenes.out"),
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes.tsv"),
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_failed.tsv"),
            os.path.join(self._lib_type_dir, "picrust", "link.tsv")
        )
        job.name = "picrust_generate_feature_table_closed"
        job.subname = "add_taxonomy"
        jobs.append(job)

        # collapse duplicates
        job = microbial_ecology.collapse_duplicates(
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes.tsv"),
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups.tsv")
        )
        job.name = "picrust_collapse_duplicates"
        job.subname = "collapse_duplicates"
        jobs.append(job)
        
        job = microbial_ecology.convert_tsv_to_biom_hdf5(
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups.tsv"),
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups.biom")
        )
        job.name = "picrust_convert_tsv_to_biom_hdf5"
        job.subname = "convert_tsv_to_biom"
        jobs.append(job)
        
        job = microbial_ecology.normalize_by_copy_number(
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups.biom"),
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups_nbcn.biom"),
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups_nbcn.tsv")
        )
        job.name = "picrust_normalize_by_copy_number"
        job.subname = "normalize_by_copy_number"
        jobs.append(job)
        
        job = microbial_ecology.predict_traits(
            config.param('picrust', 'traits', 1, 'filepath'),
            config.param('picrust', 'tree', 1, 'filepath'),
            config.param('picrust', 'counts', 1, 'filepath'),
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups_nbcn.tsv"),
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups_nbcn_picrust_KO_occurence.tsv")
        )
        job.name = "predict_traits"
        job.subname = "predict_traits"
        jobs.append(job)
        
        job = microbial_ecology.predict_metagenomes(
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups_nbcn.biom"),
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups_nbcn_picrust_KO_sum.biom"),
            os.path.join(self._lib_type_dir, "picrust", "feature_table_greengenes_nodups_nbcn_picrust_KO_sum.tsv")
        )
        job.name = "predict_metagenomes"
        job.subname = "predict_metagenomes"
        jobs.append(job)

        return jobs

    def cleanup(self):
        jobs = []
        
        # Compute date to generate report folder.
        i = datetime.datetime.now()
        year = "%04d" %i.year
        day = "%02d" %i.day
        month = "%02d" %i.month
        prefix = year + "_" + month + "_" + day + "_" + self._project_name + "_complete_report"
        
        job = microbial_ecology.cleanup(
            self._root_dir,
            self._root_dir + "/" + prefix + ".tar.gz"
        )
        job.name = "cleanup"
        job.subname = "cleanup"
        jobs.append(job)
     
        return jobs

    # Override illumina.py 
    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_nanuq_readset_file(self.args.readsets.name)
        return self._readsets

    @property
    def steps(self):
        
        return [
            self.validate, #1
            self.remove_contam, #2
            self.split_barcodes, #3
            self.qscores, #4
            self.qc, #5
            self.generate_clusters, #6
            self.classify, #7
            self.generate_feature_table, #8
            self.filter_feature_table, #9
            # AmpliconTagger diversity
            self.prepare_tables_for_diversity,#10
            self.align, #11
            self.summarize_taxonomy, #12
            self.beta_diversity, #13
            self.alpha_diversity, #14
            self.heatmap, #15
            self.blast, #16
            self.statistics, #17
            self.deliverables, #18
            self.cleanup #19
        ]

    def set_local_variables(self):
        self._parser_local = self.argparser#ArgumentParser(description='Process options.')

        # for AmpliconTagger
        self._args_local = self._parser_local.parse_args()
        config.parse_files(self._args_local.config)
        self._barcodes = self._args_local.barcodes

        if self._args_local.external_infile:
            self._external_infile = self._args_local.external_infile
            self._external_infile = os.path.abspath(self._external_infile.name)

        self._config = self._args_local.config[0]
        self._config = os.path.abspath(self._config.name)
        self._root_dir = self._args_local.output_dir + "/amplicontagger"

        if not os.path.exists(self._root_dir):
            os.makedirs(self._root_dir)

        self._duk_dir = self._root_dir + "/duk"
        self._project_name = config.param('default', 'project_name', type='string')
        
        if not self._project_name:
            self._project_name = "AmpliconTagger-project" 

        self._organism_type = config.param('default', 'organism_type', type='string')
        self._lib_type = config.param('default', 'library_type', type='string')
 
        # Make directories
        lib_type = ""
        if self._lib_type == "ex":
            lib_type = "reads_12"
            self._infile = os.path.join("./", "raw_reads", "preprocessed.fastq.gz")
            self._lib_type_dir = self._root_dir + "/" + lib_type
            self._qc_passed_reads = self._lib_type_dir + "/fastqs/ncontam_nphix_trimmed.extendedFrags_QCpassed.fastq";
        elif self._lib_type == "nc1":
            lib_type = "reads_1"
            self._infile = os.path.join("./", "raw_reads", "preprocessed.fastq.gz")
            self._lib_type_dir = self._root_dir + "/" + lib_type
            self._qc_passed_reads = self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq";
        elif self._lib_type == "nc2":
            lib_type = "reads_2"
            self._infile = os.path.join("./", "raw_reads", "preprocessed.fastq.gz")
            self._lib_type_dir = self._root_dir + "/" + lib_type
            self._qc_passed_reads = self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_QCpassed.fastq";
        elif self._lib_type == "nc1nc2":
            lib_type = "reads_12"
            self._infile = os.path.join("./", "raw_reads", "preprocessed.fastq.gz")
            self._lib_type_dir = self._root_dir + "/" + lib_type
            self._qc_passed_reads_1 = self._lib_type_dir + "/fastqs/ncontam_nphix_1_trimmed_QCpassed.fastq";
            self._qc_passed_reads_2 = self._lib_type_dir + "/fastqs/ncontam_nphix_2_trimmed_QCpassed.fastq";
        #sys.stderr.write('Project root dir : ' + self._root_dir + '\n')
        self.make_directories(self._root_dir, lib_type)
        rarefaction_values_string = config.param('rarefaction', 'n_normalization', 1, 'string')
        rarefaction_values = rarefaction_values_string.split(':')
        self._rarefaction_value = rarefaction_values[len(rarefaction_values)-1]
  
    # Define and make directories. Also desing initial infile.
    def make_directories(self, root_dir, lib_type):
          
        if not os.path.exists(root_dir):
            os.makedirs(root_dir)
        if not os.path.exists(root_dir + "/duk"):
            os.makedirs(os.path.join(root_dir, "duk"))
        if not os.path.exists(root_dir + "/" + lib_type):
            os.makedirs(os.path.join(root_dir, lib_type))
        if not os.path.exists(root_dir + "/" + lib_type + "/obs"):
            os.makedirs(os.path.join(root_dir, lib_type, "obs"))
        if not os.path.exists(root_dir + "/" + lib_type + "/feature_tables"):
            os.makedirs(os.path.join(root_dir, lib_type, "feature_tables"))
        if not os.path.exists(root_dir + "/" + lib_type + "/rdp"):
            os.makedirs(os.path.join(root_dir, lib_type, "rdp"))
        if not os.path.exists(root_dir + "/" + lib_type + "/qscores"):
            os.makedirs(os.path.join(root_dir, lib_type, "qscores"))
        if not os.path.exists(root_dir + "/" + lib_type + "/blast"):
            os.makedirs(os.path.join(root_dir, lib_type, "blast"))
        if not os.path.exists(root_dir + "/" + lib_type + "/fastqs"):
            os.makedirs(os.path.join(root_dir, lib_type, "fastqs"))

        ## AmpliconTagger_diversity
        if not os.path.exists(root_dir + "/" + lib_type):
            os.makedirs(root_dir, lib_type)
        if not os.path.exists(root_dir + "/" + lib_type + "/heatmap"):
            os.makedirs(os.path.join(root_dir, lib_type, "heatmap"))
        if not os.path.exists(root_dir + "/" + lib_type + "/blast"):
            os.makedirs(os.path.join(root_dir, lib_type, "blast"))
        if not os.path.exists(root_dir + "/" + lib_type + "/feature_tables"):
            os.makedirs(os.path.join(root_dir, lib_type, "feature_tables"))
        if not os.path.exists(root_dir + "/" + lib_type + "/mafft"):
            os.makedirs(os.path.join(root_dir, lib_type, "mafft"))
        if not os.path.exists(root_dir + "/" + lib_type + "/fasttree"):
            os.makedirs(os.path.join(root_dir, lib_type, "fasttree"))
        if not os.path.exists(root_dir + "/" + lib_type + "/tax_summary"):
            os.makedirs(os.path.join(root_dir, lib_type, "tax_summary"))
        if not os.path.exists(root_dir + "/" + lib_type + "/tax_summary" + "/absolute"):
            os.makedirs(os.path.join(root_dir, lib_type, "tax_summary", "absolute"))
        if not os.path.exists(root_dir + "/" + lib_type + "/tax_summary" + "/relative"):
            os.makedirs(os.path.join(root_dir, lib_type, "tax_summary", "relative"))
        if not os.path.exists(root_dir + "/" + lib_type + "/tax_summary" + "/absolute/plots"):
            os.makedirs(os.path.join(root_dir, lib_type, "tax_summary", "absolute", "plots"))
        if not os.path.exists(root_dir + "/" + lib_type + "/tax_summary" + "/relative/plots"):
            os.makedirs(os.path.join(root_dir, lib_type, "tax_summary", "relative", "plots"))
        if not os.path.exists(root_dir + "/" + lib_type + "/beta_div"):
            os.makedirs(os.path.join(root_dir, lib_type, "beta_div"))
        if not os.path.exists(root_dir + "/" + lib_type + "/alpha_div"):
            os.makedirs(os.path.join(root_dir, lib_type, "alpha_div"))
        if not os.path.exists(root_dir + "/" + lib_type + "/alpha_div/plots"):
            os.makedirs(os.path.join(root_dir, lib_type, "alpha_div", "plots"))
        if not os.path.exists(root_dir + "/" + lib_type + "/beta_div/plots"):
            os.makedirs(os.path.join(root_dir, lib_type, "beta_div", "plots"))
        if not os.path.exists(root_dir + "/" + lib_type + "/alpha_div/saturation"):
            os.makedirs(os.path.join(root_dir, lib_type, "alpha_div", "saturation"))
        if not os.path.exists(root_dir + "/" + lib_type + "/alpha_div/saturation/plots"):
            os.makedirs(os.path.join(root_dir, lib_type, "alpha_div", "saturation", "plots"))
        #if not os.path.exists(root_dir + "/" + lib_type + "/beta_div/upgma"):
        #    os.makedirs(os.path.join(root_dir, lib_type, "beta_div/upgma"))
        #if not os.path.exists(root_dir + "/" + lib_type + "/picrust"):
        #    os.makedirs(os.path.join(root_dir, lib_type, "picrust"))

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_nanuq_readset_file(self.args.readsets.name)
        return self._readsets


    def __init__(self):
        version = open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "..", "VERSION"), 'r').read().split('\n')[0]
        amplicontagger_string = """
###############################################################################
     _                    _ _               _____                           
    / \   _ __ ___  _ __ | (_) ___ ___  _ _|_   _|_ _  __ _  __ _  ___ _ __ 
   / _ \ | '_ ` _ \| '_ \| | |/ __/ _ \| '_ \| |/ _` |/ _` |/ _` |/ _ \ '__|
  / ___ \| | | | | | |_) | | | (_| (_) | | | | | (_| | (_| | (_| |  __/ |   
 /_/   \_\_| |_| |_| .__/|_|_|\___\___/|_| |_|_|\__,_|\__, |\__, |\___|_|   
                   |_|                                |___/ |___/         

               Support: jtremblay514@gmail.com
             Home page: jtremblay.github.io/amplicontagger.html
               Version: """ + version + """

###############################################################################"""
        sys.stderr.write(amplicontagger_string + '\n')
        time.sleep(1)
        # Add pipeline specific arguments
        self.argparser.add_argument("-b", "--barcodes", help="barcodes in fasta format", type=argparse.FileType('r'))
        self.argparser.add_argument("-e", "--external-infile", help="external fastq", type=argparse.FileType('r'), required=False)
        self.set_local_variables()
        super(AmpliconTagger, self).__init__()
                
AmpliconTagger().submit_jobs()
