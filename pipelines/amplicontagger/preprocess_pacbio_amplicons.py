#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys
import errno

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.readset import *
#from bio.utils import *

#from bio import shotgun_metagenomics
#from bio import shotgun_metatranscriptomics
#from bio import microbial_ecology
from bio import pacbio_tools
#from bio import utils

from pipelines.illumina import illumina

# Global scope variables
log = logging.getLogger(__name__)

class Preprocesspacbio(illumina.Illumina):

    def preprocess(self):
        jobs = []
        outdir = self._root_dir
        
        for readset in self.readsets:
            prefix = os.path.join(outdir, "pacbio_preprocessed_reads", readset.name + "_phred33")
            prefix2 = os.path.join(outdir, "pacbio_preprocessed_reads", readset.name + "_phred33_nic")
            #prefix2 = os.path.join(outdir, "pacbio_preprocessed_reads", readset.name + "_phred33_1400")
            #prefix3 = os.path.join(outdir, "pacbio_preprocessed_reads", readset.name + "_phred33_1400_nic")

            if readset.run_type == "PAIRED_END":
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (PAIRED_END not implemented yet)!")

            elif readset.run_type == "SINGLE_END":
                #sys.stderr.write('outdir: ' + readset.fastq1 + '\n')
                if not os.path.exists(os.path.join(outdir, "pacbio_preprocessed_reads")):
                    os.makedirs(os.path.join(outdir, "pacbio_preprocessed_reads"))
                
                job = pacbio_tools.gunzip(
                    readset.fastq1,
                    os.path.join(outdir, "pacbio_preprocessed_reads", readset.name + ".fastq")
                )
                job.name = readset.name
                job.subname = "gunzip_" + readset.name
                jobs.append(job)

                job = pacbio_tools.fix_pacbio_qscores(
                    os.path.join(outdir, "pacbio_preprocessed_reads", readset.name + ".fastq"),
                    prefix + ".fastq"
                )
                job.name = "fix_qscores_" + readset.name
                job.subname = "trim"
                jobs.append(job)
                
                #job = pacbio_tools.cut_fastq_to_specific_length(
                #    prefix + ".fastq",
                #    prefix2 + ".fastq"
                #)
                #job.name = "cut_" + readset.name
                #job.subname = "trim"
                #jobs.append(job)
                
                job = pacbio_tools.replace_illegal_characters(
                    prefix + ".fastq",
                    prefix2 + ".fastq"
                )
                job.name = "replace_illegal_char_" + readset.name
                job.subname = "trim"
                jobs.append(job)
                
                job = pacbio_tools.gzip_to_raw_reads_directory(
                    os.path.join(outdir, "pacbio_preprocessed_reads", readset.name + "_phred33_nic.fastq"),
                    os.path.join(outdir, "raw_reads", readset.name + "_R1.fastq.gz")
                )
                job.name = "gzip_to_raw_reads_directory_" + readset.name
                job.subname = "gzip"
                jobs.append(job)
            
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
            
        
        
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
            self.preprocess#1
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
        if not os.path.exists(root_dir + "/raw_reads"):
            os.makedirs(os.path.join(root_dir, "raw_reads"))
#        if not os.path.exists(root_dir + "/qced_reads"):
#            os.makedirs(os.path.join(root_dir, "qced_reads"))
#        if not os.path.exists(root_dir + "/abundance"):
#            os.makedirs(os.path.join(root_dir, "abundance"))
#        if not os.path.exists(root_dir + "/taxonomy"):
#            os.makedirs(os.path.join(root_dir, "taxonomy"))
#        if not os.path.exists(root_dir + "gene_annotation/taxonomy_genes/all/absolute/plots"):
#            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/all/absolute/plots"))
#            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/all/relative/plots"))
#        if not os.path.exists(root_dir + "/gene_annotation/taxonomy_genes/others/absolute/plots"):
#            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/others/absolute/plots"))
#            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/others/relative/plots"))
#        if not os.path.exists(root_dir + "gene_annotation/taxonomy_genes/bacteriaArchaea/absolute/plots"):
#            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/bacteriaArchaea/absolute/plots"))
#            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/bacteriaArchaea/relative/plots"))
#        if not os.path.exists(root_dir + "/statistics"):
#            os.makedirs(os.path.join(root_dir, "statistics"))


    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            #self._readsets = parse_nanuq_readset_file(self.args.readsets.name)
            self._readsets = parse_readset_file(self.args.readsets.name)
        return self._readsets


    def __init__(self):
        # Add pipeline specific arguments
        version = open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "..", "VERSION"), 'r').read().split('\n')[0]
        welcome_string = """
################################################################################
Preprocessing of PacBio amplicon data type for the AmpliconTagger pipeline.
Version: """ + version + """
################################################################################
"""
        sys.stderr.write(welcome_string + '\n')
        time.sleep(1)
        self.set_local_variables()
        sys.stderr.write('Pacbio preprocessing pipeline\n')
        super(Preprocesspacbio, self).__init__()
                
Preprocesspacbio().submit_jobs()
