#!/usr/bin/env python

#LICENSE AND COPYRIGHT

#Copyright (C) 2019 INRS - Centre Armand-Frappier

#This license does not grant you the right to use any trademark, service
#mark, tradename, or logo of the Copyright Holder.

#This license includes the non-exclusive, worldwide, free-of-charge
#patent license to make, have made, use, offer to sell, sell, import and
#otherwise transfer the Package with respect to any patent claims
#licensable by the Copyright Holder that are necessarily infringed by the
#Package. If you institute patent litigation (including a cross-claim or
#counterclaim) against any party alleging that the Package constitutes
#direct or contributory patent infringement, then this Artistic License
#to you shall terminate on the date that such litigation is filed.

#Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
#AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
#THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
#PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
#YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
#CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
#CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
#EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#Author: Julien Tremblay - jtremblay514@gmail.com

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys
import errno
import time

# Append pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# GenPipes/CAF Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.readset import *

from bio import shotgun_metagenomics
from pipelines import common

# Global scope variables
log = logging.getLogger(__name__)

class Metagenomics(common.NRCPipeline):
    """
    Shotgun Metagenomics Pipeline
    Written by Julien Tremblay
    ========================

    Pipeline that performs reads QC, de novo co-assembly (Megahit), read mapping against co-assembly,
    contigs coverage profiles, contigs abundance matrix by sample. From co-assembly are also done: gene prediction,
    gene coverage profiles, gene abundance matrix by sample. MAG or metagenome binning can also be performed.

    """

    def trim(self):
        
        """
        Step trim(): Raw fastqs will be trimmed using Trimmomatic. Interleaved fastqs will be generated after trimming. 
        """
        
        jobs = []
        
        trim_stats = []
        
        for readset in self.readsets:
            trim_file_prefix = os.path.join("qced_reads", readset.sample.name, readset.name + ".trim.")
            trim_file_prefix_long = os.path.join("qced_reads", readset.sample.name, readset.name + ".trim.long.")
            if not os.path.exists(os.path.join("qced_reads", readset.sample.name)):
                os.makedirs(os.path.join("qced_reads", readset.sample.name))
            
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
                job.name = "trimmomatic_" + readset.sample.name
                job.subname = "trim"
                jobs.append(job) 

                trim_stats.append(trim_file_prefix + "stats.csv")
            
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

            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END or SINGLE_END_LONG)!")
        
        return jobs
            
    def remove_contam(self):
        """
        Step remove_contam(): Trimmed fastqs will be filtered for contaminant sequences (e.g. Illumina adapters,
                              known primer sequences, etc). A second round of contaminant filtering will be done 
                              to filter out PhiX sequences which are usually spiked-in in Illumina sequencing runs.
        """
        jobs=[]
        
        logs = []
        readset_ids = []

        for readset in self.readsets:
            trim_file_prefix = os.path.join("qced_reads", readset.sample.name, readset.name + ".trim.")
            trim_file_prefix_long = os.path.join( "qced_reads", readset.sample.name, readset.name + ".trim.long.")
            outfile_prefix = os.path.join("qced_reads", readset.sample.name, readset.name + ".")
            outfile_prefix_long = os.path.join("qced_reads", readset.sample.name, readset.name + ".long.")
            out1up = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_unpaired_R1.fastq.gz")
            out2up = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_unpaired_R2.fastq.gz")
            out1p = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
            out2p = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
            
            if readset.run_type == "PAIRED_END":
                log = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".duk_contam_interleaved_log.txt")
                readset_id = readset.name

                job = shotgun_metagenomics.bbduk(
                    trim_file_prefix + "interleaved.fastq.gz",
                    outfile_prefix + "contam.fastq",
                    outfile_prefix + "ncontam.fastq.gz",
                    log,
                    config.param('DB', 'contaminants', 1, 'filepath'),
                    trim_file_prefix + "interleaved.fastq.gz"
                )
                job.name = "bbduk_interleaved_" + readset.sample.name
                job.subname = "duk"
                jobs.append(job)
            
        return jobs
     
    def assembly(self):
        
        """
        Step assembly(): Perform assembly using Megahit (short reads) or miniasm (long reads).
        """

        jobs = []
        
        fastq_list = []
        fastq_gz_list = []
     
        for readset in self.readsets:
            
            #Here we assume all libs inside readset file are either PAIRED or SINGLE, but not both.
            if readset.run_type == "PAIRED_END":
                reads_type = "pe"
                fastq_gz_list.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz"))
                fastq_list.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq"))
       
        if config.param("DEFAULT", "assembler", 1, "string") == "megahit":
            job = shotgun_metagenomics.megahit(
                fastq_gz_list,
                os.path.join(root, "assembly"),
                reads_type
            )
            job.name = "megahit"
            job.subname = "megahit"
            jobs.append(job)

        
        
        else:
            raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END and assembler=megahit OR SINGLE_END and assembler=miniasm)!")

        job = shotgun_metagenomics.compile_assembly_stats(
            os.path.join(root, "assembly", "Contigs.fasta"),
            os.path.join(root, "assembly", "assembly_stats.txt")
        )
        job.name = "compile_assembly_stats"
        job.subname = "compile_assembly_stats"
        jobs.append(job)

        return jobs

    
    def cleanup(self):
        
        """
        Step cleanup() : Delete unnecessary intermediate files. Compress end results/files.
        """

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
        
        return [
            self.trim,
            self.remove_contam,
            self.assembly
        ]

    def set_local_variables(self):
        self._parser_local = self.argparser

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
        def mkdir_p(path):
            try:
                os.makedirs(path)
            except OSError as exc: # Python >2.5
                if exc.errno == errno.EEXIST and os.path.isdir(path):
                    pass
                else: raise
         
        mkdir_p(root_dir)
        mkdir_p(os.path.join(root_dir, "qced_reads"))
        mkdir_p(os.path.join(root_dir, "assembly"))

        #prefixes = ["5S", "8S", "16S", "18S", "23S", "28S"]
        prefixes = ["16S"]
        for prefix in prefixes:
            mkdir_p(os.path.join(self._root_dir, "annotations", "blastn_nt_contigs"))
            mkdir_p(os.path.join(self._root_dir, "annotations", "taxonomy_consensus"))
       
        prefixes = ["consensus"]
        for prefix in prefixes:
            mkdir_p(os.path.join(root_dir, "annotations/taxonomy_" + prefix + "/all/absolute/plots"))
            mkdir_p(os.path.join(root_dir, "annotations/taxonomy_" + prefix + "/all/relative/plots"))
            mkdir_p(os.path.join(root_dir, "annotations/taxonomy_" + prefix + "/others/absolute/plots"))
            mkdir_p(os.path.join(root_dir, "annotations/taxonomy_" + prefix + "/others/relative/plots"))
            mkdir_p(os.path.join(root_dir, "annotations/taxonomy_" + prefix + "/bacteriaArchaea/absolute/plots"))
            mkdir_p(os.path.join(root_dir, "annotations/taxonomy_" + prefix + "/bacteriaArchaea/relative/plots"))
    
    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_readset_file(self.args.readsets.name)
        return self._readsets


    def __init__(self):
        metagenomics_string = """
###############################################################################
              _____ _           _                    __  __  _____ 
             / ____| |         | |                  |  \/  |/ ____|
            | (___ | |__   ___ | |_ __ _ _   _ _ __ | \  / | |  __ 
             \___ \| '_ \ / _ \| __/ _` | | | | '_ \| |\/| | | |_ |
             ____) | | | | (_) | || (_| | |_| | | | | |  | | |__| |
            |_____/|_| |_|\___/ \__\__, |\__,_|_| |_|_|  |_|\_____|
                                    __/ |                          
                                   |___/                           
   
               Support: jtremblay514@gmail.com
             Home page: jtremblay.github.io/pipelines.html

###############################################################################"""
        sys.stderr.write(metagenomics_string + '\n')
        time.sleep(1)
        # Add pipeline specific arguments
        self.argparser.add_argument("-n", "--normalize-before-assembly", help="Nomalize before assembly (default: false)", action="store_true", required=False)
        self.argparser.add_argument("-r", "--readsets", help="readset file", type=file, required=False)
        self.set_local_variables()
        sys.stderr.write('Metagenomics pipeline\n')
        super(Metagenomics, self).__init__()
                
Metagenomics().submit_jobs()
