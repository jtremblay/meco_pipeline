#!/usr/bin/env python

#LICENSE AND COPYRIGHT

#Copyright (C) 2023 INRS - Centre Armand-Frappier

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

# Append caf_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# CAF Modules
from core.config import *
from core.job import *
from core.pipeline import *
#from bio.design import *
from bio.readset import *

from bio import shotgun_metagenomics
from bio import microbial_ecology

from pipelines import common

# Global scope variables
log = logging.getLogger(__name__)

class MAGrefiner(common.CAFPipeline):

    """
    Start by processing raw reads for QC (trim + duk) and align them on the metagenome
    """

    def map_reads_on_bins(self):
    
        jobs = []
    
        min_contig_length = str(config.param("spades_mag", "min_contig_length", 1, "posint"))
        number_of_rounds = config.param('DEFAULT', 'number_of_rounds', 1, 'posint')
        for j in range(1,number_of_rounds):
            for bin in self.magsets:   
                log.info("MAG name: " + bin.name)
                log.info("MAG file: " + bin.fasta_file)
                if not os.path.exists(os.path.join("refined_bins", bin.name)):
                    os.makedirs(os.path.join("refined_bins", bin.name))
                
                # Will make index for bwa. and also bed file for computing reads spanning later.
                if(j == 1):
                    reference_contigs = os.path.join(bin.fasta_file)
                    #bwt_contigs = os.path.join(bin.fasta_file + ".bwt")
                else:
                    reference_contigs = os.path.join("refined_bins", bin.name, "gapfilling_round_" + str(j-1), "outgf.gapfilled.final.fa")

                #job = shotgun_metagenomics.make_index(
                #    reference_contigs,
                #    bwt_contigs
                #)
                #job.name = "make_index_contigs" + "_" + bin.name + "_round-" + str(j) 
                #job.subname = "make_index"
                #jobs.append(job)

                bam_list = []
                fastqs_R1 = []
                fastqs_R2 = []

                for readset in self.readsets:
                    #if(j == 1):
                    fastq_R1 = os.path.join("qced_reads", readset.name + "_R1.fastq.gz")
                    fastq_R2 = os.path.join("qced_reads", readset.name + "_R2.fastq.gz")
                    #else:
                    #    sys.stderr.write("Not implemented yet")

                    if not os.path.exists(os.path.join("refined_bins", bin.name, "bams_round_" + str(j))):
                        os.makedirs(os.path.join("refined_bins", bin.name, "bams_round_" + str(j)))
                    if not os.path.exists(os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j))):
                        os.makedirs(os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j)))
                    if not os.path.exists(os.path.join("refined_bins", bin.name, "assembly_round_" + str(j))):
                        os.makedirs(os.path.join("refined_bins", bin.name, "assembly_round_" + str(j)))
                    if not os.path.exists(os.path.join("refined_bins", bin.name, "scaffolding_round_" + str(j))):
                        os.makedirs(os.path.join("refined_bins", bin.name, "scaffolding_round_" + str(j)))
                    if not os.path.exists(os.path.join("refined_bins", bin.name, "gapfilling_round_" + str(j))):
                        os.makedirs(os.path.join("refined_bins", bin.name, "gapfilling_round_" + str(j)))
                    
                    bam_file = os.path.join("refined_bins", bin.name, "bams_round_" + str(j), readset.name + ".bam") 
                    bam_list.append(bam_file)
                    fastqs_R1.append(fastq_R1)
                    fastqs_R2.append(fastq_R2)
                    
                ## map reads against contigs of current bin/mag
                job = shotgun_metagenomics.bbmap_mag(
                    fastqs_R1,
                    fastqs_R2,
                    os.path.join("refined_bins", bin.name, "bams_round_" + str(j), bin.name + "_mapped.fq"),
                    reference_contigs,
                    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R1.fastq.gz"),
                    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R2.fastq.gz")
                )
                job.name = "bbmap-contigs" + "_" + bin.name + "_round-" + str(j)
                job.subname = "bbmap_mag"
                jobs.append(job)

                #job = shotgun_metagenomics.convert_bams_to_fastqs(
                #    os.path.join("refined_bins", bin.name, "bams_round_" + str(j), bin.name + "_mapped.bam"),
                #    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R1.fastq.gz"),
                #    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R2.fastq.gz")
                #)
                #job.name = "bams_to_fastqs" + "_" + bin.name + "_round-" + str(j)
                #job.subname = "bams_to_fastq"
                #jobs.append(job)
                
                job = shotgun_metagenomics.spades_mag(
                    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R1.fastq.gz"),
                    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R2.fastq.gz"),
                    reference_contigs,
                    os.path.join("refined_bins", bin.name, "assembly_round_" + str(j))
                )
                job.name = "spades_mag" + "_" + bin.name + "_round-" + str(j)
                job.subname = "spades_mag"
                jobs.append(job)

                #OPERA-LG + GapFiller
                
                # first map fastqs on new assembly
                #job = shotgun_metagenomics.bbmap_paired(
                #    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R1.fastq.gz"),
                #    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R2.fastq.gz"),
                #    os.path.join("refined_bins", bin.name, "assembly_round_" + str(j), "scaffolds_gt" + min_contig_length + ".fasta"),
                #    os.path.join("refined_bins", bin.name, "bams_round_" + str(j), bin.name + "_reassembled.bam")
                #)
                #job.name = "bbmap_paired-contigs_reassembled" + "_" + bin.name + "_round-" + str(j)
                #job.subname = "bbmap_paired"
                #jobs.append(job)

                job = shotgun_metagenomics.sspace(
                    os.path.join("refined_bins", bin.name, "assembly_round_" + str(j), "scaffolds_gt" + min_contig_length + ".fasta"),
                    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R1.fastq.gz"),
                    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R2.fastq.gz"),
                    os.path.join("refined_bins", bin.name, "scaffolding_round_" + str(j))
                )
                job.name = "sspace" + "_" + bin.name + "_round-" + str(j)
                job.subname = "sspace"
                jobs.append(job)

#            job = shotgun_metagenomics.opera_lg(
#                os.path.join("refined_bins", bin.name, "assembly_round_" + str(j), "scaffolds.fasta"),
#                os.path.join("refined_bins", bin.name, "bams_round_" + str(j), bin.name + "_reassembled.bam"),
#                os.path.join("refined_bins", bin.name, "scaffolding_round_" + str(j))
#            )
#            job.name = "opera-lg" + "_" + bin.name + "_round-" + str(j)
#            job.subname = "opera-lg"
#            jobs.append(job)
                
                job = shotgun_metagenomics.baseclear_gapfiller(
                    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R1.fastq.gz"),
                    os.path.join("refined_bins", bin.name, "fastqs_round_" + str(j), bin.name + "_mapped_R2.fastq.gz"),
                    os.path.join("refined_bins", bin.name, "scaffolding_round_" + str(j), "standard_output.final.scaffolds.fasta"),
                    os.path.join("refined_bins", bin.name, "gapfilling_round_" + str(j)),
                    bin.name + "_gapfilling_round_" + str(j)
                )
                job.name = "gapfiller" + "_" + bin.name + "_round-" + str(j)
                job.subname = "gapfiller"
                jobs.append(job)

        #TODO: Implement quast
        
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
    

    @property
    def steps(self):
        
        return [
            # Core steps.
            self.map_reads_on_bins#1
        ]

    def set_local_variables(self):
        self._parser_local = self.argparser

        
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
        def mkdir_p(path):
            try:
                os.makedirs(path)
            except OSError as exc: # Python >2.5
                if exc.errno == errno.EEXIST and os.path.isdir(path):
                    pass
                else: raise
         
        mkdir_p(root_dir)
        mkdir_p(os.path.join(root_dir, "refined_bins"))

    
    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_readset_file(self.args.readsets.name)
        return self._readsets
    
    @property
    def magsets(self):
        if not hasattr(self, "_magsets"):
            self._magsets = parse_magset_file(self.args.magsets.name)
        return self._magsets


    def __init__(self):
        metagenomics_string = """
###############################################################################
	      __  __          _____             __ _                 
	     |  \/  |   /\   / ____|           / _(_)                
	     | \  / |  /  \ | |  __   _ __ ___| |_ _ _ __   ___ _ __ 
	     | |\/| | / /\ \| | |_ | | '__/ _ \  _| | '_ \ / _ \ '__|
	     | |  | |/ ____ \ |__| | | | |  __/ | | | | | |  __/ |   
	     |_|  |_/_/    \_\_____| |_|  \___|_| |_|_| |_|\___|_|   
								 
               Support: jtremblay514@gmail.com
             Home page: jtremblay.github.io/pipelines.html

###############################################################################"""
        sys.stderr.write(metagenomics_string + '\n')
        time.sleep(1)
        # Add pipeline specific arguments
        self.argparser.add_argument("-m", "--magsets", help="magset file", type=file, required=True)
        self.argparser.add_argument("-r", "--readsets", help="readset file", type=file, required=True)
        self.set_local_variables()
        sys.stderr.write('MAG refiner pipeline\n')
        super(MAGrefiner, self).__init__()
                
MAGrefiner().submit_jobs()
