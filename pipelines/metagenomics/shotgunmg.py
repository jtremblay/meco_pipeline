#!/usr/bin/env python

#LICENSE AND COPYRIGHT

#Copyright (C) 2023 Julien Tremblay

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

# GenPipes/MECO Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.readset import *

from bio import shotgun_metagenomics
from bio import microbial_ecology

from pipelines import common
#from pipelines.illumina import illumina
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# Global scope variables
log = logging.getLogger(__name__)

class Metagenomics(common.MECOPipeline):
    """
    Shotgun Metagenomics Pipeline
    Written by Julien Tremblay, Genomics and Microbiomes, Julien Tremblay
    ========================

    Pipeline that performs reads QC, de novo co-assembly (Megahit), read mapping against co-assembly,
    contigs coverage profiles, contigs abundance matrix by sample. From co-assembly are also done: gene prediction,
    gene coverage profiles, gene abundance matrix by sample. MAG or metagenome binning can also be performed.

    """

    def trim(self):
        
        """
        Step trim(): Raw fastqs will be trimmed using Trimmomatic. Interleaved fastqs will be generated after trimming. 
        """
        
        mkdir_p(os.path.join(self._root_dir, "qced_reads"))
        
        jobs = []
        trim_stats = []
        
        for readset in self.readsets:
            trim_file_prefix = os.path.join("qced_reads", readset.sample.name, readset.name + ".trim.")
            trim_file_prefix_long = os.path.join("qced_reads", readset.sample.name, readset.name + ".trim.long.")
            if not os.path.exists(os.path.join("qced_reads", readset.sample.name)):
                os.makedirs(os.path.join("qced_reads", readset.sample.name))
            
            if readset.run_type == "PAIRED_END":

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
            
            elif readset.run_type == "SINGLE_END" :
                if config.param("DEFAULT", "qc_methods", 1, "string") == "trimmomatic":

                    job = shotgun_metagenomics.trimmomatic_se(
                        readset.fastq1,
                        trim_file_prefix + "pair1.fastq.gz",
                        readset.quality_offset,
                        trim_file_prefix + "out",
                        trim_file_prefix + "stats.csv"
                    )
                    job.name = "trimmomatic_" + readset.sample.name
                    job.subname = "trim"
                    jobs.append(job) 
    
                    trim_stats.append(trim_file_prefix + "stats.csv")
                
                # For nanopore using trimmomatic
                elif config.param("DEFAULT", "qc_methods", 1, "string") == "porechop,trimmomatic":
                    
                    job = shotgun_metagenomics.porechop(
                        readset.fastq1,
                        trim_file_prefix + "porechop.fastq.gz"
                    )
                    job.name = "porechop_" + readset.sample.name
                    job.subname = "porechop"
                    jobs.append(job) 
                    
                    job = shotgun_metagenomics.trimmomatic_se(
                        trim_file_prefix + "porechop.fastq.gz",
                        trim_file_prefix + "pair1.fastq.gz",
                        readset.quality_offset,
                        trim_file_prefix + "out",
                        trim_file_prefix + "stats.csv"
                    )
                    job.name = "trimmomatic_" + readset.sample.name
                    job.subname = "trim"
                    jobs.append(job) 
                
            # SPE_LSE = Short Paired End and Long Single End
            elif readset.run_type == "SPE_LSE":
                #First do SPE 
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
            
                job = shotgun_metagenomics.create_interleaved_fastq(
                    trim_file_prefix + "pair1.fastq.gz",
                    trim_file_prefix + "pair2.fastq.gz",
                    trim_file_prefix + "interleaved.fastq",
                    trim_file_prefix + "interleaved.fastq.gz"
                )
                job.name = "create_interleaved_fastq_" + readset.sample.name
                job.subname = "interleaved_fastq"
                jobs.append(job) 
                
                if config.param("DEFAULT", "qc_methods", 1, "string") == "trimmomatic":
                    job = shotgun_metagenomics.trimmomatic_se(
                        readset.fastq_long,
                        trim_file_prefix_long + ".pair1.fastq.gz",
                        readset.quality_offset,
                        trim_file_prefix_long + "out",
                        trim_file_prefix_long + "stats.csv"
                    )
                    job.name = "trimmomatic_" + readset.sample.name
                    job.subname = "trim"
                    jobs.append(job) 

                    trim_stats.append(trim_file_prefix_long + "stats.csv")
                
                # For nanopore using trimmomatic
                elif config.param("DEFAULT", "qc_methods", 1, "string") == "porechop,trimmomatic":
                    
                    job = shotgun_metagenomics.porechop(
                        readset.fastq1,
                        trim_file_prefix_long + "porechop.fastq.gz"
                    )
                    job.name = "porechop_" + readset.sample.name
                    job.subname = "porechop"
                    jobs.append(job) 
                    
                    job = shotgun_metagenomics.trimmomatic_se(
                        trim_file_prefix_long + "porechop.fastq.gz",
                        trim_file_prefix_long + "pair1.fastq.gz",
                        readset.quality_offset,
                        trim_file_prefix_long + "out",
                        trim_file_prefix_long + "stats.csv"
                    )
                    job.name = "trimmomatic_" + readset.sample.name
                    job.subname = "trim"
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
        jobs = []
        logs = []
        logs_sub = []
        readset_ids = []

        for readset in self.readsets:
            trim_file_prefix = os.path.join("qced_reads", readset.sample.name, readset.name + ".trim.")
            trim_file_prefix_long = os.path.join("qced_reads", readset.sample.name, readset.name + ".trim.long.")
            outfile_prefix = os.path.join("qced_reads", readset.sample.name, readset.name + ".")
            outfile_prefix_long = os.path.join("qced_reads", readset.sample.name, readset.name + ".long.")
            out1up = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_unpaired_R1.fastq.gz")
            out2up = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_unpaired_R2.fastq.gz")
            out1p = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
            out2p = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
            
            if readset.run_type == "PAIRED_END":
                log = os.path.join("qced_reads", readset.sample.name, readset.name + ".duk_contam_interleaved_log.txt")

                job = shotgun_metagenomics.bbduk_paired(
                    os.path.join("qced_reads", readset.sample.name, readset.name + ".trim." + "pair1.fastq.gz"),
                    os.path.join("qced_reads", readset.sample.name, readset.name + ".trim." + "pair2.fastq.gz"),
                    os.path.join("qced_reads", readset.sample.name, readset.name + "." + "contam_R1.fastq"),
                    os.path.join("qced_reads", readset.sample.name, readset.name + "." + "contam_R2.fastq"),
                    os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz"),
                    os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz"),
                    log,
                    config.param('DB', 'contaminants', 1, 'filepath')
                )
                job.name = "bbduk_paired_" + readset.sample.name
                job.subname = "duk"
                jobs.append(job)
            
                logs.append(log)
                readset_ids.append(readset.name)
                
                ref_genome = config.param('DB', 'ref_genome', 0, 'string')
                #sys.stderr.write("type(ref_genome): " + str(type(ref_genome)) + "\n")
                # if reference genome provided for reads substraction
                if isinstance(ref_genome, str) and ref_genome != "":
                    job = shotgun_metagenomics.bbmap_subtract(
                        os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz"),
                        os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz"),
                        os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_R1.fastq.gz"),
                        os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_R2.fastq.gz"),
                        ref_genome,
                        os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_log.txt")
                    )
                    job.name = "bbmap_subtract_" + readset.sample.name
                    job.subname = "bbmap_sub"
                    jobs.append(job)

                    logs_sub.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_log.txt"))
            
            elif readset.run_type == "SINGLE_END":
                if isinstance(ref_genome, str) and ref_genome != "":
                    raise Exception("Error: run type SINGLE_END not yet implemented with ref_genome subtracting option..\n")
                
                log = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".duk_contam_pair1_log.txt")

                job = shotgun_metagenomics.bbduk(
                    trim_file_prefix + "pair1.fastq.gz",
                    outfile_prefix + "contam.fastq",
                    outfile_prefix + "ncontam.fastq.gz",
                    log,
                    config.param('DB', 'contaminants', 1, 'filepath'),
                    trim_file_prefix + "pair1.fastq.gz"
                )
                job.name = "bbduk_single_end_reads_" + readset.sample.name
                job.subname = "duk"
                jobs.append(job)
                
                logs.append(log)
                readset_ids.append(readset.name)
            
            elif readset.run_type == "SPE_LSE":
                if isinstance(ref_genome, str) and ref_genome != "":
                    raise Exception("Error: run type SPE_LSE not yet implemented with ref_genome subtracting option..\n")

                log = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".SPE.duk_contam_interleaved_log.txt")
                readset_id = readset.name

                job = shotgun_metagenomics.bbduk(
                    trim_file_prefix + "interleaved.fastq.gz",
                    outfile_prefix + "contam.fastq",
                    outfile_prefix + "ncontam.fastq.gz",
                    log,
                    config.param('DB', 'contaminants', 1, 'filepath'),
                    trim_file_prefix + "interleaved.fastq.gz"
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
                
                logs.append(log)
                readset_ids.append(readset.name)
    
                ## Then do long SE reads.
                log = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".long_duk_contam_pair1_log.txt")

                job = shotgun_metagenomics.bbduk(
                    trim_file_prefix_long + "pair1.fastq.gz",
                    outfile_prefix_long + "contam.fastq",
                    outfile_prefix_long + "ncontam.fastq.gz",
                    log,
                    config.param('DB', 'contaminants', 1, 'filepath'),
                    trim_file_prefix_long + "pair1.fastq.gz"
                )
                job.name = "duk_single_end_long_reads_" + readset.sample.name
                job.subname = "duk"
                jobs.append(job)
                
                logs.append(log)
                
                sys.stderr.write("[DEBUG] " + outfile_prefix_long + "ncontam.fastq.gz\n")
        
        # Compile duk logs.
        if isinstance(ref_genome, str) and ref_genome != "":
            job = shotgun_metagenomics.merge_duk_logs_interleaved(
                logs,
                readset_ids,
                os.path.join("qced_reads", "duk_merged_logs.tsv")
            )
            job.name = "merge_duk_logs"
            job.subname = "merge_duk_logs"
            jobs.append(job)
           
            job = shotgun_metagenomics.merge_duk_sub_logs_interleaved(
                logs_sub,
                readset_ids,
                os.path.join("qced_reads", "duk_merged_sub_logs.tsv")
            )
            job.name = "merge_duk_logs_sub"
            job.subname = "merge_duk_logs"
            jobs.append(job)

            
        else:    
            job = shotgun_metagenomics.merge_duk_logs_interleaved(
                logs,
                readset_ids,
                os.path.join("qced_reads", "duk_merged_logs.tsv")
            )
            job.name = "merge_duk_logs"
            job.subname = "merge_duk_logs"
            jobs.append(job)
            
        return jobs

    def subsample(self):
        
        """
        Step subsample(): Subsample prior to assembly. As discussed in https://doi.org/10.1093/bib/bbac443
        """

        jobs = []
        sample_list = []
        fastq_gz_list_R1 = []
        fastq_gz_list_R2 = []
        fastq_gz_out_list_R1 = []
        fastq_gz_out_list_R2 = []
        
        subsample = config.param('DEFAULT', 'subsample_reads_prior_to_assembly', required=False, type='string') 
        ref_genome = config.param('DB', 'ref_genome', required=False, type='string') 
     
        if isinstance(subsample, str) and subsample == "yes":
            for readset in self.readsets:
                #Here we assume all libs inside readset file are either PAIRED or SINGLE, but not both.
                if readset.run_type == "PAIRED_END":
                    if isinstance(ref_genome, str) and ref_genome != "":
                        sample_list.append(readset.name)
                        fastq_gz_list_R1.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_R1.fastq.gz"))
                        fastq_gz_list_R2.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_R2.fastq.gz"))
                        fastq_gz_out_list_R1.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_subsampled_R1.fastq.gz"))
                        fastq_gz_out_list_R2.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_subsampled_R2.fastq.gz"))
                    else:
                        sample_list.append(readset.name)
                        fastq_gz_list_R1.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz"))
                        fastq_gz_list_R2.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz"))
                        fastq_gz_out_list_R1.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_subsampled_R1.fastq.gz"))
                        fastq_gz_out_list_R2.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_subsampled_R2.fastq.gz"))
          
            for i in range(0, len(fastq_gz_list_R1)):
                job = shotgun_metagenomics.subsample(
                    fastq_gz_list_R1[i],
                    fastq_gz_list_R2[i],
                    fastq_gz_out_list_R1[i],
                    fastq_gz_out_list_R2[i]
                )
                job.name = "subsample_" + sample_list[i]
                job.subname = "subsample"
                jobs.append(job)
        
        return jobs

    def assembly(self):
        
        """
        Step assembly(): Perform assembly using Megahit or meta SPAdes (short reads). miniasm (long reads) in development.
        """
        mkdir_p(os.path.join(self._root_dir, "assembly"))
        ref_genome = config.param('DB', 'ref_genome', required=False, type='string') 
        subsample = config.param('DEFAULT', 'subsample_reads_prior_to_assembly', required=False, type='string') 

        root = self._root_dir
        jobs = []
        fastq_gz_list_R1 = []
        fastq_gz_list_R2 = []
        fastq_gz_list_merged = []
        short_fastq_gz_list = []
        long_fastq_gz_list = []
        long_fasta_list = []
        long_corrected_fasta_list = []
        paf_list = []
     
        for readset in self.readsets:
            #Here we assume all libs inside readset file are either PAIRED or SINGLE, but not both.
            if readset.run_type == "PAIRED_END":
                reads_type = "pe"
                if config.param("DEFAULT", "assembler", 1, "string") == "megahit_merged":
                    fastq_gz_list_R1.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".notCombined_1.fastq.gz"))
                    fastq_gz_list_R2.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".notCombined_2.fastq.gz"))
                    fastq_gz_list_merged.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".extendedFrags.fastq.gz"))
                else:
                    if isinstance(ref_genome, str) and ref_genome != "":
                        if isinstance(subsample, str) and subsample == "yes":
                            fastq_gz_list_R1.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_subsampled_R1.fastq.gz"))
                            fastq_gz_list_R2.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_subsampled_R2.fastq.gz"))
                        else:
                            fastq_gz_list_R1.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_R1.fastq.gz"))
                            fastq_gz_list_R2.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_R2.fastq.gz"))
                    else:
                        if isinstance(subsample, str) and subsample == "yes":
                            fastq_gz_list_R1.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_subsampled_R1.fastq.gz"))
                            fastq_gz_list_R2.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_subsampled_R2.fastq.gz"))
                        else:
                            fastq_gz_list_R1.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz"))
                            fastq_gz_list_R2.append(os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz"))

            elif readset.run_type == "SINGLE_END":
                fastq_gz_list.append(os.path.join(root, "qced_reads", readset.sample.name, readset.name + ".ncontam.fastq.gz"))
                reads_type = "se"
            elif readset.run_type == "SPE_LSE":
                short_fastq_gz_list.append(os.path.join(root, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz"))
                long_fastq_gz_list.append(os.path.join(root, "qced_reads", readset.sample.name, readset.name + ".long.ncontam.fastq.gz"))
                long_fasta_list.append(os.path.join(root, "qced_reads", readset.sample.name, readset.name + ".long.ncontam.fasta"))
                long_corrected_fasta_list.append(os.path.join(root, "qced_reads", readset.sample.name, readset.name + ".corrected.fasta"))
                reads_type = "spe_lse"
       
        if config.param("DEFAULT", "assembler", 1, "string") == "megahit":
            
            job = shotgun_metagenomics.megahit_R1_R2(
                fastq_gz_list_R1,
                fastq_gz_list_R2,
                os.path.join("assembly"),
                reads_type
            )
            job.name = "megahit"
            job.subname = "megahit"
            jobs.append(job)

        elif config.param("DEFAULT", "assembler", 1, "string") == "megahit_merged":
            job = shotgun_metagenomics.megahit_R1_R2_merged(
                fastq_gz_list_R1,
                fastq_gz_list_R2,
                fastq_gz_list_merged,
                os.path.join("assembly"),
                reads_type
            )
            job.name = "megahit"
            job.subname = "megahit"
            jobs.append(job)
        
        elif config.param("DEFAULT", "assembler", 1, "string") == "spades" and readset.run_type == "PAIRED_END":
            
            job = shotgun_metagenomics.spades_R1_R2(
                fastq_gz_list_R1,
                fastq_gz_list_R2,
                os.path.join("assembly")
            )
            job.name = "meta-SPAdes"
            job.subname = "spades"
            jobs.append(job)

        elif config.param("DEFAULT", "assembler", 1, "string") == "miniasm" and readset.run_type == "SINGLE_END":
            
            job = shotgun_metagenomics.cat_fastqs(
                fastq_gz_list,
                os.path.join("assembly", "all_qced_reads.fastq.gz")
            )
            job.name = "cat_fastq"
            job.subname = "cat_fastq"
            jobs.append(job)

            job = shotgun_metagenomics.minimap2_all_vs_all(
                os.path.join("assembly", "all_qced_reads.fastq.gz"),
                os.path.join("assembly", "all_qced_reads.paf.gz")
            )
            job.name = "minimap2"
            job.subname = "minimap2_all_vs_all"
            jobs.append(job)

            job = shotgun_metagenomics.miniasm(
                os.path.join("assembly", "all_qced_reads.fastq.gz"),
                os.path.join("assembly", "all_qced_reads.paf.gz"),
                os.path.join("assembly", "Contigs.gfa"),
                os.path.join("assembly", "Contigs.fasta")
            )
            job.name = "miniasm"
            job.subname = "miniasm"
            jobs.append(job)
        
        elif config.param("DEFAULT", "assembler", 1, "string") == "canu" and readset.run_type == "SINGLE_END":
            
            job = shotgun_metagenomics.canu(
                fastq_gz_list,
                os.path.join(root, "assembly")
            )
            job.name = "canu"
            job.subname = "canu"
            jobs.append(job)
            
        elif readset.run_type == "SPE_LSE" and config.param("DEFAULT", "assembler", 1, "string") == "fmlrc,miniasm":

            if not os.path.exists(os.path.join(root,"assembly", "bwt")):
                os.makedirs(os.path.join(root,"assembly", "bwt"))
            
            job = shotgun_metagenomics.fastqsgz_to_fastas(
                long_fastq_gz_list,
                long_fasta_list,
                os.path.join(root, "assembly", "bwt", "convert_fastqs_to_fastas.done")
            )
            job.name = "fastqsgz_to_fastas"
            job.subname = "fastqsgz_to_fastas"
            jobs.append(job)
            
            for readset in self.readsets:
            
                job = shotgun_metagenomics.ropebwt2(
                    os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz"),
                    os.path.join("assembly", "bwt", readset.name + ".bwt2"),
                    os.path.join("assembly", "bwt", readset.name + ".npy")
                )
                job.name = "ropebwt2"
                job.subname = "ropebwt2"
                jobs.append(job)

                job = shotgun_metagenomics.fmlrc(
                    os.path.join("assembly", "bwt", readset.name + ".npy"),
                    os.path.join("qced_reads", readset.sample.name, readset.name + ".long.ncontam.fasta"),
                    os.path.join("qced_reads", readset.sample.name, readset.name + ".corrected.fasta")
                )
                job.name = "fmlrc"
                job.subname = "fmlrc"
                jobs.append(job)
            
            # Then with corrected reads in hand, perform the minimap2 + miniasm.
            job = shotgun_metagenomics.cat_fastqs(
                long_corrected_fasta_list,
                os.path.join("assembly", "all_qced_reads.fasta.gz")
            )
            job.name = "cat_fastq"
            job.subname = "cat_fastq"
            jobs.append(job)
            
            job = shotgun_metagenomics.minimap2_all_vs_all(
                os.path.join("assembly", "all_qced_reads.fasta.gz"),
                os.path.join("assembly", "all_qced_reads.paf.gz")
            )
            job.name = "minimap2"
            job.subname = "minimap2_all_vs_all"
            jobs.append(job)

            job = shotgun_metagenomics.miniasm(
                os.path.join("assembly", "all_qced_reads.fasta.gz"),
                os.path.join("assembly", "all_qced_reads.paf.gz"),
                os.path.join("assembly", "Contigs.gfa"),
                os.path.join("assembly", "Contigs.fasta")
            )
            job.name = "miniasm"
            job.subname = "miniasm"
            jobs.append(job)
        
        else:
            raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END and assembler=megahit OR SINGLE_END and assembler=miniasm)!")

        job = shotgun_metagenomics.compile_assembly_stats(
            os.path.join("assembly", "Contigs.fasta"),
            os.path.join("assembly", "assembly_stats.txt")
        )
        job.name = "compile_assembly_stats"
        job.subname = "compile_assembly_stats"
        jobs.append(job)
        
        job = shotgun_metagenomics.quast(
            os.path.join("assembly", "Contigs.fasta"),
            os.path.join("assembly", "quast")
        )
        job.name = "quast"
        job.subname = "quast"
        jobs.append(job)
        
        job = shotgun_metagenomics.get_contigs_length_and_gc(
            os.path.join("assembly", "Contigs.fasta"),
            os.path.join("assembly", "Contigs_length_gc.tsv")
        )
        job.name = "get_contigs_length_and_gc"
        job.subname = "get_contigs_length_and_gc"
        jobs.append(job)

        return jobs

    def gene_prediction(self):
        
        """
        Step gene_prediction(): Use Prodigal to predict ORFs of bacterial genomes.
        """
        
        mkdir_p(os.path.join(self._root_dir, "gene_prediction"))
        
        jobs = []
        
        infile = os.path.join("assembly", "Contigs.fasta")
        outdir = os.path.join("gene_prediction")
        
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
        job.name = "prodigal"
        job.subname = "prodigal"
        jobs.append(job)
 
        return jobs

    def abundance(self):
        
        """
        Step abundance(): map qced reads on contigs + generate contigs and gene abundance matrices.
        """
        mkdir_p(os.path.join(self._root_dir, "gene_abundance"))
        mkdir_p(os.path.join(self._root_dir, "contig_abundance"))
        
        jobs = []
        cov_list_contigs = []
        cov_list_genes = []
        flagstats_contigs_list = []
        trimmomatic_list = []
        bbduk_list = []
        bbduk_sub_list = []
       
        ref_genome = config.param('DB', 'ref_genome', 0, 'string')

        # Will make index for bwa. and also bed file for computing reads spanning later.
        reference_contigs = os.path.join("assembly", "Contigs.fasta")
        bed_contigs = os.path.join("assembly", "Contigs.bed")
        bwt_contigs = os.path.join("assembly", "Contigs.fasta.bwt")
        bed_genes = os.path.join("gene_prediction", "Contigs_genes.bed")

        job = shotgun_metagenomics.fasta_to_bed(
            reference_contigs,
            bed_contigs
        )
        job.name = "fasta_to_bed"
        job.subname = "fasta_to_bed"
        jobs.append(job)
        
        job = shotgun_metagenomics.gff_to_bed(
            os.path.join("gene_prediction", "Contigs_renamed.gff"),
            bed_genes
        )
        job.name = "gff_to_bed"
        job.subname = "gff_to_bed"
        jobs.append(job)
                    
        if config.param('DEFAULT', 'mapper', 1, "string") == "bwa":
            job = shotgun_metagenomics.make_index(
                reference_contigs,
                bwt_contigs
            )
            job.name = "make_index_contigs_bwa"
            job.subname = "make_index"
            jobs.append(job)
            
        elif config.param('DEFAULT', 'mapper', 1, "string") == "bbmap":
            job = shotgun_metagenomics.bbmap_index(
                reference_contigs,
                os.path.join("assembly", "bbmap_ref")
                #with hidden bbmap_index.done generated as outfile
            )
            job.name = "make_index_contigs_bbmap"
            job.subname = "make_index"
            jobs.append(job)
        

        for readset in self.readsets:
            bam_contigs = os.path.join("contig_abundance", readset.sample.name, readset.name + ".bam")
            flagstats_contigs = os.path.join("contig_abundance", readset.sample.name, readset.name + ".flagstats")
            flagstats_contigs_list.append(flagstats_contigs)
            trimmomatic = os.path.join("qced_reads", readset.sample.name, readset.name + ".trim.stats.csv")
            trimmomatic_list.append(trimmomatic)
            bbduk = os.path.join("qced_reads", readset.sample.name, readset.name + ".duk_contam_interleaved_log.txt")
            bbduk_list.append(bbduk)
            bbduk_sub = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_log.txt")
            bbduk_sub_list.append(bbduk_sub)
            cov_contigs = os.path.join("contig_abundance", readset.sample.name, readset.name + ".cov")
            cov_list_contigs.append(cov_contigs)
            cov_genes = os.path.join("gene_abundance", readset.sample.name, readset.name + ".cov")
            cov_list_genes.append(cov_genes)
    
            outdir = os.path.join("contig_abundance", readset.sample.name)
            outdir_genes = os.path.join("gene_abundance", readset.sample.name)
            if not os.path.exists(outdir):
                os.makedirs(os.path.join(outdir))
            if not os.path.exists(outdir_genes):
                os.makedirs(os.path.join(outdir_genes))
       
            if( readset.run_type == "PAIRED_END"):
                #infile = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
                flag = "0x2"
                if isinstance(ref_genome, str) and ref_genome != "":
                    out1 = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_R1.fastq.gz")
                    out2 = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_sub_R2.fastq.gz")
                else:
                    out1 = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
                    out2 = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")

                ## map against contigs
                if config.param('DEFAULT', 'mapper', 1, "string") == "bwa":
                    job = shotgun_metagenomics.bwa_mem_samtools(
                        reference_contigs,
                        bwt_contigs,
                        out1,
                        out2,
                        bam_contigs
                    )
                    job.name = "bwa_mem-contigs-" + readset.sample.name
                    job.subname = "bwa"
                    jobs.append(job)
                    
                elif config.param('DEFAULT', 'mapper', 1, "string") == "bbmap":
                    job = shotgun_metagenomics.bbmap_paired(
                        out1,
                        out2,
                        reference_contigs,
                        bam_contigs,
                        os.path.join("assembly", "bbmap_ref")
                    )
                    job.name = "bbmap-contigs-" + readset.sample.name
                    job.subname = "bbmap"
                    jobs.append(job)
                else:
                    raise Exception("Error: mapper should be set to either bbmap or bwa ")

                job = shotgun_metagenomics.flagstats(
                    bam_contigs,
                    flagstats_contigs
                )
                job.name = "flagstats-" + readset.sample.name
                job.subname = "flagstats"
                jobs.append(job)
                
            elif(readset.run_type == "SINGLE_END"):

                infile = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam.fastq.gz")
                flag = "0x0"
                
                if config.param("DEFAULT", "assembler", 1, "string") == "megahit":
                
                    job = shotgun_metagenomics.bwa_mem_samtools_se(
                        reference_contigs,
                        bwt_contigs,
                        infile,
                        bam_contigs
                    )
                    job.name = "bwa_mem-contigs-" + readset.sample.name
                    job.subname = "bwa"
                    jobs.append(job)
                
                elif config.param("DEFAULT", "assembler", 1, "string") == "miniasm" or config.param("DEFAULT", "assembler", 1, "string") == "canu":
                    job = shotgun_metagenomics.minimap2_samtools_se(
                        reference_contigs,
                        #bwt_contigs,
                        infile,
                        bam_contigs
                    )
                    job.name = "bwa_mem-contigs-" + readset.sample.name
                    job.subname = "bwa"
                    jobs.append(job)
                
                #-------------------------------------#
                job = shotgun_metagenomics.flagstats(
                    bam_contigs,
                    flagstats_contigs
                )
                job.name = "flagstats-" + readset.sample.name
                job.subname = "flagstats"
                jobs.append(job)

            job = shotgun_metagenomics.coverage_bed_v2_24(
                bam_contigs,
                bed_contigs,
                cov_contigs,
                flag
            )
            job.name = "bedtoolsCov-contigs-" + readset.sample.name
            job.subname = "bedtools"
            jobs.append(job)
            
            job = shotgun_metagenomics.coverage_bed_v2_24(
                bam_contigs,
                bed_genes,
                cov_genes,
                flag
            )
            job.name = "bedtoolsCov-genes-" + readset.sample.name
            job.subname = "bedtools"
            jobs.append(job)
            
        # Once all coverage has been computed, merge all tables.
        job = shotgun_metagenomics.merge_counts(
            cov_list_contigs,
            os.path.join("contig_abundance", "merged_contig_abundance.tsv"),
            "contigs"
        )
        job.name = "merge_contig_abundance"
        job.subname = "merge_abundance"
        jobs.append(job)
        
        job = shotgun_metagenomics.normalize_counts(
            os.path.join("contig_abundance", "merged_contig_abundance.tsv"),
            os.path.join("contig_abundance", "merged_contig_abundance_cpm.tsv")
        )
        job.name = "normalize_contig_abundance"
        job.subname = "normalization"
        jobs.append(job)
        
        job = shotgun_metagenomics.merge_counts(
            cov_list_genes,
            os.path.join("gene_abundance", "merged_gene_abundance.tsv"),
            "genes"
        )
        job.name = "merge_gene_abundance"
        job.subname = "merge_abundance"
        jobs.append(job)
        
        job = shotgun_metagenomics.normalize_counts(
            os.path.join("gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join("gene_abundance", "merged_gene_abundance_cpm.tsv")
        )
        job.name = "normalize_gene_abundance"
        job.subname = "normalization"
        jobs.append(job)
        
        if isinstance(ref_genome, str) and ref_genome != "":
            job = shotgun_metagenomics.merge_flagstats(
                flagstats_contigs_list,
                trimmomatic_list,
                bbduk_list,
                os.path.join("contig_abundance", "qc_mapping_stats.tsv"),
                bbduk_sub_list
            )
            job.name = "flagstats_merge"
            job.subname = "flagstats"
            jobs.append(job)
        else:
            job = shotgun_metagenomics.merge_flagstats(
                flagstats_contigs_list,
                trimmomatic_list,
                bbduk_list,
                os.path.join("contig_abundance", "qc_mapping_stats.tsv"),
                False
            )
            job.name = "flagstats_merge"
            job.subname = "flagstats"
            jobs.append(job)

        #Generate files for potential gam-ngs merge
        job = shotgun_metagenomics.get_insert_size(
            os.path.join("contig_abundance"),
            os.path.join("contig_abundance", "lib_stats.tsv"),
            os.path.join("contig_abundance", "qc_mapping_stats.tsv")
        )
        job.name = "get_insert_size"
        job.subname = "get_insert_size"
        jobs.append(job)
         
        return jobs 
     
    def exonerate(self):
        
        """
        Step exonerate(): Split co-assembly and predicted genes fasta files in smaller chunks to perform downstream annotations.
        """
        jobs = []

        infile_contigs = os.path.join("assembly", "Contigs.fasta")
        chunks_dir_contigs = os.path.join("assembly", "fasta_chunks")
        number_of_chunks_file_contigs = os.path.join("assembly", "estimated_number_of_chunks_contigs.txt")
        infile_faa = os.path.join("gene_prediction", "Contigs_renamed.faa")
        chunks_dir = os.path.join("gene_prediction", "fasta_chunks")
        number_of_chunks_file_genes = os.path.join("gene_prediction", "estimated_number_of_chunks_genes.txt")
        
        #FNA contigs assembly (i.e. contigs)
        job = shotgun_metagenomics.estimate_number_of_chunks(
            infile_contigs,
            number_of_chunks_file_contigs,
            str(config.param('exonerate', 'targeted_chunk_file_size_contigs', 1, 'int'))
        )
        job.name = "estimate_chunks_file_size"
        jobs.append(job)

        job = shotgun_metagenomics.exonerate(
            infile_contigs,
            chunks_dir_contigs,
            number_of_chunks_file_contigs,
            "Contigs.fasta" 
        )
        job.name = "exonerate_contigs_fna"
        job.subname = "exonerate"
        jobs.append(job)
        
        #FAA genes (i.e. predicted genes)
        job = shotgun_metagenomics.estimate_number_of_chunks(
            infile_faa,
            number_of_chunks_file_genes,
            str(config.param('exonerate', 'targeted_chunk_file_size_genes', 1, 'int'))
        )
        job.name = "estimate_chunks_file_size"
        jobs.append(job)
        
        job = shotgun_metagenomics.exonerate(
            infile_faa,
            chunks_dir,
            number_of_chunks_file_genes,
            "Contigs_renamed.faa" 
        )
        job.name = "exonerate_genes_fna"
        job.subname = "exonerate"
        jobs.append(job)
        
        return jobs    
    
    def kegg_annotation(self):
        """
        Step kegg_annotation(): Assignment of a KEGG ortholog (KO) to each gene.
                                One of the two following methods can be used:
                                1) kofamscan (free resource from GenomeNet) - Should be the default method in the .ini file.
                                   This method was modified however, because of the astronomical numbers of tmp files it generated.
                                2) hmmsearch against the KOfam concatenated hmm models.
                                3) DIAMOND BLASTp of predicted genes vs KEGG db.
        """

        jobs = []
        
        if config.param('DEFAULT', 'KO_method', 'string') == "kofamscan":
            
            chunks_dir = os.path.join("gene_prediction", "fasta_chunks")
            blast_dir = os.path.join("annotations", "kofamscan")
            number_chunks_file = os.path.join("gene_prediction", "estimated_number_of_chunks_genes.txt")
            if not os.path.exists(os.path.join("annotations", "kofamscan")):
                os.makedirs(os.path.join("annotations", "kofamscan"))
            
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
            
            num_chunks = num_chunks 
            for i in range(0, num_chunks):
                infiles.append(os.path.join(chunks_dir, "Contigs_renamed.faa_chunk_{:07d}".format(i)))
                outfiles.append(os.path.join(blast_dir, "kofamscan_chunk_{:07d}.tsv".format(i)))
                dones.append(os.path.join(blast_dir, "kofamscan_chunk_{:07d}.done".format(i)))
         
            job = shotgun_metagenomics.kofamscan_array_job(
                os.path.join(chunks_dir),
                "Contigs_renamed.faa_chunk_",
                os.path.join(blast_dir),
                "kofamscan_chunk_",
                "kofamscan",
                infiles,
                outfiles,
                dones,
                self._curr_scheduler
            )
            job.name = "kofamscan_array_job"
            job.subname = "kofamscan"
            job.job_array_num_task = num_chunks
            jobs.append(job)
            
            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks(
                blast_dir,
                os.path.join("annotations", "kofamscan.tsv"),
                num_chunks,
                "kofamscan"
            )
            job.name = "merge_kofamscan_kegg"
            job.subname = "merge"
            jobs.append(job)
            
            # Generate a clean table of module/pathways.
            job = shotgun_metagenomics.parse_kofam(
                os.path.join("annotations", "kofamscan.tsv"),
                os.path.join("annotations", "KOs_parsed.tsv")
            )
            job.name = "parse_kegg"
            job.subname = "merge"
            jobs.append(job)
        
        elif config.param('DEFAULT', 'KO_method', 'string') == "diamond_blastp":
            # Do blastp on KEGG for co-assembly 
            chunks_dir = os.path.join("gene_prediction", "fasta_chunks")
            blast_dir = os.path.join("annotations", "blastp_kegg")
            number_chunks_file = os.path.join("gene_prediction", "estimated_number_of_chunks_genes.txt")
            
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
            
            num_chunks = num_chunks 
            for i in range(0, num_chunks):
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
                config.param('blastp_kegg', 'db_kegg', 1, 'filepath'),
                dones,
                self._curr_scheduler
            )
            job.name = "DIAMOND-blastp_array_kegg"
            job.subname = "blastp_kegg"
            job.job_array_num_task = num_chunks
            jobs.append(job)
            
            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks(
                blast_dir,
                os.path.join("annotations", "blastp_kegg.tsv"),
                num_chunks,
                "blastp"
            )
            job.name = "merge_DIAMOND-blastp_kegg"
            job.subname = "merge"
            jobs.append(job)
            
            # Generate a clean table of module/pathways.
            job = shotgun_metagenomics.parse_kegg(
                os.path.join("annotations", "blastp_kegg.tsv"),
                os.path.join("annotations", "KOs_parsed.tsv")
            )
            job.name = "parse_kegg"
            job.subname = "merge"
            jobs.append(job)
        
        elif config.param('DEFAULT', 'KO_method', 'string') == "hmmsearch_kofam":
            # Do blastp on KEGG for co-assembly 
            chunks_dir = os.path.join("gene_prediction", "fasta_chunks")
            hmmsearch_outdir = os.path.join("annotations", "hmmsearch_kofam")
            number_chunks_file = os.path.join("gene_prediction", "estimated_number_of_chunks_genes.txt")
            if not os.path.exists(os.path.join("annotations", "hmmsearch_kofam")):
                os.makedirs(os.path.join("annotations", "hmmsearch_kofam"))
            
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
                #infiles.append(os.path.join(chunks_dir, "Contigs_renamed.faa_chunk_{:07d}".format(i)))
                #outfiles.append(os.path.join(blast_dir, "hmmsearch_chunk_{:07d}.tsv".format(i)))
                #dones.append(os.path.join(blast_dir, "hmmsearch_chunk_{:07d}.done".format(i)))
                infiles.append(os.path.join(chunks_dir, "Contigs_renamed.faa_chunk_{:07d}".format(i)))
                tblouts.append(os.path.join(hmmsearch_outdir, "hmmsearch_chunk_{:07d}.tblout".format(i)))
                domtblouts.append(os.path.join(hmmsearch_outdir, "hmmsearch_chunk_{:07d}.domtblout".format(i)))
                pfamtblouts.append(os.path.join(hmmsearch_outdir, "hmmsearch_chunk_{:07d}.pfamtblout".format(i)))
                dones.append(os.path.join(hmmsearch_outdir, "hmmsearch_chunk_{:07d}.done".format(i)))
         
            #job = shotgun_metagenomics.diamond_hmmsearch_kofam_array_job(
            #    os.path.join(chunks_dir),
            #    "Contigs_renamed.faa_chunk_",
            #    os.path.join(blast_dir),
            #    "hmmsearch_chunk_",
            #    "blastp",
            #    infiles,
            #    outfiles,
            #    config.param('hmmsearch_kofam', 'db_kegg', 1, 'filepath'),
            #    dones,
            #    self._curr_scheduler
            #)
            #job.name = "hmmsearch_array_kofam"
            #job.subname = "hmmsearch"
            #job.job_array_num_task = num_chunks
            #jobs.append(job)
            
            # Merge output chunks
            #job = shotgun_metagenomics.merge_chunks(
            #    blast_dir,
            #    os.path.join("annotations", "hmmsearch_kofam.tsv"),
            #    num_chunks,
            #    "blastp"
            #)
            #job.name = "merge_hmmsearch_kofam"
            #job.subname = "merge"
            #jobs.append(job)
        
            job = shotgun_metagenomics.hmmsearch_array_job(
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
                os.path.join("annotations"),
                num_chunks,
                "hmmsearch",
                "hmmsearch_kofam"
            )
            job.name = "hmmsearch_kofam_merge"
            job.subname = "merge"      
            jobs.append(job)
            
            # Generate a clean table of module/pathways.
            #job = shotgun_metagenomics.parse_kegg(
            #    os.path.join("annotations", "hmmsearch_kofam_tblout.tsv"),
            #    os.path.join("annotations", "KOs_parsed.tsv")
            #)
            
            job = shotgun_metagenomics.parse_kofam(
                os.path.join("annotations", "hmmsearch_kofam_tblout.tsv"),
                os.path.join("annotations", "KOs_parsed.tsv"),
                hmmsearch=True
            )
            job.name = "parse_kegg"
            job.subname = "merge"
            jobs.append(job)
        
        else:
            raise Exception("Error: In the .ini file, please chose between one of the following two values under the [DEFAULT] section: KO_method=kofamscan OR KO_method=diamond_blastp OR KO_method=hmmsearch_kofam")

        return jobs 
    
    def ncrna(self):
        """
        Step ncrna(): scan for non-coding RNA in contigs: rRNA and tRNA sequence.
        """

        jobs = []
        
        # Do tRNA of co-assembly 
        chunks_dir = os.path.join("assembly", "fasta_chunks")
        outdir = os.path.join("annotations", "trna")
        number_chunks_file = os.path.join("assembly", "estimated_number_of_chunks_contigs.txt")
        
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
        
        num_chunks = num_chunks 
        for i in range(0, num_chunks):
            #sys.stderr.write("loop i: " + str(i) + "\n")
            infiles.append(os.path.join(chunks_dir, "Contigs.fasta_chunk_{:07d}".format(i)))
            outfiles.append(os.path.join(outdir, "trna_chunk_{:07d}.tsv".format(i)))
            dones.append(os.path.join(outdir, "trna_chunk_{:07d}.done".format(i)))
     
        #job = shotgun_metagenomics.trnascanse_array_job(
        #    os.path.join(chunks_dir),
        #    "Contigs.fasta_chunk_",
        #    os.path.join(outdir),
        #    "trna_chunk_",
        #    "trnascanse",
        #    infiles,
        #    outfiles,
        #    dones,
        #    self._curr_scheduler
        #)
        #job.name = "trna_array_job"
        #job.subname = "trnascanse"
        #job.job_array_num_task = num_chunks
        #jobs.append(job)
        
        # Merge output chunks
        #job = shotgun_metagenomics.merge_chunks(
        #    outdir,
        #    os.path.join("annotations", "trna.tsv"),
        #    num_chunks,
        #    "trna"
        #)
        #job.name = "merge_tRNAScan-SE"
        #job.subname = "merge"
        #jobs.append(job)

        # Do rRNA prediction from co-assembly with barrnap
        prefixes = ["bac", "arc", "euk", "mito"]
        for prefix in prefixes:
            job = shotgun_metagenomics.barrnap(
                os.path.join("assembly", "Contigs.fasta"),
                prefix,
                os.path.join("annotations", "rrna", "barrnap", prefix + ".fna"),
                os.path.join("annotations", "rrna", "barrnap", prefix + ".tsv")
            )
            job.name = "barrnap_" + prefix
            job.subname = "barrnap"
            jobs.append(job)
       
            job = shotgun_metagenomics.split_barrnap_fna(
                os.path.join("annotations", "rrna", "barrnap", prefix + ".fna"),
                prefix,
                os.path.join("annotations", "rrna", "barrnap", prefix + "_5S.fna"), 
                os.path.join("annotations", "rrna", "barrnap", prefix + "_SSU.fna"), 
                os.path.join("annotations", "rrna", "barrnap", prefix + "_LSU.fna"), 
                os.path.join("annotations", "rrna", "barrnap", prefix + "_12S.fna"), 
                os.path.join("annotations", "rrna", "barrnap", prefix + "_5-8S.fna") 
            )
            job.name = "split_barrnap_fna_" + prefix
            jobs.append(job)
    
        prefixes = ["arc_SSU", "arc_LSU", "bac_SSU", "bac_LSU", "euk_SSU", "euk_LSU"]
        for prefix in prefixes:
            # create bed files for abundance - we have to only get the rea
            # Just take taxonomy file which contains coord for each rRNA genes, then generate bed file.
            job = shotgun_metagenomics.barrnap_fna_to_bed(
                os.path.join("annotations", "rrna", "barrnap", prefix + ".fna"),
                os.path.join("annotations", "rrna", "barrnap", prefix + ".bed")
            )
            job.name = "barrnap_to_bed_" + prefix
            job.subname = "barrnap_to_bed"
            jobs.append(job)

            # Then from this new bed file, process all bams to just get reads that falls into coords of new bed file
            # which corresponds to rRNA genes. Abundance will be used for both blastn and rdp taxonomy.
            cov_list = []
            for readset in self.readsets:
            
                if(readset.run_type == "SINGLE_END"):
                    flag = "0x0"
                elif(readset.run_type == "PAIRED_END"):
                    flag = "0x2"

                bam_contigs = os.path.join("contig_abundance", readset.sample.name, readset.name + ".bam")
            
                job = shotgun_metagenomics.coverage_bed_v2_24(
                    bam_contigs,
                    os.path.join("annotations", "rrna", "barrnap", prefix + ".bed"),
                    os.path.join("annotations", "rrna", "abundance", prefix + "_" + readset.name + ".cov"),
                    flag
                    #prefix
                )
                job.name = "bedtoolsCov-contigs-rrna_" + readset.sample.name + "_" + prefix
                job.subname = "bedtools_rrna"
                jobs.append(job)
                
                job = shotgun_metagenomics.rename_rrna_bed(
                    os.path.join("annotations", "rrna", "abundance", prefix + "_" + readset.name + ".cov"),
                    os.path.join("annotations", "rrna", "abundance", prefix + "_" + readset.name + "_renamed.cov"),
                )
                job.name = "rename_bed_" + readset.sample.name + "_" + prefix
                job.subname = "rename_bed"
                jobs.append(job)

                cov_list.append(os.path.join("annotations", "rrna", "abundance", prefix + "_" + readset.name + "_renamed.cov"))
        
            job = shotgun_metagenomics.merge_counts(
                cov_list,
                os.path.join("annotations", "rrna", "abundance", prefix + "_merged_abundance.tsv"),
                os.path.join("annotations", "rrna", "abundance", prefix + "_merged_abundance_cpm.tsv"),
                "contigs",
                True
            )
            job.name = "merge_contig_abundance_" + prefix
            job.subname = "merge_gene_abundance"
            jobs.append(job)
            
            # RDP classifier with rrna results.
            job = microbial_ecology.rdp_wrapper(
                os.path.join("annotations", "rrna", "barrnap", prefix + ".fna"),
                os.path.join("annotations", "rrna", "rdp", prefix + "_rdp.tsv")
            )
            job.name = "classify_" + prefix
            job.subname = "RDP"
            jobs.append(job)
           
            # Convert RDP table to in-house taxonomy format.
            job = shotgun_metagenomics.rdp_to_taxonomy(
                os.path.join("annotations", "rrna", "rdp", prefix + "_rdp.tsv"),
                os.path.join("annotations", "rrna", "rdp", prefix + "_rdp_taxonomy.tsv")
            )
            job.name = "rdp_to_taxonomy_rrna_" + prefix
            job.subname = "rdp_to_taxonomy"
            jobs.append(job)

            # Then with rdp output, generate otu table
            job = shotgun_metagenomics.generate_feature_table_rrna(
                os.path.join("annotations", "rrna", "rdp", prefix + "_rdp_taxonomy.tsv"),
                os.path.join("annotations", "rrna", "abundance", prefix + "_merged_abundance.tsv"),
                os.path.join("annotations", "rrna", "rdp", prefix +  "_feature_table.tsv")
            )
            job.name = "generate_feature_table_" + prefix
            job.subname = "generate_feature_table"
            jobs.append(job)

        return jobs
    
    def rpsblast_cog(self):
        """
        Step rpsblast_cog(): RPSBLAST of predicted genes vs COG db
        """
        
        jobs = []
        
        # Do rpsblast on COG for big assembly  
        if(config.param('DEFAULT', 'skip_cog', 1, 'string') == 'yes'):
            fname = os.path.join("annotations", "rpsblast_cog.tsv")
            open(fname, 'a').close()
            os.utime(fname, None)
        else:
            number_chunks_file = os.path.join("gene_prediction", "estimated_number_of_chunks_genes.txt")
        
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
            
            chunks_dir = os.path.join("gene_prediction", "fasta_chunks")
            blast_dir = os.path.join("annotations", "rpsblast_cog")
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
                config.param('cog', 'db', 1, 'string'),
                self._curr_scheduler
            )
            job.name = "rpsblast_cog"
            job.subname = "rpsblast"
            job.job_array_num_task = num_chunks
            jobs.append(job)
    
            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks(
                blast_dir,
                os.path.join("annotations", "rpsblast_cog.tsv"),
                num_chunks,
                "rpsblast"
            )
            job.name = "merge_rpsblast_cog"
            job.subname = "merge"
            jobs.append(job)
        
        return jobs 
    
    def rpsblast_kog(self):
        """
        Step rpsblast_kog(): RPSBLAST of predicted genes vs KOG db
        """
        
        jobs = []
        
        # Do rpsblast on KOG?
        if(config.param('DEFAULT', 'skip_kog', 1, 'string') == 'yes'):
            fname = os.path.join("annotations", "rpsblast_kog.tsv")
            open(fname, 'a').close()
            os.utime(fname, None)
        else:
            number_chunks_file = os.path.join("gene_prediction", "estimated_number_of_chunks_genes.txt")
        
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
            
            chunks_dir = os.path.join("gene_prediction", "fasta_chunks")
            blast_dir = os.path.join("annotations", "rpsblast_kog")
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
                config.param('kog', 'db', 1, 'string'),
                self._curr_scheduler
            )
            job.name = "rpsblast_kog"
            job.subname = "rpsblast"
            job.job_array_num_task = num_chunks
            jobs.append(job)
    
            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks(
                blast_dir,
                os.path.join("annotations", "rpsblast_kog.tsv"),
                num_chunks,
                "rpsblast"
            )
            job.name = "merge_rpsblast_kog"
            job.subname = "merge"
            jobs.append(job)
        
        return jobs 

    def hmmsearch_cazy(self):
        
        """
        Step hmmsearch_cazy(): HMMSCAN of predicted genes vs CAZy       
        """
        
        jobs = []
        if(config.param('DEFAULT', 'skip_cazy', 1, 'string') == 'yes'):
            fname = os.path.join("annotations", "cazy_parsed.tsv")
            open(fname, 'a').close()
            os.utime(fname, None)
        else:
            chunks_dir = os.path.join("gene_prediction", "fasta_chunks")
            hmmsearch_out_dir = os.path.join("annotations", "hmmsearch_cazy")
            number_chunks_file = os.path.join("gene_prediction", "estimated_number_of_chunks_genes.txt")
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
                tblouts.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.tblout".format(i)))
                domtblouts.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.domtblout".format(i)))
                pfamtblouts.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.pfamtblout".format(i)))
                dones.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.done".format(i)))

            job = shotgun_metagenomics.hmmsearch_array_job(
                os.path.join(chunks_dir),
                "Contigs_renamed.faa_chunk_",
                os.path.join(hmmsearch_out_dir),
                "hmmsearch_chunk_",
                infiles,
                tblouts, domtblouts, pfamtblouts,
                dones,
                config.param('cazy', 'db', required=True),
                self._curr_scheduler
            )
            job.name = "hmmsearch_cazy_array_job"
            job.subname = "hmmsearch"
            job.job_array_num_task = num_chunks
            jobs.append(job)

            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks_hmms(
                hmmsearch_out_dir,
                os.path.join("annotations"),
                num_chunks,
                "hmmsearch",
                "hmmsearch_cazy"
            )
            job.name = "hmmsearch_cazy_merge"
            job.subname = "merge"      
            jobs.append(job)
           
            job = shotgun_metagenomics.parse_hmms(
                os.path.join("annotations", "hmmsearch_cazy_tblout.tsv"),
                os.path.join("annotations", "hmmsearch_cazy_tblout_parsed.tsv")
            )
            job.name = "hmmsearch_cazy_filter"
            job.subname = "merge"      
            jobs.append(job)
            
            job = shotgun_metagenomics.parse_cazy(
                os.path.join("annotations", "hmmsearch_cazy_tblout_parsed.tsv"),
                config.param('cazy', 'ref_db', 1, type="filepath"),
                os.path.join("annotations", "cazy_parsed.tsv")
            )
            job.name = "hmmsearch_cazy_parse"
            job.subname = "parse"      
            jobs.append(job)
        
        return jobs

    def hmmsearch_pfam(self):
        
        """
        Step hmmsearch_pfam(): HMMSCAN of predicted genes vs PFAM-A DB.
        """
        
        jobs = []
        
        chunks_dir = os.path.join("gene_prediction", "fasta_chunks")
        hmmsearch_out_dir = os.path.join("annotations", "hmmsearch_pfam")
        number_chunks_file = os.path.join("gene_prediction", "estimated_number_of_chunks_genes.txt")
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
            tblouts.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.tblout".format(i)))
            domtblouts.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.domtblout".format(i)))
            pfamtblouts.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.pfamtblout".format(i)))
            dones.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.done".format(i)))

        job = shotgun_metagenomics.hmmsearch_array_job(
            os.path.join(chunks_dir),
            "Contigs_renamed.faa_chunk_",
            os.path.join(hmmsearch_out_dir),
            "hmmsearch_chunk_",
            infiles,
            tblouts, domtblouts, pfamtblouts,
            dones,
            config.param('pfam', 'db', required=True),
            self._curr_scheduler
        )
        job.name = "hmmsearch_pfam_array_job"
        job.subname = "hmmsearch"
        job.job_array_num_task = num_chunks
        jobs.append(job)

        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks_hmms(
            hmmsearch_out_dir,
            os.path.join("annotations"),
            num_chunks,
            "hmmsearch",
            "hmmsearch_pfam"
        )
        job.name = "hmmsearch_pfam_merge"
        job.subname = "merge"      
        jobs.append(job)
       
        job = shotgun_metagenomics.parse_hmms(
            os.path.join("annotations", "hmmsearch_pfam_tblout.tsv"),
            os.path.join("annotations", "hmmsearch_pfam_tblout_parsed.tsv")
        )
        job.name = "hmmsearch_pfam_parse"
        job.subname = "merge"      
        jobs.append(job)
        
        return jobs

    def diamond_blastp_nr(self):
        
        """
        Step diamond_blastp_nr(): DIAMOND BLASTp of predicted genes vs NCBI nr
        """
        
        jobs = []
       
        fname = os.path.join("annotations", "blastp_nr_annotated.tsv")
        open(fname, 'a').close()
        os.utime(fname, None)
        
        # Do blastp on KEGG for big assembly  
        chunks_dir = os.path.join("gene_prediction", "fasta_chunks")
        blast_dir = os.path.join("annotations", "blastp_nr")
        number_chunks_file = os.path.join("gene_prediction", "estimated_number_of_chunks_genes.txt")
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
            config.param('blastp_nr', 'db_nr', 1, 'filepath'),
            dones,
            self._curr_scheduler
        )
        job.name = "DIAMOND_blastp_nr_array_job"
        job.subname = "blastp_nr"
        job.job_array_num_task = num_chunks
        jobs.append(job)
        
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join("annotations", "blastp_nr.tsv"),
            num_chunks,
            "blastp"
        )
        job.name = "merge_DIAMOND_blastp_nr"
        job.subname = "merge"
        jobs.append(job)
        
        return jobs 

    
    def blastn_nt_contigs(self):
        """
        Step blastn_nt_contigs(): BLASTn contigs against NCBI nt
        """
        mkdir_p(os.path.join(self._root_dir, "annotations", "blastn_nt_contigs"))
        
        jobs = []
        
        # Do blastn on nt for big assembly  
        chunks_dir = os.path.join("assembly", "fasta_chunks")
        blast_dir = os.path.join("annotations", "blastn_nt_contigs")
        number_chunks_file = os.path.join("assembly", "estimated_number_of_chunks_contigs.txt")
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
            dones
        )
        job.name = "blastn_nt_contigs"
        job.subname = "blastn"
        job.job_array_num_task = num_chunks
        jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join("annotations", "blastn_nt_contigs.tsv"),
            num_chunks,
            "blastn" 
        )
        job.name = "blastn_nt_big_contigs_merge"
        job.subname = "blastn"
        jobs.append(job)
        
        return jobs
    
    def blastn_ncbi_genomes(self):
        """
        Step blastn_nt_genomes(): BLASTn genomes against NCBI genomes
        """
        jobs = []
        
        # Do blastn on nt for big assembly  
        chunks_dir = os.path.join("assembly", "fasta_chunks")
        blast_dir = os.path.join("annotations", "blastn_ncbi_genomes")
        number_chunks_file = os.path.join("assembly", "estimated_number_of_chunks_contigs.txt")
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
            dones
        )
        job.name = "blastn_nt_contigs_array_job"
        job.subname = "blastn"
        job.job_array_num_task = num_chunks
        jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join("annotations", "blastn_nt_contigs.tsv"),
            num_chunks,
            "blastn" 
        )
        job.name = "blastn_nt_contigs_merge"
        job.subname = "blastn"
        jobs.append(job)
        
        return jobs
    
    #Deprecated. Still there for reference.
    def hmmsearch_rrna(self):
        jobs = []
        
        rrna_dir = os.path.join("annotations", "rrna")
        job = shotgun_metagenomics.hmmscan(
            os.path.join("assembly", "Contigs.fasta"),
            os.path.join("annotations", "rrna", "rrna"),
            rrna_dir,
            config.param('rrna', 'db', 1, "filepath") 
        )
        job.name = "hmmsearch_rrna"
        job.subname = "hmmsearch_rrna"
        jobs.append(job)
        
        job = shotgun_metagenomics.split_rrna(
            os.path.join("annotations", "rrna", "rrna.domtblout"),
            os.path.join("assembly", "Contigs.fasta"),
            os.path.join("annotations", "rrna", "arc_5S.fasta"),
            os.path.join("annotations", "rrna", "bac_5S.fasta"),
            os.path.join("annotations", "rrna", "euk_5S.fasta"),
            os.path.join("annotations", "rrna", "arc_16S.fasta"),
            os.path.join("annotations", "rrna", "bac_16S.fasta"),
            os.path.join("annotations", "rrna", "euk_18S.fasta"),
            os.path.join("annotations", "rrna", "arc_23S.fasta"),
            os.path.join("annotations", "rrna", "bac_23S.fasta"),
            os.path.join("annotations", "rrna", "euk_28S.fasta"),
            os.path.join("annotations", "rrna", "arc_bac_16S.fasta")
        )
        job.name = "split_rrna"
        job.subname = "split_rrna"
        jobs.append(job)
        
        prefixes = ["arc_bac_16S", "bac_23S", "euk_28S", "euk_18S"]
            
        # get ncbi id with blastn against nt AND perform RDP classifier.
        # Note that blastn otu tables will be generated in a later step: taxonomy_annotation.
        for prefix in prefixes:
            # create bed files for abundance - we have to only get the rea
            # Just take taxonomy file which contains coord for each rRNA genes, then generate bed file.
            job = shotgun_metagenomics.rrna_to_bed(
                os.path.join("annotations", "rrna", prefix + ".fasta"),
                os.path.join("annotations", "rrna", prefix + ".bed")
            )
            job.name = "fasta_to_bed_" + prefix
            job.subname = "fasta_to_bed"
            jobs.append(job)

            # Then from this new bed file, process all bams to just get reads that falls into coords of new bed file
            # which corresponds to rRNA genes. Abundance will be used for both blastn and rdp taxonomy.
            cov_list = []
            for readset in self.readsets:
            
                if(readset.run_type == "SINGLE_END"):
                    flag = "0x0"
                elif(readset.run_type == "PAIRED_END"):
                    flag = "0x2"

                bam_contigs = os.path.join("contig_abundance", readset.sample.name, readset.name + ".bam")
            
                job = shotgun_metagenomics.coverage_bed(
                    bam_contigs,
                    os.path.join("annotations", "rrna", prefix + ".bed"),
                    os.path.join("annotations", "rrna", "abundance", prefix + "_" + readset.name + ".cov"),
                    flag
                )
                job.name = "bedtoolsCov-contigs-rrna_" + readset.sample.name + "_" + prefix
                job.subname = "bedtools_rrna"
                jobs.append(job)
                cov_list.append(os.path.join("annotations", "rrna", "abundance", prefix + "_" + readset.name + ".cov"))
        
            job = shotgun_metagenomics.merge_counts(
                cov_list,
                os.path.join("annotations", "rrna", prefix + "_merged_abundance.tsv"),
                os.path.join("annotations", "rrna", prefix + "_merged_abundance_cpm.tsv"),
                "contigs",
                True
            )
            job.name = "merge_gene_abundance_contigs_" + prefix
            job.subname = "merge_gene_abundance"
            jobs.append(job)
            
            # RDP classifier with rrna results.
            job = microbial_ecology.rdp_wrapper(
                os.path.join("annotations", "rrna", prefix + ".fasta"),
                os.path.join("annotations", "rrna", "rrna_rdp", prefix + "_rdp.tsv")
            )
            job.name = "classify_" + prefix
            job.subname = "RDP"
            jobs.append(job)
           
            # Convert RDP table to in-house taxonomy format.
            job = shotgun_metagenomics.rdp_to_taxonomy(
                os.path.join("annotations", "rrna", "rrna_rdp", prefix + "_rdp.tsv"),
                os.path.join("annotations", "rrna", "rrna_rdp", prefix + "_rdp_taxonomy.tsv")
            )
            job.name = "rdp_to_taxonomy_rrna_" + prefix
            job.subname = "rdp_to_taxonomy"
            jobs.append(job)

            # Then with rdp output, generate otu table
            job = shotgun_metagenomics.generate_feature_table(
                os.path.join("annotations", "rrna", "rrna_rdp", prefix + "_rdp_taxonomy.tsv"),
                os.path.join("annotations", "rrna", prefix + "_merged_abundance.tsv"),
                os.path.join("annotations", "rrna", "rrna_rdp", prefix +  "_feature_table.tsv")
            )
            job.name = "generate_feature_table"
            job.subname = "generate_feature_table"
            jobs.append(job)
                                
            prefixes = ["arc_bac_16S", "bac_23S", "euk_28S", "euk_18S"]

            for k in range(0, len(prefixes)):
                for m in range(1, 6):

                    job = microbial_ecology.summarize_taxonomy_absolute(
                        os.path.join("annotations", "rrna", "rrna_rdp", prefixes[k] + "_feature_table.tsv"),
                        m,
                        os.path.join("annotations", "rrna", "rrna_rdp", prefixes[k], "absolute") 
                    )
                    job.name = "tax_summary" + "_absolute_" + prefixes[k] + "_" + str(m)
                    job.subname = "summarize_taxonomy"
                    jobs.append(job)
                    
                    job = microbial_ecology.summarize_taxonomy_relative(
                        os.path.join("annotations", "rrna", "rrna_rdp", prefixes[k] + "_feature_table.tsv"),
                        m,
                        os.path.join("annotations", "rrna", "rrna_rdp", prefixes[k], "relative") 
                    )
                    job.name = "tax_summary" + "_relative_" + prefixes[k] + "_" + str(m)
                    job.subname = "summarize_taxonomy"
                    jobs.append(job)
            
            ## Then, once abundance is done, BLASTN on NCBI genomes. OTU table generation will be done later in taxonomic_annotation step!
            job = shotgun_metagenomics.blastn(
                os.path.join("annotations", "rrna", prefix + ".fasta"),
                os.path.join("annotations", "rrna", "rrna_blastn", "blastn_nt_" + prefix + ".tsv"),
                os.path.join("annotations", "rrna", "rrna_blastn"),
                "blastn_genomes"
            )
            job.name = "blastn_nt_" + prefix
            job.subname = "blastn_genomes"
            jobs.append(job)
            
            job = shotgun_metagenomics.extract_silva_taxonomy(
                os.path.join("annotations", "rrna", "rrna_blastn", "blastn_nt_" + prefix + ".tsv"),
                os.path.join("annotations", "rrna", prefix + "_merged_abundance.tsv"),
                config.param('silva_tax', 'accession_to_tax', required=True),
                os.path.join("annotations", "rrna", "rrna_blastn", "feature_table_" + prefix + ".tsv")
            )
            job.name = "extract_taxonomy_rrna_blastn_" + prefix
            job.subname = "silva_tax"
            jobs.append(job)

        return jobs
   
    def reads_centric_taxonomy(self):
        """
        Step reads_centric_taxonomy(): Perform classification based on jgi itagger paper approach for MG data.
        """

        jobs = []
        rdps = []
        names = []
        
        for readset in self.readsets:
            if(readset.run_type == "PAIRED_END"):
                in1p = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
                in2p = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
                joined = os.path.join("qced_reads", readset.sample.name, readset.name + ".joined.fasta.gz")
                joined_matched = os.path.join("qced_reads", readset.sample.name, readset.name + ".joined_rRNA.fasta")
            
            elif(readset.run_type == "SINGLE_END"):
                in1p = os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam.fastq.gz")
                joined = os.path.join("qced_reads", readset.sample.name, readset.name + ".fasta.gz")
                joined_matched = os.path.join("qced_reads", readset.sample.name, readset.name + "_rRNA.fasta")
                
            rdp = os.path.join("qced_reads", readset.sample.name, readset.name + ".joined_rRNA_rdp.tsv")
            log = os.path.join("qced_reads", readset.sample.name, readset.name + ".joined.duk_log.txt")
            rdps.append(rdp)
            names.append(readset.sample.name)

            if(readset.run_type == "PAIRED_END"):
                # join reads.
                job = shotgun_metagenomics.join_reads(
                    in1p,
                    in2p,
                    joined
                )
                job.name = "join_reads_" + readset.sample.name
                job.subname = "join_reads"
                jobs.append(job)
            
            elif(readset.run_type == "SINGLE_END"):
                job = shotgun_metagenomics.fastq_to_fasta(
                    in1p,
                    joined
                )
                job.name = "fastq_to_fasta_" + readset.sample.name
                job.subname = "fastq_to_fasta"
                jobs.append(job)

            # get potential rRNA reads with duk against 99% dereplicated greengenes db.
            job = shotgun_metagenomics.duk_gz_matched_only(
                joined,
                joined_matched,
                log,
                config.param('DB', 'rrna', 1, 'filepath') 
            )
            job.name = "duk_rrna_" + readset.sample.name
            job.subname = "duk_gz"
            jobs.append(job)

            # For each of these potential rRNA reads, perform rdp classification
            job = microbial_ecology.rdp_wrapper(
                joined_matched,
                rdp
            )
            job.name = "classify_" + readset.sample.name
            job.subname = "RDP"
            jobs.append(job)

        # with rdp results in hand, generate otu table
        job = shotgun_metagenomics.generate_feature_table_reads_centric(
            rdps,
            names,
            os.path.join("annotations", "taxonomy_reads_centric", "feature_table.tsv")
        )
        job.name = "read_centric_feature_table"
        job.subname = "generate_feature_table"
        jobs.append(job)
            
        infiles = [os.path.join("annotations", "taxonomy_reads_centric", "feature_table.tsv")]
        prefixes = ["reads_centric"]
        directories = [os.path.join("annotations", "taxonomy_reads_centric")]
        normalizations = [""]
        types = ["absolute", "relative"]
        organisms = [""]
        organisms2 = ["all"]
        
        for i in range(0, len(infiles)):
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
                                job.name = "tax_summary_" + normalizations[n] + "_" + types[j] + "_" + organisms2[k]  + "_L" + str(m) + "_" + prefixes[i]
                                job.subname = "summarize_taxonomy"
                                jobs.append(job)
                            elif types[j] == "relative":
                                job = microbial_ecology.summarize_taxonomy_relative(
                                    os.path.join(directories[i], "feature_table" + normalizations[n] + organisms[k] + ".tsv"),
                                    m,
                                    os.path.join(directories[i], organisms2[k], types[j]) 
                                )
                                job.name = "tax_summary_" + normalizations[n] + "_" + types[j] + "_" + organisms2[k]  + "_L" + str(m) + "_" + prefixes[i]
                                job.subname = "summarize_taxonomy"
                                jobs.append(job)
            
                
            # Plot taxa - ALL
            for n in range(0, len(normalizations)):
                for j in range(0, len(types)):
                    for k in range(0, len(organisms)):
                        for m in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
                            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                                os.path.join(directories[i], organisms2[k], types[j], "feature_table" + normalizations[n] + organisms[k] + "_L" + str(m) + ".txt"),
                                os.path.join(directories[i], organisms2[k], types[j], "plots"),
                                "taxonomy_L" + str(m) + "_" + types[j]
                            )
                            job.name = "plot_taxa_single_" + normalizations[n] + "_" + types[j]  + "_L" + str(m) + "_" + prefixes[i] + "_" + organisms2[k]
                            job.subname = "plot_taxa"
                            jobs.append(job)

        return jobs
    
    # Generate taxonomic summary from various otu tables generated with various methods.
    def taxonomic_annotation(self):
        """
        Step taxonomic_annotation(): Using DIAMOND blastp results generated in the diamond_blastp vs nr step, each contig
                                     will be assigned a taxonimic lineage using the CAT package. Taxonomic summaries and 
                                     figures are then generated.
        """
        
        prefixes = ["consensus"]
        for prefix in prefixes:
            mkdir_p(os.path.join("annotations/taxonomy_" + prefix + "/all/absolute/plots"))
            mkdir_p(os.path.join("annotations/taxonomy_" + prefix + "/all/relative/plots"))
            mkdir_p(os.path.join("annotations/taxonomy_" + prefix + "/others/absolute/plots"))
            mkdir_p(os.path.join("annotations/taxonomy_" + prefix + "/others/relative/plots"))
            mkdir_p(os.path.join("annotations/taxonomy_" + prefix + "/bacteriaArchaea/absolute/plots"))
            mkdir_p(os.path.join("annotations/taxonomy_" + prefix + "/bacteriaArchaea/relative/plots"))
            mkdir_p(os.path.join("annotations/taxonomy_" + prefix + "/normalized_bacteriaArchaea/relative/plots"))
            mkdir_p(os.path.join("annotations/taxonomy_" + prefix + "/normalized_bacteriaArchaea/absolute/plots"))
        
        jobs = []
        
        job = shotgun_metagenomics.convert_orf_ids_for_cat(
           os.path.join("gene_prediction", "Contigs_renamed.gff"),
           os.path.join("annotations", "blastp_nr.tsv"),
           os.path.join("annotations", "blastp_nr_alt_orf_ids.tsv")
        )
        job.name = "convert_ord_ids_for_CAT"
        job.subname = "convert_ord_ids_for_CAT"
        jobs.append(job)
        
        job = shotgun_metagenomics.CAT(
           os.path.join("assembly", "Contigs.fasta"),
           os.path.join("gene_prediction", "Contigs.faa"),
           os.path.join("annotations", "blastp_nr_alt_orf_ids.tsv"),
           os.path.join("annotations", "taxonomy_consensus", "out")
        )
        job.name = "CAT"
        job.subname = "CAT"
        jobs.append(job)
        
        job = shotgun_metagenomics.generate_feature_table_consensus(
           os.path.join("annotations", "taxonomy_consensus", "out.contig2classification_with_names.tsv"),
           os.path.join("contig_abundance", "merged_contig_abundance.tsv"),
           os.path.join("annotations", "taxonomy_consensus", "taxonomy.tsv"),
           os.path.join("annotations", "taxonomy_consensus", "feature_table.tsv")
        )
        job.name = "generate_feature_table_consensus"
        job.subname = "generate_feature_table"
        jobs.append(job)
   
        prefixes = [
            "consensus"
        ]
        infiles = [
            os.path.join("annotations", "taxonomy_consensus", "feature_table.tsv")
        ]
        directories = [
            os.path.join("annotations", "taxonomy_consensus")
        ]
        abundances = [
            os.path.join("contig_abundance", "merged_contig_abundance_cpm.tsv") 
        ] 
  
        # Then proces all otu tables including the extended ones if -e option
        for i in range(0, len(infiles)):
            # Generate tsv  table of raw counts (unnormalized)
            job = microbial_ecology.split_feature_table(
                os.path.join(directories[i], "feature_table.tsv"),
                os.path.join(directories[i], "feature_table_bacteriaArchaea.tsv"),
                os.path.join(directories[i], "feature_table_others.tsv"),
                "bacteriaArchaea"
            )
            job.name = "split_feature_table_bactArch_unnormalized_" + prefixes[i]
            job.subname = "split_feature_table"
            jobs.append(job)
    
            # This job doesnt really do any normalization at all. It takes cpm values and creates an otu table.
            job = microbial_ecology.normalize_cpm_feature_table(
                abundances[i],
                os.path.join(directories[i], "feature_table.tsv"),
                os.path.join(directories[i], "feature_table_normalized.tsv")
            )
            job.name = "normalize_feature_table_all_" + prefixes[i]
            job.subname = "normalization"
            jobs.append(job)
    
            job = microbial_ecology.split_feature_table(
                os.path.join(directories[i], "feature_table_normalized.tsv"),
                os.path.join(directories[i], "feature_table_normalized_bacteriaArchaea.tsv"),
                os.path.join(directories[i], "feature_table_normalized_others.tsv"),
                "bacteriaArchaea"
            )
            job.name = "split_feature_table_bactArch_" + prefixes[i]
            job.subname = "split_feature_table"
            jobs.append(job)
           
            normalizations = ["_normalized", ""]
            types = ["absolute", "relative"]
            organisms = ["", "_others", "_bacteriaArchaea"]
            organisms2 = ["all", "others", "bacteriaArchaea"]

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
                
            # Plot taxa - ALL
            for n in range(0, len(normalizations)):
                for j in range(0, len(types)):
                    for k in range(0, len(organisms)):
                        for m in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
                            job = microbial_ecology.plot_taxa_single_with_mapping_file(
                                os.path.join(directories[i], organisms2[k], types[j], "feature_table" + normalizations[n] + organisms[k] + "_L" + str(m) + ".tsv"),
                                os.path.join(directories[i], organisms2[k], types[j], "plots"),
                                "taxonomy_L" + str(m) + "_" + types[j]
                            )
                            job.name = "plot_taxa_single" + normalizations[n] + "_" + types[j]  + "_L" + str(m) + "_" + prefixes[i] + "_" + organisms2[k]
                            job.subname = "plot_taxa"
                            jobs.append(job)
                    
        return jobs


    # Here generate final GFF (for viewing data in a genome browser). 
    # And generate final DDA sheets. with logFC, gene_name and actual normalized values.
    def finalize(self):
        jobs = []
        job = shotgun_metagenomics.generate_gff(
            # infiles
            os.path.join("gene_prediction", "Contigs_renamed.gff"),
            os.path.join("assembly", "Contigs.fasta"),
            os.path.join("annotations", "KOs_parsed.tsv"),
            os.path.join("annotations", "hmmsearch_pfam_tblout.tsv"),
            os.path.join("annotations", "rpsblast_cog.tsv"),
            os.path.join("annotations", "rpsblast_kog.tsv"),
            os.path.join("annotations", "taxonomy_consensus", "taxonomy.tsv"),
            os.path.join("annotations", "blastp_nr_annotated.tsv"),
            os.path.join("annotations", "cazy_parsed.tsv"),
            # outfiles
            os.path.join("annotations", "Contigs.gff"),
            os.path.join("annotations", "Contigs.fasta"),
            os.path.join("annotations", "annotations.tsv")
        )
        job.name = "generate_gff_and_annotations_file"
        job.subname = "generate_gff"
        jobs.append(job)
        
        return jobs
    
    def differential_abundance(self):
        jobs = []
        mkdir_p(os.path.join(self._root_dir, "DDA"))
        
        # Do not necessarily do pairwise.
        if(config.param('DDA', 'do_deg_pairwise', 1, 'string') == 'yes'):
            design_files = config.param('DEFAULT', 'designFile', 1, 'string')
            designs = design_files.split(":")
            for design_file in designs:
                job = shotgun_metagenomics.edger(
                    os.path.join("gene_abundance", "merged_gene_abundance.tsv"),
                    os.path.join(design_file),
                    os.path.join("DDA")
                )
                job.name = "DDA"
                job.subname = "DDA"
                jobs.append(job)
        
        # Do not necessarily do GLMs.
        if config.param('DDA', 'do_deg_glm', 1, 'string') == 'yes':
            mapping_files = config.param('DDA', 'mapping_files', 1, 'string')
            mappings = mapping_files.split(":")
            for map in mappings:
                job = shotgun_metagenomics.edger_glm(
                    os.path.join("gene_abundance", "merged_gene_abundance.tsv"),
                    map,
                    os.path.join("DDA_GLM")
                )
                job.name = "DDA"
                job.subname = "DDA"
                jobs.append(job)
        

        return jobs
    
    def alpha_diversity(self):
        
        """
        Step alpha_diversity(): Compute alpha-diversity on contig abundance and gene abundance matrices.
        """
        
        #def mkdir_p(path):
        #    try:
        #        os.makedirs(path)
        #    except OSError as exc: # Python >2.5
        #        if exc.errno == errno.EEXIST and os.path.isdir(path):
        #            pass
        #        else: raise
        mkdir_p(os.path.join("alpha_div", "contig_abundance"))
        mkdir_p(os.path.join("alpha_div", "gene_abundance"))
        mkdir_p(os.path.join("alpha_div", "rpob"))
        mkdir_p(os.path.join("alpha_div", "reca"))
        jobs = []
        
        job = microbial_ecology.rtk(
            os.path.join("contig_abundance", "merged_contig_abundance.tsv"),
            os.path.join("alpha_div", "contig_abundance"),
            "contig_id",
            remove_last_col = False
        )
        job.name = "alpha_diversity_contig_abundance"
        job.subname = "alphadiv"
        jobs.append(job)
        
        job = microbial_ecology.rtk(
            os.path.join("gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join("alpha_div", "gene_abundance"),
            "gene_id",
            remove_last_col = False
        )
        job.name = "alpha_diversity_gene_abundance"
        job.subname = "alphadiv"
        jobs.append(job)

        job = microbial_ecology.generate_cog_abundance(
            os.path.join("gene_abundance", "merged_gene_abundance.tsv"),
            "COG0085",
            os.path.join("annotations", "rpsblast_cog.tsv"),
            os.path.join("alpha_div", "rpob", "rpob_abundance.tsv")
        )
        job.name = "generate_rpob_matrix_for_alpha_div"
        job.subname = "alphadiv"
        jobs.append(job)
        
        job = microbial_ecology.generate_cog_abundance(
            os.path.join("gene_abundance", "merged_gene_abundance.tsv"),
            "COG0468",
            os.path.join("annotations", "rpsblast_cog.tsv"),
            os.path.join("alpha_div", "reca", "reca_abundance.tsv")
        )
        job.name = "generate_reca_matrix_for_alpha_div"
        job.subname = "alphadiv"
        jobs.append(job)
        
        job = microbial_ecology.rtk(
            os.path.join("alpha_div", "rpob", "rpob_abundance.tsv"),
            os.path.join("alpha_div", "rpob"),
            "gene_id",
            remove_last_col = False
        )
        job.name = "alpha_diversity_rpob"
        job.subname = "alphadiv"
        jobs.append(job)
        
        job = microbial_ecology.rtk(
            os.path.join("alpha_div", "reca", "reca_abundance.tsv"),
            os.path.join("alpha_div", "reca"),
            "gene_id",
            remove_last_col = False
        )
        job.name = "alpha_diversity_reca"
        job.subname = "alphadiv"
        jobs.append(job)

        return jobs
    
    def beta_diversity(self):
        
        """
        Step beta_diversity(): Compute Bray beta-diversity + PCoA ordinations on normalized all + normalized bacteria and archaea only.
        """
        
        #mkdir_p(os.path.join("beta_div", "bray_curtis_feature_table_normalized"))
        #mkdir_p(os.path.join("beta_div", "contig_abundance_normalized_bacteriaArchaea"))
        #mkdir_p(os.path.join("beta_div", "contig_abundance_normalized"))
        #mkdir_p(os.path.join("beta_div", "gene_abundance_normalized"))
        
        # Beta diversity.
        #def mkdir_p(path):
        #    try:
        #        os.makedirs(path)
        #    except OSError as exc: # Python >2.5
        #        if exc.errno == errno.EEXIST and os.path.isdir(path):
        #            pass
        #        else: raise
        mkdir_p(os.path.join("beta_div"))
        fname = os.path.join("beta_div", "tree.fasttree")
        if not os.path.isfile(fname):
            open(fname, 'a').close()
            os.utime(fname, None)
            sys.stderr.write("tree file not does exists...\n")
        
        jobs = []
       
        # define files.
        #feature_table_consensus_tsv     = os.path.join("annotations", "taxonomy_consensus", "feature_table_normalized.tsv")
        feature_table_consensus_ba_tsv  = os.path.join("annotations", "taxonomy_consensus", "feature_table_normalized_bacteriaArchaea.tsv")
        contig_abundance = os.path.join("contig_abundance", "merged_contig_abundance_cpm.tsv")
        gene_abundance = os.path.join("gene_abundance", "merged_gene_abundance_cpm.tsv")

        types_tsv = [feature_table_consensus_ba_tsv, contig_abundance, gene_abundance]
        prefixes = ["bray_curtis_contig_bacteriaArchaea", "bray_curtis_contig_abundance", "bray_curtis_gene_abundance"]

        for i in range(len(types_tsv)):
            mkdir_p(self._root_dir + "beta_div/" + prefixes[i] + "/")

            curr_table_tsv = types_tsv[i]

            job = microbial_ecology.beta_diversity(
                curr_table_tsv,
                "bray-curtis",
                self._root_dir + "beta_div/" + prefixes[i] + "/",
                self._root_dir + "beta_div/tree.fasttree",
                prefixes[i]
            )
            job.name = "beta_diversity_" + prefixes[i]
            job.subname = "beta_diversity"
            jobs.append(job)
            
            job = microbial_ecology.principal_coordinates(
                self._root_dir + "beta_div/" + prefixes[i] + "/" + prefixes[i] + ".tsv",
                self._root_dir + "beta_div/" + prefixes[i] + "/" + prefixes[i] + "_coords.tsv"
            )
            job.name = "beta_diversity_coords_" + prefixes[i]
            job.subname = "pc"
            jobs.append(job)
            
            job = microbial_ecology.pca_3d_plot(
                self._root_dir + "beta_div/" + prefixes[i] + "/" + prefixes[i] + "_coords.tsv",
                self._root_dir + "beta_div/" + prefixes[i] + "/" + "3d_bray_curtis_plot"
            )
            job.name = "beta_diversity_3d_plot_" + prefixes[i]
            job.subname = "pca_plot"
            jobs.append(job)
            
        return jobs

    def binning(self):
        """
        step binning(): Generates MAGs/bins with the metabat2 or maxbin2 software.
        """
        
        jobs = []
        bams = []
        abundance_files = []
       
        binners = config.param('DEFAULT', 'binner', 1, 'string').split(":")

        for binner in binners:
            if binner == "metabat2" or binner == "maxbin2":
                #binner = "metabat2"
                if not os.path.exists(os.path.join("binning", binner)):
                    os.makedirs(os.path.join("binning", binner))
                if not os.path.exists(os.path.join("annotations", binner, "parsed_bins")):
                    os.makedirs(os.path.join("annotations", binner, "parsed_bins"))
            else:
                raise Exception("Error: mapper should be set to either metabat2, maxbin2 or metabat2:maxbin2 ")

        for readset in self.readsets:
            bam_contigs = os.path.join("contig_abundance", readset.sample.name, readset.name + ".bam")
            bams.append(bam_contigs)
            abundance_file = os.path.join("contig_abundance", readset.sample.name, readset.name + "_abundance.txt")
            abundance_files.append(abundance_file)
            
        job = shotgun_metagenomics.get_abundance_for_metabat(
            bams,
            os.path.join("binning","abundance.txt"),
            os.path.join("binning", "paired_contigs.txt")
        )
        job.name = "metabat2_abundance"
        job.subname = "metabat_abundance"
        jobs.append(job)
      
        for binner in binners:

            if binner == "metabat2":
                job = shotgun_metagenomics.metabat2(
                    os.path.join("assembly", "Contigs.fasta"),
                    os.path.join("binning", "abundance.txt"),
                    os.path.join("binning", binner, "out"),
                    os.path.join("binning", binner, "metabat2.done") 
                )
                job.name = "metabat2"
                job.subname = "metabat2"
                jobs.append(job)

            elif binner == "maxbin2":
                job = shotgun_metagenomics.maxbin2(
                    os.path.join("assembly", "Contigs.fasta"),
                    os.path.join("binning", "abundance.txt"),
                    os.path.join("binning", binner),
                    os.path.join("binning", binner, "out"),
                    os.path.join("binning", binner, "maxbin2.done") 
                )
                job.name = "maxbin2"
                job.subname = "maxbin2"
                jobs.append(job)
       
            # checkm on raw bins
            job = shotgun_metagenomics.checkm(
                os.path.join("binning", binner, binner + ".done"),
                os.path.join("binning", binner), 
                os.path.join("binning", binner, "out_checkm"),
                os.path.join("binning", binner, "out_checkm.txt")
            )
            job.name = "checkm_raw_bins_" + binner
            job.subname = "checkm"
            jobs.append(job)
           
            # parse bins by contigs taxonomy
            job = shotgun_metagenomics.parse_bins(
                os.path.join("binning", binner, "out_checkm.txt"),
                os.path.join("binning", binner), 
                os.path.join("binning", binner, "parsed_bins"),
                os.path.join("binning", binner, "parsed_bins.tsv"),
                os.path.join("binning", binner, "raw_bins.tsv"),
                os.path.join("annotations", "taxonomy_consensus", "taxonomy.tsv")
            )
            job.name = "parse_bins_" + binner
            job.subname = "parse_bins"
            jobs.append(job)
            
            # checkm on parsed bins
            job = shotgun_metagenomics.checkm(
                os.path.join("binning", binner, "parsed_bins.tsv"),
                os.path.join("binning", binner, "parsed_bins"), 
                os.path.join("binning", binner, "parsed_bins", "out_checkm"),
                os.path.join("binning", binner, "parsed_bins", "out_checkm.txt")
            )
            job.name = "checkm_parsed_bins_" + binner
            job.subname = "checkm"
            jobs.append(job)

            job = shotgun_metagenomics.summarize_bins(
                os.path.join("binning", binner, "out_checkm.txt"),
                os.path.join("binning", binner, "raw_bins.tsv"),
                os.path.join("binning", "abundance.txt"),
                os.path.join("binning", binner),
                os.path.join("binning", binner, "summarized_bins.tsv"),
                os.path.join("binning", binner, "bins_list.txt")
            )
            job.name = "summarize_raw_bins_" + binner
            job.subname = "summarize_bins"
            jobs.append(job)
            
            job = shotgun_metagenomics.generate_link_file(
                os.path.join("binning", binner, "raw_bins.tsv"),
                os.path.join("annotations", "annotations.tsv"),
                os.path.join("binning", binner, "bins_list.txt"),
                os.path.join("binning", binner, "link.tsv")
            )
            job.name = "generate_link_raw_bins_" + binner
            job.subname = "generate_link"
            jobs.append(job)
            
            job = shotgun_metagenomics.bins_feature_table(
                os.path.join("binning", binner, "summarized_bins.tsv"),
                os.path.join("binning", binner, "raw_bins.tsv"),
                os.path.join("contig_abundance", "merged_contig_abundance_cpm.tsv"),
                os.path.join("annotations", binner, "bins", "feature_table_normalized.tsv")
            )
            job.name = "raw_bins_feature_table_normalized_" + binner
            job.subname = "bins_feature_table"
            jobs.append(job)
            
            job = shotgun_metagenomics.bins_feature_table(
                os.path.join("binning", binner, "summarized_bins.tsv"),
                os.path.join("binning", binner, "raw_bins.tsv"),
                os.path.join("contig_abundance", "merged_contig_abundance.tsv"),
                os.path.join("annotations", binner, "bins", "feature_table.tsv")
            )
            job.name = "raw_bins_feature_table_raw_counts_" + binner
            job.subname = "bins_feature_table"
            jobs.append(job)
            
            ## parsed bins
            job = shotgun_metagenomics.summarize_bins(
                os.path.join("binning", binner, "parsed_bins", "out_checkm.txt"),
                os.path.join("binning", binner, "parsed_bins.tsv"),
                os.path.join("binning", "abundance.txt"),
                os.path.join("binning", binner, "parsed_bins"),
                os.path.join("binning", binner, "parsed_bins", "summarized_bins.tsv"),
                os.path.join("binning", binner, "parsed_bins", "bins_list.txt")
            )
            job.name = "summarize_bins_" + binner
            job.subname = "summarize_bins"
            jobs.append(job)
            
            job = shotgun_metagenomics.generate_link_file(
                os.path.join("binning", binner, "parsed_bins.tsv"),
                os.path.join("annotations", "annotations.tsv"),
                os.path.join("binning", binner, "parsed_bins", "bins_list.txt"),
                os.path.join("binning", binner, "parsed_bins", "link.tsv")
            )
            job.name = "generate_link_" + binner
            job.subname = "generate_link"
            jobs.append(job)

            # create abundance table of bins (otu table format). 
            # First, using normalized abundance
            job = shotgun_metagenomics.bins_feature_table(
                os.path.join("binning", binner, "parsed_bins", "summarized_bins.tsv"),
                os.path.join("binning", binner, "parsed_bins.tsv"),
                os.path.join("contig_abundance", "merged_contig_abundance_cpm.tsv"),
                os.path.join("annotations", binner, "parsed_bins", "feature_table_normalized.tsv")
            )
            job.name = "bins_feature_table_normalized_" + binner
            job.subname = "bins_feature_table"
            jobs.append(job)
            
            # Then raw counts.
            job = shotgun_metagenomics.bins_feature_table(
                os.path.join("binning", binner, "parsed_bins", "summarized_bins.tsv"),
                os.path.join("binning", binner, "parsed_bins.tsv"),
                os.path.join("contig_abundance", "merged_contig_abundance.tsv"),
                os.path.join("annotations", binner, "parsed_bins", "feature_table.tsv")
            )
            job.name = "bins_feature_table_raw_counts_" + binner
            job.subname = "bins_feature_table"
            jobs.append(job)

            bins_feature_tables = [
                os.path.join("annotations", binner, "bins", "feature_table"),
                os.path.join("annotations", binner, "parsed_bins", "feature_table"),
                os.path.join("annotations", binner, "bins", "feature_table_normalized"),
                os.path.join("annotations", binner, "parsed_bins", "feature_table_normalized")
            ]
            
            # Alpha and Beta diversity.
            #def mkdir_p(path):
            #    try:
            #        os.makedirs(path)
            #    except OSError as exc: # Python >2.5
            #        if exc.errno == errno.EEXIST and os.path.isdir(path):
            #            pass
            #        else: raise
            mkdir_p(os.path.join("annotations", binner, "bins", "beta_div"))
            fname = os.path.join("annotations", binner, "bins", "beta_div", "tree.fasttree")
            if not os.path.isfile(fname):
                open(fname, 'a').close()
                os.utime(fname, None)
                sys.stderr.write("tree file not does exists...\n")
           
            # define files.
            # unnormalized for alpha div.
            feature_table_bins_tsv         = os.path.join("annotations", binner, "bins", "feature_table.tsv")
            # normalized for beta_div
            feature_table_bins_tsv_norm    = os.path.join("annotations", binner, "bins", "feature_table_normalized.tsv")  # same thing for bins...
            
            types_tsv  = [feature_table_bins_tsv]
            types_tsv_filtered  = [feature_table_bins_tsv_norm]
            prefixes = ["bray_curtis_feature_table_normalized"]
            prefixes_alpha_div = ["bins"]

            # Alpha div 
            mkdir_p(os.path.join("annotations", binner, "bins", "alpha_div"))
            mkdir_p(os.path.join("annotations", binner, "bins", "beta_div"))
            
            job = microbial_ecology.rtk(
                os.path.join("annotations", binner, "bins", "feature_table.tsv"),
                "./" + os.path.join("annotations", binner, "bins", "alpha_div"),
                "#FEATURE_ID",
                remove_last_col = True
            )
            job.name = "alpha_diversity_rtk_" + binner
            job.subname = "alpha_diversity"
            jobs.append(job)

            # Beta diversity.
            job = microbial_ecology.beta_diversity(
                os.path.join("annotations", binner, "bins", "feature_table_normalized.tsv"),
                "bray-curtis",
                os.path.join("annotations", binner, "bins", "beta_div"), 
                os.path.join("annotations", binner, "bins", "beta_div", "tree.fasttree"),
                "bray_curtis_feature_table_normalized"
            )
            job.name = "beta_diversity_" + "bray_curtis_feature_table_normalized_" + binner
            job.subname = "beta_diversity"
            jobs.append(job)
            
            job = microbial_ecology.principal_coordinates(
                os.path.join("annotations", binner, "bins", "beta_div", "bray_curtis_feature_table_normalized.tsv"), 
                os.path.join("annotations", binner, "bins", "beta_div", "bray_curtis_feature_table_normalized_coords.tsv") 
            )
            job.name = "beta_diversity_bray_curtis_coords_" + "feature_table_normalized_" + binner
            job.subname = "pcoa"
            jobs.append(job)
        
            job = microbial_ecology.pca_3d_plot(
                os.path.join("annotations", binner, "bins", "beta_div", "bray_curtis_feature_table_normalized_coords.tsv"), 
                os.path.join("annotations", binner, "bins", "beta_div", "3d_bray_curtis_plot") 
            )
            job.name = "beta_diversity_bray_curtis_bins_3d_plot_" + binner
            job.subname = "pca_plot"
            jobs.append(job)
        
        return jobs
    
    def anvio(self):
        """
        step anvio(): Generate files to use in the anvio software.
        """
        
        mkdir_p(os.path.join("anvio"))
        
        jobs = []
        
        job = shotgun_metagenomics.anvio_simplify_fasta_headers(
            os.path.join("assembly", "Contigs.fasta"),
            os.path.join("anvio", "Contigs_simplified.fasta")
        )
        job.name = "anvio_simplify_headers"
        job.subname = "anvio_simplify_headers"
        jobs.append(job)
        
        job = shotgun_metagenomics.anvio_gene_calling(
            os.path.join("gene_prediction", "Contigs_renamed.faa"),
            os.path.join("gene_prediction", "Contigs_renamed.gff"),
            os.path.join("annotations", "rrna", "barrnap", "bac.tsv"),
            os.path.join("anvio", "Contigs_renamed_for_anvio.tsv")
        )
        job.name = "anvio_gene_calling"
        job.subname = "anvio_gene_calling"
        jobs.append(job)
       
        job = shotgun_metagenomics.anvio_make_contigs_db(
            os.path.join("anvio", "Contigs_simplified.fasta"),
            os.path.basename(os.getcwd()),
            os.path.join("anvio", "Contigs_renamed_for_anvio.tsv"),
            os.path.join("anvio", "contigs.db"),
            os.path.join("anvio", "make_contigs_db.done")
        )
        job.name = "anvio_make_contigs_db"
        job.subname = "anvio_make_contigs_db"
        jobs.append(job)
        
        KO_assignment_file = ""
        if config.param('DEFAULT', 'KO_method', 'string') == "hmmsearch_kofam":
            KO_assignment_file = os.path.join("annotations", "hmmsearch_kofam_tblout.tsv")
        elif config.param('DEFAULT', 'KO_method', 'string') == "kofamscan":
            KO_assignment_file = os.path.join("annotations", "kofamscan.tsv")
        elif config.param('DEFAULT', 'KO_method', 'string') == "diamond_blastp":
            KO_assignment_file = os.path.join("annotations", "blastp_kegg.tsv")
        else:
            raise Exception("Error: In the .ini file, please chose between one of the following two values under the [DEFAULT] section: KO_method=kofamscan OR KO_method=diamond_blastp OR KO_method=hmmsearch_kofam")

        job = shotgun_metagenomics.anvio_generate_kegg_file(
            KO_assignment_file,
            os.path.join("annotations", "KOs_parsed.tsv"),
            os.path.join("anvio", "kegg_annotations.tsv")
        )
        job.name = "anvio_kegg"
        job.subname = "anvio_kegg"
        jobs.append(job)
        
        job = shotgun_metagenomics.anvio_generate_taxonomy_file(
            os.path.join("anvio", "Contigs_renamed_for_anvio.tsv"),
            os.path.join("annotations", "taxonomy_consensus", "taxonomy.tsv"), 
            os.path.join("anvio", "taxonomy_annotations.tsv")
        )
        job.name = "anvio_taxonomy"
        job.subname = "anvio_taxonomy"
        jobs.append(job)
        
        job = shotgun_metagenomics.anvio_generate_cog_file(
            os.path.join("annotations", "rpsblast_cog.tsv"), 
            os.path.join("anvio", "cog_annotations.tsv")
        )
        job.name = "anvio_cog"
        job.subname = "anvio_cog"
        jobs.append(job)
        
        job = shotgun_metagenomics.anvio_import_annotations(
            os.path.join("anvio", "make_contigs_db.done"),
            os.path.join("anvio", "cog_annotations.tsv"),
            os.path.join("anvio", "taxonomy_annotations.tsv"),
            os.path.join("anvio", "kegg_annotations.tsv"),
            os.path.join("anvio", "contigs.db"),
            os.path.join("anvio", "import.done")
        )
        job.name = "anvio_import"
        job.subname = "anvio_import"
        jobs.append(job)

        profiles = []
        for readset in self.readsets:
            bam_contigs = os.path.join("contig_abundance", readset.sample.name, readset.name + ".bam")
            job = shotgun_metagenomics.anvio_import_profiles(
                os.path.join("anvio", "import.done"),
                bam_contigs,
                os.path.join("anvio", "contigs.db"),
                os.path.join("anvio", "profiles", readset.name),
                os.path.join("anvio", "profiles", readset.name, "PROFILE.db")
            )
            job.name = "anvio_profile_" + readset.name
            job.subname = "anvio_profile"
            jobs.append(job)

            profiles.append(os.path.join("anvio", "profiles", readset.name, "PROFILE.db"))
        
        job = shotgun_metagenomics.anvio_merge(
            profiles,
            os.path.join("anvio", "contigs.db"),
            os.path.join("anvio", "profiles", "SAMPLES-MERGED"),
            os.path.join("anvio", "merge_profiles.done")
        )
        job.name = "anvio_merge"
        job.subname = "anvio_profile"
        jobs.append(job)
        
        job = shotgun_metagenomics.anvio_run_hmm(
            os.path.join("anvio", "merge_profiles.done"),
            os.path.join("anvio", "profiles", "SAMPLES-MERGED", "PROFILE.db"),
            os.path.join("anvio", "contigs.db"),
            os.path.join("anvio", "run_hmm.done")
        )
        job.name = "anvio_run_hmm"
        job.subname = "anvio_run_hmm"
        jobs.append(job)
        
        job = shotgun_metagenomics.anvio_import_metadata(
            os.path.join("anvio", "run_hmm.done"),
            os.path.join("anvio", "profiles", "SAMPLES-MERGED", "PROFILE.db"),
            os.path.join("anvio", "import_metadata.done")
        )
        job.name = "anvio_import_metadata"
        job.subname = "anvio_import_metadata"
        jobs.append(job)
        
        binners = config.param('DEFAULT', 'binner', 1, 'string').split(":")
        for binner in binners:
            if binner != "metabat2" and binner != "maxbin2":
                raise Exception("Error: mapper should be set to either metabat2, maxbin2 or metabat2:maxbin2 ")
        
            job = shotgun_metagenomics.anvio_generate_mag_file(
                os.path.join("binning", binner, binner + ".done"),
                os.path.join("binning", binner),
                os.path.join("anvio", binner + "_contigs.tsv")
            )
            job.name = "anvio_generate_mag_file_" + binner
            job.subname = "anvio_metabat2"
            jobs.append(job)
            
            job = shotgun_metagenomics.anvio_import_mags(
                os.path.join("anvio", "import_metadata.done"),
                os.path.join("anvio", "contigs.db"),
                os.path.join("anvio", "profiles", "SAMPLES-MERGED", "PROFILE.db"),
                "MAGs_" + binner,
                os.path.join("anvio", binner + "_contigs.tsv"),
                os.path.join("anvio", "import_mags_" + binner + ".done")
            )
            job.name = "anvio_import_mags_" + binner
            job.subname = "anvio_import_mags"
            jobs.append(job)
            
        return jobs

    def tsne(self):
        jobs = []
        
        # All contigs... TODO select by size
        job = shotgun_metagenomics.calc_kmer_freq(
            os.path.join("assembly", "Contigs.fasta"),
            os.path.join("assembly", "Contigs_tf5mer.tsv")
        
        )
        job.name = "calc_kmer_freq"
        job.subname = "calc_kmer_freq"
        jobs.append(job)
        
        job = shotgun_metagenomics.clr_transform(
            os.path.join("assembly", "Contigs_tf5mer.tsv"),
            os.path.join("assembly", "Contigs_tf5mer_clr.tsv")
        )
        job.name = "clr_transform"
        job.subname = "clr_transform"
        jobs.append(job)
        
        job = shotgun_metagenomics.bhtsne(
            os.path.join("assembly", "Contigs_tf5mer_clr.tsv"),
            os.path.join("assembly", "Contigs_tf5mer_clr_bhtsne.tsv")
        )
        job.name = "bhtsne"
        job.subname = "bhtsne"
        jobs.append(job)
        
        # Maxbin2 and metabat2
        binners = config.param('DEFAULT', 'binner', 1, 'string').split(":")

        for binner in binners:
            job = shotgun_metagenomics.cat_bins_into_single_file(
                os.path.join("binning", binner, binner + ".done"),
                os.path.join("binning", binner),
                os.path.join("binning", binner, "all.fa")
            )
            job.name = "cat_bins_" + binner
            job.subname = "cat_bins"
            jobs.append(job)
            
            job = shotgun_metagenomics.calc_kmer_freq(
                os.path.join("binning", binner, "all.fa"),
                os.path.join("binning", binner, "all_tf5mer.tsv")
            )
            job.name = "calc_kmer_freq_bins_" + binner
            job.subname = "calc_kmer_freq"
            jobs.append(job)
            
            job = shotgun_metagenomics.clr_transform(
                os.path.join("binning", binner, "all_tf5mer.tsv"),
                os.path.join("binning", binner, "all_tf5mer_clr.tsv")
            )
            job.name = "clr_transform_bins_" + binner
            job.subname = "clr_transform"
            jobs.append(job)
            
            job = shotgun_metagenomics.bhtsne(
                os.path.join("binning", binner, "all_tf5mer_clr.tsv"),
                os.path.join("binning", binner, "all_tf5mer_clr_bhtsne.tsv")
            )
            job.name = "bhtsne_bins_" + binner
            job.subname = "bhtsne"
            jobs.append(job)

            return jobs
    
    def virome(self):
        """
        Step virome(): HMMSCAN of predicted genes vs core viral protein models.
        """
        
        mkdir_p(os.path.join("annotations", "virome"))
        
        jobs = []
        
        chunks_dir = os.path.join("gene_prediction", "fasta_chunks")
        hmmsearch_out_dir = os.path.join("annotations", "virome")
        number_chunks_file = os.path.join("gene_prediction", "estimated_number_of_chunks_genes.txt")
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
            tblouts.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.tblout".format(i)))
            domtblouts.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.domtblout".format(i)))
            pfamtblouts.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.pfamtblout".format(i)))
            dones.append(os.path.join(hmmsearch_out_dir, "hmmsearch_chunk_{:07d}.done".format(i)))

        job = shotgun_metagenomics.hmmsearch_array_job(
            os.path.join(chunks_dir),
            "Contigs_renamed.faa_chunk_",
            os.path.join(hmmsearch_out_dir),
            "hmmsearch_chunk_",
            infiles,
            tblouts, domtblouts, pfamtblouts,
            dones,
            config.param('virome', 'db', required=True),
            self._curr_scheduler
        )
        job.name = "hmmsearch_pfam_virome"
        job.subname = "hmmscan"
        job.job_array_num_task = num_chunks
        jobs.append(job)

        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks_hmms(
            hmmsearch_out_dir,
            os.path.join("annotations"),
            num_chunks,
            "hmmscan",
            "virome/virome"
        )
        job.name = "hmmsearch_virome_merge"
        job.subname = "merge"      
        jobs.append(job)
        
        job = shotgun_metagenomics.parse_viral_contigs(
            os.path.join("annotations", "virome", "virome_tblout.tsv"),
            os.path.join("assembly", "Contigs.fasta"),
            os.path.join("assembly", "Contigs_length_gc.tsv"),
            os.path.join("annotations", "annotations.tsv"),
            os.path.join("annotations", "virome", "virome_master.tsv"),
            os.path.join("annotations", "virome", "viral_contigs_list.txt"),
            os.path.join("annotations", "virome", "viral_contigs.fasta")
        )
        job.name = "parse_viral_contigs"
        job.subname = "parse_viral_contigs"      
        jobs.append(job)
        
        mkdir_p(os.path.join("annotations", "virome", "chunks"))

        job = shotgun_metagenomics.exonerate_direct(
            os.path.join("annotations", "virome", "viral_contigs.fasta"),
            os.path.join("annotations", "virome", "chunks"),
            config.param("virome", "num_chunks", 1, "posint"),
            "viral_contigs.fasta"
        )
        job.name = "exonerate_viral_contigs_fna"
        job.subname = "exonerate"
        jobs.append(job)

        # Do blastn on nt of viral contigs
        mkdir_p(os.path.join("annotations", "virome", "blastn_nt_contigs"))
        
        chunks_dir = os.path.join("annotations", "virome", "chunks")
        blast_dir = os.path.join("annotations", "virome", "blastn_nt_contigs")
        #number_chunks_file = os.path.join("assembly", "estimated_number_of_chunks_contigs.txt")
        infiles = []
        outfiles = []
        dones = []
        
        num_chunks = config.param("virome", "num_chunks", 1, "posint")
        for i in range(num_chunks):
            infiles.append(os.path.join(chunks_dir, "viral_contigs.fasta_chunk_{:07d}".format(i)))
            outfiles.append(os.path.join(blast_dir, "blastn_chunk_{:07d}.tsv".format(i)))
            dones.append(os.path.join(blast_dir, "blastn_chunk_{:07d}.done".format(i)))
     
        job = shotgun_metagenomics.blastn_array_job(
            os.path.join(chunks_dir),
            "viral_contigs.fasta_chunk_",
            os.path.join(blast_dir),
            "blastn_chunk_",
            "blastn",
            infiles,
            outfiles,
            dones,
            self._curr_scheduler

        )
        job.name = "blastn_nt_viral_contigs"
        job.subname = "blastn"
        job.job_array_num_task = num_chunks
        jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join("annotations", "virome", "blastn_nt_contigs.tsv"),
            num_chunks,
            "blastn" 
        )
        job.name = "merge_blastn_nt_viral_contigs"
        job.subname = "merge_blastn"
        jobs.append(job)
        
        return jobs
        

        return jobs


    def overrepresentation(self):
        jobs = []
        
        #types = ["pathways", "modules", "KO"]
        types = ["KO"]
        for type in types:
            job = shotgun_metagenomics.kegg_overrep(
                os.path.join("gene_abundance", "merged_gene_abundance_cpm.tsv"),
                os.path.join("annotations", "KOs_parsed.tsv"),
                type,
                os.path.join("annotations", type + "_matrix_cpm.tsv")
            )
            job.name = "kegg_overrep_cpm_" + type
            job.subname = "kegg_overrep"
            jobs.append(job)
            
            job = shotgun_metagenomics.kegg_overrep(
                os.path.join("gene_abundance", "merged_gene_abundance.tsv"),
                os.path.join("annotations", "KOs_parsed.tsv"),
                type,
                os.path.join("annotations", type + "_matrix.tsv")
            )
            job.name = "kegg_overrep_read_counts_" + type
            job.subname = "kegg_overrep"
            jobs.append(job)
            
        # COG and KOG
        job = shotgun_metagenomics.cog_overrep(
            os.path.join("gene_abundance", "merged_gene_abundance_cpm.tsv"),
            os.path.join("annotations", "rpsblast_cog.tsv"),
            os.path.join("annotations", "cog_matrix_cpm.tsv")
        )
        job.name = "cog_overrep_cpm"
        job.subname = "cog_overrep"
        jobs.append(job)
        
        job = shotgun_metagenomics.cog_overrep(
            os.path.join("gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join("annotations", "rpsblast_cog.tsv"),
            os.path.join("annotations", "cog_matrix.tsv")
        )
        job.name = "COG_overrep_read_counts"
        job.subname = "cog_overrep"
        jobs.append(job)
        
        job = shotgun_metagenomics.pfam_overrep(
            os.path.join("gene_abundance", "merged_gene_abundance_cpm.tsv"),
            os.path.join("annotations", "hmmsearch_pfam_tblout_parsed.tsv"),
            os.path.join("annotations", "pfam_matrix_cpm.tsv")
        )
        job.name = "pfam_overrep_cpm"
        job.subname = "cog_overrep"
        jobs.append(job)
        
        job = shotgun_metagenomics.pfam_overrep(
            os.path.join("gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join("annotations", "hmmsearch_pfam_tblout_parsed.tsv"),
            os.path.join("annotations", "pfam_matrix.tsv")
        )
        job.name = "pfam_overrep_read_counts"
        job.subname = "cog_overrep"
        jobs.append(job)
        
        return jobs
    
    def cleanup(self):
        #Here, compress all .fastq files into .fastq.gz.
        jobs = []
        #job = shotgun_metagenomics.mymethod(
        #)
        #job.name = "myjobname"
        #jobs.append(job)
        sys.stderr.write('[DEBUG] not implemented yet\n')
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
            # Core steps.
            self.trim,
            self.remove_contam,
            self.subsample,
            self.assembly,
            self.gene_prediction,
            self.abundance,
            self.exonerate,
            self.kegg_annotation,
            self.rpsblast_cog,
            self.rpsblast_kog,
            self.hmmsearch_pfam,
            self.hmmsearch_cazy,
            self.diamond_blastp_nr,
            self.ncrna,
            self.taxonomic_annotation,
            self.beta_diversity,
            self.alpha_diversity,
            self.overrepresentation,
            self.finalize,
            self.binning,
            # Steps in development. Not officially supported.
            self.anvio,
            self.tsne,
            self.virome,
            self.differential_abundance,
            self.blastn_nt_contigs,
            self.blastn_ncbi_genomes,
            self.hmmsearch_rrna,
            self.reads_centric_taxonomy,
            self.cleanup
        ]

    def set_local_variables(self):
        self._parser_local = self.argparser

        # barcodes
        self._args_local = self._parser_local.parse_args()
        config.parse_files(self._args_local.config)
        self._config = self._args_local.config[0]
        self._config = os.path.abspath(self._config.name)
        self._curr_scheduler = self._args_local.job_scheduler
        #sys.stderr.write("[DEBUG] current job scheduler: " + str(self._curr_scheduler) + "\n")
        #self._extended_taxonomy = self._args_local.extended_taxonomy
        
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
        mkdir_p(os.path.join(root_dir, "annotations"))
        mkdir_p(os.path.join(root_dir, "annotations", "taxonomy_consensus"))
        mkdir_p(os.path.join(root_dir, "annotations", "rrna", "rdp"))
        mkdir_p(os.path.join(root_dir, "annotations", "rrna", "barrnap"))
        mkdir_p(os.path.join(root_dir, "annotations", "rrna", "abundance"))
        mkdir_p(os.path.join(root_dir, "annotations", "trna"))

    
    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_readset_file(self.args.readsets.name)
        return self._readsets
        


    def __init__(self):
        version = open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "..", "VERSION"), 'r').read().split('\n')[0]
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
               Version: """ + version + """

###############################################################################"""
        sys.stderr.write(metagenomics_string + '\n')
        time.sleep(1)
        # Add pipeline specific arguments
        self.argparser.add_argument("-r", "--readsets", help="readset file", type=argparse.FileType('r'), required=False)
        self.set_local_variables()
        sys.stderr.write('Shotgun metagenomics pipeline\n')
        super(Metagenomics, self).__init__()
                
Metagenomics().submit_jobs()
