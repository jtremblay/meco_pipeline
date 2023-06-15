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

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
#from bio.design import *
from bio.readset import *

from bio import shotgun_metagenomics
from bio import microbial_ecology

from pipelines.illumina import illumina

# Global scope variables
log = logging.getLogger(__name__)

class Metapangenomics(illumina.Illumina):

    """
    Start by processing raw reads for QC (trim + duk) and align them on the metagenome
    """

    def trim(self):
        jobs = []
        # Merge all demultiplexed fastq files in one file. One file for reads1 and one file for reads2 of Illumina paired-end.
        # If library is Illumina single end or 454 or PacBio, only one file will be generated.
        outdir = self._root_dir
        trim_stats = []
        #sys.stderr.write('outdir: ' + outdir + '\n')
        
        for readset in self.readsets:
            trim_file_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".trim.")
            trim_file_prefix_long = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".trim.long.")
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
            
    def duk(self):
        jobs=[]
        outdir = self._root_dir
        logs = []
        readset_ids = []

        for readset in self.readsets:
            trim_file_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".trim.")
            trim_file_prefix_long = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".trim.long.")
            outfile_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".")
            outfile_prefix_long = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".long.")
            out1up = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_unpaired_R1.fastq.gz")
            out2up = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_unpaired_R2.fastq.gz")
            out1p = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
            out2p = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
            
            if readset.run_type == "PAIRED_END":
                log = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".duk_contam_interleaved_log.txt")
                readset_id = readset.name

                job = shotgun_metagenomics.duk_gz(
                    trim_file_prefix + "interleaved.fastq.gz",
                    outfile_prefix + "contam.fastq",
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
                
                logs.append(log)
                readset_ids.append(readset.name)
            
            elif readset.run_type == "SINGLE_END":
                log = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".duk_contam_pair1_log.txt")

                job = shotgun_metagenomics.duk_gz(
                    trim_file_prefix + "pair1.fastq.gz",
                    outfile_prefix + "contam.fastq",
                    outfile_prefix + "ncontam.fastq.gz",
                    log,
                    config.param('DB', 'contaminants', 1, 'filepath')
                )
                job.name = "duk_single_end_reads_" + readset.sample.name
                job.subname = "duk"
                jobs.append(job)
                
                logs.append(log)
                readset_ids.append(readset.name)
            
            elif readset.run_type == "SPE_LSE":

                log = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".SPE.duk_contam_interleaved_log.txt")
                readset_id = readset.name

                job = shotgun_metagenomics.duk_gz(
                    trim_file_prefix + "interleaved.fastq.gz",
                    outfile_prefix + "contam.fastq",
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
                
                logs.append(log)
                readset_ids.append(readset.name)
    
                ## Then do long SE reads.
                log = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".long_duk_contam_pair1_log.txt")

                job = shotgun_metagenomics.duk_gz(
                    trim_file_prefix_long + "pair1.fastq.gz",
                    outfile_prefix_long + "contam.fastq",
                    outfile_prefix_long + "ncontam.fastq.gz",
                    log,
                    config.param('DB', 'contaminants', 1, 'filepath')
                )
                job.name = "duk_single_end_long_reads_" + readset.sample.name
                job.subname = "duk"
                jobs.append(job)
                
                logs.append(log)
                #readset_ids.append(readset.name)
                sys.stderr.write("[DEBUG] " + outfile_prefix_long + "ncontam.fastq.gz\n")
        # Compile duk logs.
        job = shotgun_metagenomics.merge_duk_logs_interleaved(
            logs,
            readset_ids,
            os.path.join(self._root_dir, "qced_reads", "duk_merged_logs.tsv")
        )
        job.name = "merge_duk_logs"
        job.subname = "merge_duk_logs"
        jobs.append(job)
            
        return jobs
     

    def gene_prediction(self):
        root = self._root_dir
        jobs = []
        
        infile = os.path.join(root, "assembly", "Contigs.fasta")
        outdir = os.path.join(root, "gene_prediction")
        
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
        jobs = []
        cov_list_contigs = []
        cov_list_genes = []
        flagstats_contigs_list = []
        trimmomatic_list = []
        
        # Will make index for bwa. and also bed file for computing reads spanning later.
        reference_contigs = os.path.join(self._root_dir, "assembly", "Contigs.fasta")
        bed_contigs = os.path.join(self._root_dir, "assembly", "Contigs.bed")
        bwt_contigs = os.path.join(self._root_dir, "assembly", "Contigs.fasta.bwt")
        bed_genes = os.path.join(self._root_dir, "gene_prediction", "Contigs_genes.bed")

        job = shotgun_metagenomics.make_index(
            reference_contigs,
            bwt_contigs
        )
        job.name = "make_index_contigs"
        job.subname = "make_index"
        jobs.append(job)
        
        job = shotgun_metagenomics.fasta_to_bed(
            reference_contigs,
            bed_contigs
        )
        job.name = "fasta_to_bed"
        job.subname = "fasta_to_bed"
        jobs.append(job)
        
        job = shotgun_metagenomics.gff_to_bed(
            os.path.join(self._root_dir, "gene_prediction", "Contigs_renamed.gff"),
            bed_genes
        )
        job.name = "gff_to_bed"
        job.subname = "gff_to_bed"
        jobs.append(job)
        

        for readset in self.readsets:
            bam_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".bam")
            flagstats_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".flagstats")
            flagstats_contigs_list.append(flagstats_contigs)
            trimmomatic = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".trim.stats.csv")
            trimmomatic_list.append(trimmomatic)
            cov_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".cov")
            cov_list_contigs.append(cov_contigs)
            cov_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".cov")
            cov_list_genes.append(cov_genes)
    
            outdir = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name)
            outdir_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name)
            if not os.path.exists(outdir):
                os.makedirs(os.path.join(outdir))
            if not os.path.exists(outdir_genes):
                os.makedirs(os.path.join(outdir_genes))
       
            if( readset.run_type == "PAIRED_END"):
                infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
                out1 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
                out2 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
                flag = "0x2"

                ## map against contigs
                job = shotgun_metagenomics.bwa_mem_samtools(
                    reference_contigs,
                    bwt_contigs,
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
                
            elif(readset.run_type == "SINGLE_END"):

                infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam.fastq.gz")
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
                job.name = "flagstats" + readset.sample.name
                job.subname = "flagstats"
                jobs.append(job)

            job = shotgun_metagenomics.coverage_bed(
                bam_contigs,
                bed_contigs,
                cov_contigs,
                flag
            )
            job.name = "bedtoolsCov-contigs-" + readset.sample.name
            job.subname = "bedtools"
            jobs.append(job)
            
            job = shotgun_metagenomics.coverage_bed(
                bam_contigs,
                bed_genes,
                cov_genes,
                flag
            )
            job.name = "bedtoolsCov-genes-" + readset.sample.name
            job.subname = "bedtools"
            jobs.append(job)
            
        # Once all coverage has been computed, merge all tables.
        # sys.stderr.write('gene abundance: ' + ','.join([str(x) for x in cov_list] ) + '\n')
        job = shotgun_metagenomics.merge_counts(
            cov_list_contigs,
            os.path.join(self._root_dir, "contigs_abundance", "merged_contigs_abundance.tsv"),
            #os.path.join(self._root_dir, "contigs_abundance", "merged_contigs_abundance_RPKM.tsv"),
            os.path.join(self._root_dir, "contigs_abundance", "merged_contigs_abundance_cpm.tsv"),
            "contigs"
        )
        job.name = "merge_gene_abundance_contigs"
        job.subname = "merge_gene_abundance"
        jobs.append(job)
        
        job = shotgun_metagenomics.merge_counts(
            cov_list_genes,
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            #os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance_RPKM.tsv"),
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance_cpm.tsv"),
            "genes"
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

        #Generate files for potential gam-ngs merge
        job = shotgun_metagenomics.get_insert_size(
            os.path.join(self._root_dir, "contigs_abundance"),
            os.path.join(self._root_dir, "contigs_abundance", "lib_stats.tsv"),
            os.path.join(self._root_dir, "contigs_abundance", "qc_mapping_stats.tsv")
        )
        job.name = "get_insert_size"
        job.subname = "get_insert_size"
        jobs.append(job)
         
        return jobs 
     
    def exonerate(self):
        jobs = []

        # Split for big assembly
        infile_contigs = os.path.join(self._root_dir, "assembly", "Contigs.fasta")
        chunks_dir_contigs = os.path.join(self._root_dir, "assembly", "fasta_chunks")
        number_of_chunks_file_contigs = os.path.join(self._root_dir, "assembly", "estimated_number_of_chunks_contigs.txt")
        #infile_fna = os.path.join(self._root_dir, "gene_prediction", "Contigs_renamed.fna")
        infile_faa = os.path.join(self._root_dir, "gene_prediction", "Contigs_renamed.faa")
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
        number_of_chunks_file_genes = os.path.join(self._root_dir, "gene_prediction", "estimated_number_of_chunks_genes.txt")
        
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
     
    def blastn_nt_contigs(self):
        jobs = []
        
        # Do blastn on nt for big assembly  
        chunks_dir = os.path.join(self._root_dir, "assembly", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "blastn_nt_contigs")
        #num_chunks = config.param('exonerate_contigs', 'num_fasta_chunks', type='posint')
        number_chunks_file = os.path.join(self._root_dir, "assembly", "estimated_number_of_chunks_contigs.txt")
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
            os.path.join(self._root_dir, "gene_annotation", "blastn_nt_contigs.tsv"),
            num_chunks,
            "blastn" 
        )
        job.name = "blastn_nt_big_contigs_merge"
        job.subname = "blastn"
        jobs.append(job)
        
        job = shotgun_metagenomics.keep_blast_best_hit(
            os.path.join(self._root_dir, "gene_annotation", "blastn_nt_contigs.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "blastn_nt_contigs_besthit.tsv")
        )
        job.name = "blastn_nt_big_contigs_best_hit"
        job.subname = "keep_best_hit"
        jobs.append(job)
            
        job = shotgun_metagenomics.extract_taxonomy(
            os.path.join(self._root_dir, "gene_annotation", "blastn_nt_contigs_besthit.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_contigs", "taxonomy.tsv"),
            config.param('ncbi_tax', 'accession_to_tax', required=True)
        )
        job.name = "extract_taxonomy_contigs"
        job.subname = "ncbi_tax"
        jobs.append(job)
   
        job = shotgun_metagenomics.generate_otu_table(
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_contigs", "taxonomy.tsv"),
            os.path.join(self._root_dir, "contigs_abundance", "merged_contigs_abundance.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_contigs", "otu_table.tsv")
        )
        job.name = "generate_otu_table_contigs"
        job.subname = "generate_otu_table"
        jobs.append(job)
        
        return jobs
    

    def blastn_nt(self):
        jobs = []
        
        # Do blastn on nt for co-assembly  
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "blastn_nt")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')

        for i in range(num_chunks):
            job = shotgun_metagenomics.blastn(
                os.path.join(chunks_dir, "Contigs_renamed.fana_chunk_{:07d}".format(i)),
                os.path.join(blast_dir, "blastn_chunk_{:07d}.tsv".format(i)),
                blast_dir,
                "blastn"
            )
            job.name = "blastn_nt_genes"
            job.subname = "blastn"
            jobs.append(job)
      
        
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "blastn_nt.tsv"),
            num_chunks,
            "blastn" 
        )
        job.name = "blastn_nt_genes_merge"
        job.subname = "blastn"
        jobs.append(job)
        
        job = shotgun_metagenomics.keep_blast_best_hit(
            os.path.join(self._root_dir, "gene_annotation", "blastn_nt.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "blastn_nt_besthit.tsv")
        )
        job.name = "blastn_nt_genes_best_hit"
        job.subname = "blastn_best_hit"
        jobs.append(job)
        
        job = shotgun_metagenomics.extract_taxonomy(
            os.path.join(self._root_dir, "gene_annotation", "blastn_nt_besthit.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "taxonomy.tsv")
        )
        job.name = "extract_taxonomy_genes"
        job.subname = "ncbi_tax"
        jobs.append(job)
   
        job = shotgun_metagenomics.generate_otu_table(
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "taxonomy.tsv"),
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "otu_table.tsv")
        )
        job.name = "generate_otu_table_genes"
        job.subname = "generate_otu_table"
        jobs.append(job)
        
        return jobs
    
    def diamond_blastp_kegg(self):
        jobs = []
            
        # TO FIX eventually
        fname = os.path.join(self._root_dir, "gene_annotation", "blastp_nr_annotated.tsv")
        open(fname, 'a').close()
        os.utime(fname, None)
        
        # Do blastp on KEGG for co-assembly 
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "blastp_kegg")
        #num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
        number_chunks_file = os.path.join(self._root_dir, "gene_prediction", "estimated_number_of_chunks_genes.txt")
        
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
            dones
        )
        job.name = "blastp_array_kegg"
        job.subname = "blastp_kegg"
        job.job_array_num_task = num_chunks
        jobs.append(job)
        
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "blastp_kegg.tsv"),
            num_chunks,
            "blastp"
        )
        job.name = "blastp_kegg"
        job.subname = "merge"
        jobs.append(job)
        
        # Generate a clean table of module/pathways.
        job = shotgun_metagenomics.parse_kegg(
            os.path.join(self._root_dir, "gene_annotation", "blastp_kegg.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "blastp_kegg_parsed.tsv")
        )
        job.name = "parse_kegg"
        job.subname = "merge"
        jobs.append(job)

        return jobs 
    
    
    def rpsblast_cog(self):
        jobs = []
        
        # Do rpsblast on COG for big assembly  
        if(config.param('blastp', 'skip_cog', 1, 'string') == 'yes'):
            fname = os.path.join(self._root_dir, "gene_annotation", "rpsblast_cog.tsv")
            open(fname, 'a').close()
            os.utime(fname, None)
        else:
            chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
            blast_dir = os.path.join(self._root_dir, "gene_annotation", "rpsblast_cog")
            num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
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
                config.param('cog', 'db') 
            )
            job.name = "rpsblast_cog"
            job.subname = "rpsblast"
            job.job_array_num_task = num_chunks
            jobs.append(job)
    
            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks(
                blast_dir,
                os.path.join(self._root_dir, "gene_annotation", "rpsblast_cog.tsv"),
                num_chunks,
                "rpsblast"
            )
            job.name = "rpsblast_cog"
            job.subname = "merge"
            jobs.append(job)
        
        return jobs 
    
    def rpsblast_kog(self):
        jobs = []
        
        # Do rpsblast on COG for big assembly  
        if(config.param('DEFAULT', 'skip_kog', 1, 'string') == 'yes'):
            fname = os.path.join(self._root_dir, "gene_annotation", "rpsblast_kog.tsv")
            open(fname, 'a').close()
            os.utime(fname, None)
        else:
            chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
            blast_dir = os.path.join(self._root_dir, "gene_annotation", "rpsblast_kog")
            num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
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
                config.param('kog', 'db') 
            )
            job.name = "rpsblast_kog"
            job.subname = "rpsblast"
            job.job_array_num_task = num_chunks
            jobs.append(job)

            #for i in range(num_chunks):
            #    job = shotgun_metagenomics.rpsblast(
            #        os.path.join(chunks_dir, "Contigs_renamed.faa_chunk_{:07d}".format(i)),
            #        os.path.join(blast_dir, "rpsblast_chunk_{:07d}.tsv".format(i)),
            #        blast_dir,
            #        config.param('kog', 'db') 
            #    )
            #    job.name = "rpsblast_kog"
            #    job.subname = "rpsblast"
            #    jobs.append(job)
          
            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks(
                blast_dir,
                os.path.join(self._root_dir, "gene_annotation", "rpsblast_kog.tsv"),
                num_chunks,
                "rpsblast"
            )
            job.name = "rpsblast_kog_merge"
            job.subname = "merge"
            jobs.append(job)
            
        return jobs 
    
    
    def hmmscan_pfam(self):
        jobs = []
        
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "hmmscan_pfam")
        #num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint', required=True)
        number_chunks_file = os.path.join(self._root_dir, "gene_prediction", "estimated_number_of_chunks_genes.txt")
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
            config.param('pfam', 'db', required=True) 
        )
        job.name = "hmmscan_pfam"
        job.subname = "hmmscan"
        job.job_array_num_task = num_chunks
        jobs.append(job)

        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks_hmms(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation"),
            num_chunks,
            "hmmscan",
            "hmmscan_pfam"
        )
        job.name = "hmmscan_pfam_merge"
        job.subname = "merge"      
        jobs.append(job)
        
        return jobs

    def diamond_blastp_nr(self):
        jobs = []
       
        # For now, create/touch dummy file. To improve/fix eventually.
        #if(config.param('blastp', 'skip_blastp_nr', 1, 'string') == 'yes'):
        #sys.stderr.write("Skipping diamond_blastp_nr step...\n")
        #fname = os.path.join(self._root_dir, "gene_annotation", "blastp_nr_annotated.tsv")
        #open(fname, 'a').close()
        #os.utime(fname, None)
        
        # Do blastp on KEGG for big assembly  
        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "blastp_nr")
        #num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
        number_chunks_file = os.path.join(self._root_dir, "gene_prediction", "estimated_number_of_chunks_genes.txt")
        
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
            dones
        )
        job.name = "blastp_array_nr"
        job.subname = "blastp_nr"
        job.job_array_num_task = num_chunks
        jobs.append(job)
        
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "blastp_nr.tsv"),
            num_chunks,
            "blastp"
        )
        job.name = "merge_diamond_blastp_nr"
        job.subname = "merge"
        jobs.append(job)
        
        job = shotgun_metagenomics.keep_blast_best_hit(
            os.path.join(self._root_dir, "gene_annotation", "blastp_nr.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "blastp_nr_besthit.tsv")
        )   
        job.name = "diamond_blastp_nr_genes_best_hit"
        job.subname = "blastp_best_hit"
        jobs.append(job)
        
        job = shotgun_metagenomics.extract_taxonomy(
            os.path.join(self._root_dir, "gene_annotation", "blastp_nr_besthit.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "taxonomy.tsv"),
            config.param('ncbi_tax', 'accession_to_tax_nr', required=True)
        )
        job.name = "extract_taxonomy_genes"
        job.subname = "ncbi_tax"
        jobs.append(job)
   
        job = shotgun_metagenomics.generate_otu_table(
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "taxonomy.tsv"),
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "otu_table.tsv")
        )
        job.name = "generate_otu_table_genes"
        job.subname = "generate_otu_table"
        jobs.append(job)
        
    
        return jobs 
    

    def hmmscan_rrna(self):
        jobs = []
        
        #chunks_dir = os.path.join(self._root_dir, "assembly", "fasta_chunks")
        rrna_dir = os.path.join(self._root_dir, "gene_annotation", "rrna")
        #num_chunks = config.param('exonerate_contigs', 'num_fasta_chunks', type='posint', required=True)
        job = shotgun_metagenomics.hmmscan(
            os.path.join(self._root_dir, "assembly", "Contigs.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna"),
            rrna_dir,
            config.param('rrna', 'db', required=True) 
        )
        job.name = "hmmscan_rrna"
        job.subname = "hmmscan"
        jobs.append(job)
        
        # Split rnammer fasta output between subunits.
        job = shotgun_metagenomics.split_rrna(
            os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna.domtblout"),
            os.path.join(self._root_dir, "assembly", "Contigs.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "arc_5S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "bac_5S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "euk_5S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "arc_16S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "bac_16S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "euk_18S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "arc_23S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "bac_23S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "euk_28S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rrna", "arc_bac_16S.fasta")
        )
        job.name = "split_rrna"
        job.subname = "split_rrna"
        jobs.append(job)
        
        #prefixes = ["rrna_5S", "rrna_8S",  "rrna_16S", "rrna_18S", "rrna_23S", "rrna_28S"]
        #prefixes = ["rrna_5S", "rrna_16S", "rrna_23S"]
        prefixes = ["arc_bac_16S"]
            
        # get ncbi id with blastn against nt AND perform RDP classifier.
        # Note that blastn otu tables will be generated in a later step: taxonomy_annotation.
        for prefix in prefixes:
            # create bed files for abundance - we have to only get the rea
            # Just take taxonomy file which contains coord for each rRNA genes, then generate bed file.
            job = shotgun_metagenomics.rrna_to_bed(
                os.path.join(self._root_dir, "gene_annotation", "rrna", prefix + ".fasta"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", prefix + ".bed")
            )
            job.name = "rnammer_to_bed_" + prefix
            job.subname = "rnammer_to_bed"
            jobs.append(job)

            # Then from this new bed file, process all bams to just get reads that falls into coords of new bed file
            # which corresponds to rRNA genes. Abundance will be used for both blastn and rdp taxonomy.
            cov_list = []
            for readset in self.readsets:
            
                if(readset.run_type == "SINGLE_END"):
                    flag = "0x0"
                elif(readset.run_type == "PAIRED_END"):
                    flag = "0x2"

                bam_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".bam")
            
                job = shotgun_metagenomics.coverage_bed(
                    bam_contigs,
                    os.path.join(self._root_dir, "gene_annotation", "rrna", prefix + ".bed"),
                    os.path.join(self._root_dir, "gene_annotation", "rrna", "abundance", prefix + "_" + readset.name + ".cov"),
                    flag
                )
                job.name = "bedtoolsCov-contigs-rrna_" + readset.sample.name + "_" + prefix
                job.subname = "bedtools"
                jobs.append(job)
                cov_list.append(os.path.join(self._root_dir, "gene_annotation", "rrna", "abundance", prefix + "_" + readset.name + ".cov"))
        
            job = shotgun_metagenomics.merge_counts(
                cov_list,
                os.path.join(self._root_dir, "gene_annotation", "rrna", prefix + "_merged_abundance.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", prefix + "_merged_abundance_RPKM.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", prefix + "_merged_abundance_cpm.tsv"),
                "genes"
            )
            job.name = "merge_gene_abundance_contigs_" + prefix
            job.subname = "merge_gene_abundance"
            jobs.append(job)
            
            # RDP classifier with rrna results.
            job = microbial_ecology.rdp_wrapper(
                os.path.join(self._root_dir, "gene_annotation", "rrna", prefix + ".fasta"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_rdp", prefix + "_rdp.tsv")
            )
            job.name = "classify_" + prefix
            job.subname = "RDP"
            jobs.append(job)
           
            # Convert RDP table to in-house taxonomy format.
            job = shotgun_metagenomics.rdp_to_taxonomy(
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_rdp", prefix + "_rdp.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_rdp", prefix + "_rdp_taxonomy.tsv")
            )
            job.name = "rdp_to_taxonomy_rrna_" + prefix
            job.subname = "rdp_to_taxonomy"
            jobs.append(job)

            # Then with rdp output, generate otu table
            job = shotgun_metagenomics.generate_otu_table(
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_rdp", prefix + "_rdp_taxonomy.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", prefix + "_merged_abundance.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_rdp", prefix +  "_otu_table.tsv")
            )
            job.name = "generate_otu_table"
            job.subname = "generate_otu_table"
            jobs.append(job)
            
            # Then, once abundance is done, BLASTN on nt NCBI. OTU table generation will be done later in taxonomic_annotation step!
            job = shotgun_metagenomics.blastn(
                os.path.join(self._root_dir, "gene_annotation", "rrna", prefix + ".fasta"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_blastn", "blastn_nt_" + prefix + ".tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_blastn"),
                "blastn"
            )
            job.name = "blastn_nt_" + prefix
            job.subname = "blastn"
            jobs.append(job)
            
            job = shotgun_metagenomics.extract_taxonomy(
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_blastn", "blastn_nt_" + prefix + ".tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_blastn", "blastn_nt_taxonomy_" + prefix + ".tsv")
            )
            job.name = "extract_taxonomy_rrna_blastn_" + prefix
            job.subname = "ncbi_tax"
            jobs.append(job)
   
            job = shotgun_metagenomics.generate_otu_table(
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_blastn", "blastn_nt_taxonomy_" + prefix + ".tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", prefix + "_merged_abundance.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rrna", "rrna_blastn", prefix + "_otu_table.tsv")
            )
            job.name = "generate_otu_table_rrna_blastn_" + prefix
            job.subname = "generate_otu_table"
            jobs.append(job)

        return jobs
    
    
    # Generate taxonomic summary from various otu tables generated with various methods.
    def taxonomic_annotation(self):
        jobs = []
        
        job = shotgun_metagenomics.generate_otu_table_consensus(
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "taxonomy.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_contigs", "taxonomy.tsv"),
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "blastn_nt_contigs_besthit.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "blastp_nr_besthit.tsv"),
            os.path.join(self._root_dir, "gene_prediction", "Contigs_renamed.gff"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_consensus", "otu_table.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_consensus", "taxonomy.tsv")
        )
        job.name = "generate_otu_table_genes_consensus"
        job.subname = "generate_otu_table"
        jobs.append(job)
   
        if self._extended_taxonomy:
            prefixes = [
                "genes",
                "contigs"
            ]
            #    "reads_centric",
            #    "rnammer_16S_rdp",
            #    "rnammer_16S_blastn",
            #    "emirge_rdp",
            #    "emirge_blastn"
        
            infiles = [
                os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "otu_table.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "taxonomy_contigs", "otu_table.tsv")
            ]
            #     os.path.join(self._root_dir, "gene_annotation", "taxonomy_reads_centric", "otu_table.tsv"),
            #     os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_rdp", "rnammer_16S_otu_table.tsv"),
            #     os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_blastn", "rnammer_16S_otu_table.tsv"),
            #     os.path.join(self._root_dir, "gene_annotation", "emirge", "otu_table_rdp.tsv"),
            #     os.path.join(self._root_dir, "gene_annotation", "emirge", "otu_table_blastn.tsv")
            
            directories = [
                os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes"),              #cpm with removal of low counts.
                os.path.join(self._root_dir, "gene_annotation", "taxonomy_contigs")            #cpm with removal of low counts.
            ]
            #    os.path.join(self._root_dir, "gene_annotation", "taxonomy_reads_centric"),      #perc
            #    os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_rdp"),      #perc
            #    os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_blastn"),   #perc
            #    os.path.join(self._root_dir, "gene_annotation", "emirge"),                      #perc (already normed.)
            #    os.path.join(self._root_dir, "gene_annotation", "emirge")                       #perc (already normed.)
            #]

        else:
            if os.path.isfile(os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "otu_table.tsv")):

                #sys.stderr.write("taxonomy_genes/otu_table.tsv found...\n")

                prefixes = [
                    "genes", 
                    "contigs",
                    "consensus"
                ]
                infiles = [
                    os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "otu_table.tsv"),
                    os.path.join(self._root_dir, "gene_annotation", "taxonomy_contigs", "otu_table.tsv"),
                    os.path.join(self._root_dir, "gene_annotation", "taxonomy_consensus", "otu_table.tsv")
                ]
                directories = [
                    os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes"),
                    os.path.join(self._root_dir, "gene_annotation", "taxonomy_contigs"),
                    os.path.join(self._root_dir, "gene_annotation", "taxonomy_consensus")
                ]
                abundances = [
                    os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance_cpm.tsv"),
                    os.path.join(self._root_dir, "contigs_abundance", "merged_contigs_abundance_cpm.tsv"),
                    os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance_cpm.tsv") # Because genes and not contigs are used in 
                ] 
            
            else:
                
                #sys.stderr.write("taxonomy_genes/otu_table.tsv found...\n")
                
                prefixes = ["contigs"]
                infiles = [os.path.join(self._root_dir, "gene_annotation", "taxonomy_contigs", "otu_table.tsv")]
                directories = [os.path.join(self._root_dir, "gene_annotation", "taxonomy_contigs")]
                abundances = [os.path.join(self._root_dir, "contigs_abundance", "merged_contigs_abundance_cpm.tsv")] 
                # Create dummy output files. Probably find a more elegant way of doing that.
                sys.stderr.write("Skipping gene centric taxonomic step...\n")
                fname = os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "taxonomy.tsv")
                open(fname, 'a').close()
                os.utime(fname, None)
                fname = os.path.join(self._root_dir, "gene_annotation", "blastn_nt_besthit.tsv")
                open(fname, 'a').close()
                os.utime(fname, None)
  
        # Then proces all otu tables including the extended ones if -e option
        for i in range(0, len(infiles)):
            job = microbial_ecology.convert_otu_to_biom_hdf5(
                infiles[i],
                os.path.join(directories[i], "otu_table.biom")
            )
            job.name = "convert_otu_to_biom_" + prefixes[i]
            job.subname = "convert_otu_table"
            jobs.append(job)

            # Generate tsv and biom table of raw counts (unnormalized)
            job = microbial_ecology.split_otu_table(
                os.path.join(directories[i], "otu_table.tsv"),
                os.path.join(directories[i], "otu_table_bacteriaArchaea.tsv"),
                os.path.join(directories[i], "otu_table_others.tsv"),
                "bacteriaArchaea"
            )
            job.name = "split_otu_table_bactArch_unnormalized_" + prefixes[i]
            job.name = "split_otu_table"
            jobs.append(job)
    
            job = microbial_ecology.convert_otu_to_biom_hdf5(
                os.path.join(directories[i], "otu_table_bacteriaArchaea.tsv"),
                os.path.join(directories[i], "otu_table_bacteriaArchaea.biom")
            )
            job.name = "convert_otu_table_to_biom_bactArch_unnormalized_" + prefixes[i]
            job.subname = "convert_otu_table"
            jobs.append(job)
            
            job = microbial_ecology.convert_otu_to_biom_hdf5(
                os.path.join(directories[i], "otu_table_others.tsv"),
                os.path.join(directories[i], "otu_table_others.biom")
            )
            job.name = "convert_otu_table_to_biom_others_" + prefixes[i]
            job.subname = "convert_otu_table"
            jobs.append(job)
    

            # This job doesnt really do any normalization at all. It takes cpm values and creates an otu table.
            job = microbial_ecology.normalize_cpm_otu_table_hdf5(
                abundances[i],
                os.path.join(directories[i], "otu_table.tsv"),
                os.path.join(directories[i], "otu_table_normalized.biom"),
                os.path.join(directories[i], "otu_table_normalized.tsv")
            )
            job.name = "normalize_otu_table_all_" + prefixes[i]
            job.subname = "normalization"
            jobs.append(job)
    
            job = microbial_ecology.split_otu_table(
                os.path.join(directories[i], "otu_table_normalized.tsv"),
                os.path.join(directories[i], "otu_table_normalized_bacteriaArchaea.tsv"),
                os.path.join(directories[i], "otu_table_normalized_others.tsv"),
                "bacteriaArchaea"
            )
            job.name = "split_otu_table_bactArch_" + prefixes[i]
            job.name = "split_otu_table"
            jobs.append(job)
            
            job = microbial_ecology.convert_otu_to_biom_hdf5(
                os.path.join(directories[i], "otu_table_normalized_bacteriaArchaea.tsv"),
                os.path.join(directories[i], "otu_table_normalized_bacteriaArchaea.biom")
            )
            job.name = "convert_otu_table_to_biom_bactArch_" + prefixes[i]
            job.subname = "convert_otu_table"
            jobs.append(job)
            
            job = microbial_ecology.convert_otu_to_biom_hdf5(
                os.path.join(directories[i], "otu_table_normalized_others.tsv"),
                os.path.join(directories[i], "otu_table_normalized_others.biom")
            )
            job.name = "convert_otu_table_to_biom_others_" + prefixes[i]
            job.subname = "convert_otu_table"
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
                                    os.path.join(directories[i], "otu_table" + normalizations[n] + organisms[k] + ".biom"),
                                    m,
                                    os.path.join(directories[i], organisms2[k], types[j]) 
                                )
                                job.name = "tax_summary" + normalizations[n] + "_" + types[j] + "_" + organisms2[k]  + "_L" + str(m) + "_" + prefixes[i]
                                job.subname = "summarize_taxonomy"
                                jobs.append(job)
                            elif types[j] == "relative":
                                job = microbial_ecology.summarize_taxonomy_relative(
                                    os.path.join(directories[i], "otu_table" + normalizations[n] + organisms[k] + ".biom"),
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
                                os.path.join(directories[i], organisms2[k], types[j], "otu_table" + normalizations[n] + organisms[k] + "_L" + str(m) + ".txt"),
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
            os.path.join(self._root_dir, "gene_prediction", "Contigs_renamed.gff"),
            os.path.join(self._root_dir, "assembly", "Contigs.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "blastp_kegg_parsed.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "hmmscan_pfam_tblout.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "rpsblast_cog.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "rpsblast_kog.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_consensus", "taxonomy.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "blastp_nr_annotated.tsv"),
            # outfiles
            os.path.join(self._root_dir, "gene_annotation", "Contigs.gff"),
            os.path.join(self._root_dir, "gene_annotation", "Contigs.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "annotations.tsv")
        )
        job.name = "generate_gff"
        job.subname = "generate_gff"
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
        
        return [
            # Core steps.
            self.trim,#1
            self.duk,#2
            self.gene_prediction,#4
            self.abundance,#5
            self.exonerate,#6
            self.blastn_nt_contigs,#7
            self.diamond_blastp_kegg,#8
            self.rpsblast_cog,#9
            self.rpsblast_kog,#10
            self.hmmscan_pfam,#11
            self.diamond_blastp_nr,#12
            self.taxonomic_annotation,#13
        ]

    def set_local_variables(self):
        self._parser_local = argparse.ArgumentParser(description='Process options.')
        self._parser_local.add_argument("-c", "--config", help="config INI-style file", nargs="+", type=file, required=True)
        self._parser_local.add_argument("-r", "--readsets", help="readset file", type=file, required=False)
        self._parser_local.add_argument("-s", "--steps", help="step range e.g. '1-5', '3,6,7', '2,4-8'", required=True)
        self._parser_local.add_argument("-o", "--output-dir", help="output directory (default: current)", default=os.getcwd())
        self._parser_local.add_argument("-j", "--job-scheduler", help="job scheduler type (default: slurm)", choices=["torque", "batch", "daemon", "slurm"], default="slurm")
        self._parser_local.add_argument("-f", "--force", help="force creation of jobs even if up to date (default: false)", action="store_true")
        self._parser_local.add_argument("-n", "--normalize", help="Nomalize before assembly (default: false)", action="store_true", required=False)
        self._parser_local.add_argument("-e", "--extended-taxonomy", help="Perform extended taxonomy analyses (default: false)", action="store_true", required=False)
        self._parser_local.add_argument("-l", "--log", help="log level (default: info)", choices=["debug", "info", "warning", "error", "critical"], default="info")
        self._parser_local.add_argument("-z", "--json", help="generate pipeline path in json format", default=sys.stdout, type=argparse.FileType('w'), required=False)

        # barcodes
        self._args_local = self._parser_local.parse_args()
        config.parse_files(self._args_local.config)
        self._config = self._args_local.config[0]
        self._config = os.path.abspath(self._config.name)
        self._normalize = self._args_local.normalize
        self._extended_taxonomy = self._args_local.extended_taxonomy
        
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
        mkdir_p(os.path.join(root_dir, "gene_prediction"))
        mkdir_p(os.path.join(root_dir, "gene_abundance"))
        mkdir_p(os.path.join(root_dir, "contigs_abundance"))
        mkdir_p(os.path.join(root_dir, "annotations"))

    
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
        self.set_local_variables()
        self.argparser.add_argument("-n", "--normalize-before-assembly", help="Nomalize before assembly (default: false)", action="store_true", required=False)
        self.argparser.add_argument("-e", "--extended-taxonomy", help="Perform extended taxonomy analyses (default: false)", action="store_true", required=False)
        sys.stderr.write('Metagenomics pipeline\n')
        super(Metagenomics, self).__init__()
                
Metagenomics().submit_jobs()
