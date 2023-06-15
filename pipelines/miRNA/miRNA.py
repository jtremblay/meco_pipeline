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

# Append pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# GenPipes/CAF Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.readset import *

from bio import shotgun_metagenomics
from bio import microbial_ecology
from bio import mirna

from pipelines import common
#from pipelines.illumina import illumina

# Global scope variables
log = logging.getLogger(__name__)

class MiRNA(common.CAFPipeline):
    """
    miRNA Pipeline
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
        # Merge all demultiplexed fastq files in one file. One file for reads1 and one file for reads2 of Illumina paired-end.
        # If library is Illumina single end or 454 or PacBio, only one file will be generated.
        #outdir = self._root_dir
        trim_stats = []
        #sys.stderr.write('outdir: ' + outdir + '\n')
        
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

            elif readset.run_type == "SINGLE_END" :
                if config.param("DEFAULT", "qc_methods", 1, "string") == "trimmomatic":

                    job = mirna.trimmomatic_se(
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
        #outdir = self._root_dir
        logs = []
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
     
    def filter_by_size(self):
        """
        Step filter_by_size(): filter fastqs by size
        """
        #Here, compress all .fastq files into .fastq.gz.
        jobs = []
        
        for readset in self.readsets:
            
            job = mirna.size_select(
                os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam.fastq.gz"),
                os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_mirna.fastq.gz"),
                os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_gtmirna.fastq.gz")
            )
            job.name = "size_select_fastq-" + readset.sample.name
            job.subname = "size_select"
            jobs.append(job)

        return jobs
    
    def align(self):
        """
        Step align(): Align reads against other ncRNAs.
        """
        # Here combine all fastqs and dereplicate at 100% id.
        jobs = []
        
        #cov_list = []
        flagstats_contigs_list = []
        trimmomatic_list = []

        for readset in self.readsets:
            if not os.path.exists(os.path.join("alignments", readset.sample.name)):
                os.makedirs(os.path.join("alignments", readset.sample.name))
            
            if self._type == "denovo_sub": 
                job = mirna.bwa_aln_samtools(
                    config.param('bwa', 'reference_rfam', required=True, type='filepath'),
                    config.param('bwa', 'bwt_rfam', required=True, type='filepath'),
                    os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_mirna.fastq.gz"),
                    os.path.join("alignments", readset.sample.name, readset.name + ".bam")
                )
                job.name = "bwa-aln_rfam_" + readset.sample.name
                job.subname = "bwa"
                jobs.append(job)
                
                job = mirna.extract_unmapped(
                    os.path.join("alignments", readset.sample.name, readset.name + ".bam"),
                    os.path.join("alignments", readset.sample.name, readset.name + "_unmapped.bam"),
                    os.path.join("alignments", readset.sample.name, readset.name + "_unmapped.fastq")
                )
                job.name = "extract_unmapped_rfam_" + readset.sample.name
                job.subname = "extract_unmapped"
                jobs.append(job)
                
                job = mirna.bwa_aln_samtools(
                    config.param('bwa', 'reference_5S', required=True, type='filepath'),
                    config.param('bwa', 'bwt_5S', required=True, type='filepath'),
                    os.path.join("alignments", readset.sample.name, readset.name + "_unmapped.fastq.gz"),
                    os.path.join("alignments", readset.sample.name, readset.name + "_unmapped2a.bam")
                )
                job.name = "bwa-aln_5S_" + readset.sample.name
                job.subname = "bwa"
                jobs.append(job)
                
                job = mirna.extract_unmapped(
                    os.path.join("alignments", readset.sample.name, readset.name + "_unmapped2a.bam"),
                    os.path.join("alignments", readset.sample.name, readset.name + "_unmapped2b.bam"),
                    os.path.join("alignments", readset.sample.name, readset.name + "_unmapped2b.fastq")
                )
                job.name = "extract_unmapped_5S_" + readset.sample.name
                job.subname = "extract_unmapped"
                jobs.append(job)
                
            
            elif self._type == "denovo": 
                job = mirna.bwa_aln_samtools(
                    config.param('bwa', 'reference_mirbase', required=True, type='filepath'),
                    config.param('bwa', 'bwt_mirbase', required=True, type='filepath'),
                    os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_mirna.fastq.gz"),
                    os.path.join("alignments", readset.sample.name, readset.name + ".bam")
                )
                job.name = "bwa-aln_rfam_" + readset.sample.name
                job.subname = "bwa"
                jobs.append(job)
                
                job = mirna.extract_mapped(
                    os.path.join("alignments", readset.sample.name, readset.name + ".bam"),
                    os.path.join("alignments", readset.sample.name, readset.name + "_mapped.bam"),
                    os.path.join("alignments", readset.sample.name, readset.name + "_mapped.fastq")
                )
                job.name = "extract_unmapped_rfam_" + readset.sample.name
                job.subname = "extract_unmapped"
                jobs.append(job)
            
            elif self._type == "reference": 
                #flag = "0x0"
                #cov = os.path.join("alignments", readset.sample.name, readset.name + ".cov")
                #cov_list.append(cov)
                flagstats_contigs = os.path.join("alignments", readset.sample.name, readset.name + ".bam.flagstat")
                flagstats_contigs_list.append(flagstats_contigs)
                trimmomatic = os.path.join("qced_reads", readset.sample.name, readset.name + ".trim.stats.csv")
                trimmomatic_list.append(trimmomatic)
                
                job = mirna.bwa_aln_samtools(
                    config.param('bwa', 'reference_host', required=True, type='filepath'),
                    config.param('bwa', 'bwt_host', required=True, type='filepath'),
                    os.path.join("qced_reads", readset.sample.name, readset.name + ".ncontam_mirna.fastq.gz"),
                    os.path.join("alignments", readset.sample.name, readset.name + ".bam")
                )
                job.name = "bwa-aln_host_" + readset.sample.name
                job.subname = "bwa"
                jobs.append(job)
                
                job = shotgun_metagenomics.flagstats(
                    os.path.join("alignments", readset.sample.name, readset.name + ".bam"),
                    os.path.join("alignments", readset.sample.name, readset.name + ".bam.flagstat")
                )
                job.name = "flagstats-" + readset.sample.name
                job.subname = "flagstats"
                jobs.append(job)

            job = mirna.coverage_bed(
                os.path.join("alignments", readset.sample.name, readset.name + ".bam"),
                config.param('bedtools', 'gff', required=True, type='filepath'),
                os.path.join("alignments", readset.sample.name, readset.name + ".cov")
                #flag
            )
            job.name = "bedtoolsCov-" + readset.sample.name
            job.subname = "bedtools"
            jobs.append(job)
             
            
        job = mirna.merge_flagstats(
            flagstats_contigs_list,
            trimmomatic_list,
            os.path.join("alignments", "qc_mapping_stats.tsv")
        )
        job.name = "flagstats_merge"
        job.subname = "flagstats"
        jobs.append(job)


        return jobs
    
    def abundance(self):
        """
        Step abundance(): Generate abundance matrix for miRNA.
        """
        # Here combine all fastqs and dereplicate at 100% id.
        jobs = []
        infiles = []
        cov_list = []
        if not os.path.exists(os.path.join("abundance")):
            os.makedirs(os.path.join("abundance"))
        
        if self._type == "denovo_sub": 
            for readset in self.readsets:
                #infiles.append(os.path.join("alignments", readset.sample.name, readset.name + "_unmapped.fastq.gz"))
                infiles.append(os.path.join("alignments", readset.sample.name, readset.name + "_unmapped2b.fastq.gz"))
        
        elif self._type == "denovo": 
            for readset in self.readsets:
                infiles.append(os.path.join("alignments", readset.sample.name, readset.name + "_mapped.fastq.gz"))
        
        elif self._type == "reference": 
            for readset in self.readsets:
                cov = os.path.join("alignments", readset.sample.name, readset.name + ".cov")
                cov_list.append(cov)
            
        if self._type == "denovo" or self._type == "denovo_sub":
            job = mirna.dereplicate(
                infiles,
                os.path.join("abundance", "derep.fna")
            )
            job.name = "dereplicate_miRNAs"
            job.subname = "dereplicate"
            jobs.append(job)
            
            job = mirna.generate_matrix(
                os.path.join("abundance", "derep.fna"),
                os.path.join("abundance")
            )
            job.name = "generate_abundance_matrix_miRNAs"
            job.subname = "dereplicate"
            jobs.append(job)
            
            job = mirna.filter_by_abundance(
                os.path.join("abundance", "obs.tsv"),
                os.path.join("abundance", "obs.fasta"),
                os.path.join("abundance", "obs_filtered.tsv"),
                os.path.join("abundance", "obs_filtered.fasta"),
                os.path.join("abundance", "obs_filtered_rc.fasta")
            )
            job.name = "filter_miRNA_by_abundance"
            job.name = "filter"
            jobs.append(job)

            job = mirna.normalize(
                os.path.join("abundance", "obs_filtered.tsv"),
                os.path.join("abundance", "obs_filtered.fasta"),
                os.path.join("abundance", "obs_filtered_rpkm.tsv")
            )
            job.name = "normalize_miRNA_matrix"
            job.name = "normalize"
            jobs.append(job)
        
        elif self._type == "reference": 
            # Once all coverage has been computed, merge all tables.
            # sys.stderr.write('gene abundance: ' + ','.join([str(x) for x in cov_list] ) + '\n')
            if not os.path.exists(os.path.join("abundance")):
                os.makedirs(os.path.join("abundance"))

            job = mirna.merge_counts(
                cov_list,
                os.path.join("abundance", "miRNA_abundance.tsv"),
                os.path.join("abundance", "miRNA_abundance_cpm.tsv"),
                os.path.join("abundance", "miRNA_abundance_rpkm.tsv")
            )
            job.name = "merge_gene_abundance"
            job.subname = "merge_gene_abundance"
            jobs.append(job)
           
            # Generate rc for alignment.
            job = mirna.generate_fasta_file(
                os.path.join("abundance", "miRNA_abundance.tsv"),
                os.path.join("abundance", "miRNA_abundance.fasta"),
                os.path.join("abundance", "miRNA_abundance_rc.fasta"),
                config.param("bwa", "reference_mirbase_for_query", 1, "filepath")
            )
            job.name = "generate_fasta_file"
            job.subname = "generate_fasta_file"
            jobs.append(job)

        return jobs
   
    def exonerate(self):
        
        """
        Step exonerate(): Split co-assembly and predicted genes fasta files in smaller chunks to perform downstream annotations.
        """
        jobs = []

        #infile_fna = os.path.join("gene_prediction", "Contigs_renamed.fna")
        infile = os.path.join("abundance", "obs_filtered.fasta")
        chunks_dir = os.path.join("abundance", "fasta_chunks")
        number_of_chunks_file = os.path.join("abundance", "estimated_number_of_chunks.txt")
        
        #FAA genes (i.e. predicted genes)
        job = shotgun_metagenomics.estimate_number_of_chunks(
            infile,
            number_of_chunks_file,
            str(config.param('exonerate', 'targeted_chunk_file_size', 1, 'int'))
        )
        job.name = "estimate_chunks_file_size"
        jobs.append(job)
        
        job = shotgun_metagenomics.exonerate(
            infile,
            chunks_dir,
            number_of_chunks_file,
            "obs_filtered.fasta" 
        )
        job.name = "exonerate_fna"
        job.subname = "exonerate"
        jobs.append(job)
        
        return jobs    
    
    def identify_mirbase_mature(self):
        """
        Step blastn mirna against DB(s).
        """
        if not os.path.exists(os.path.join("annotations", "blastn_mirbase_mature")):
            os.makedirs(os.path.join("annotations", "blastn_mirbase_mature"))
        
        jobs = []
        
        if self._type == "denovo" or self._type == "denovo_sub":
            
            # Do blastn on nt for big assembly  
            chunks_dir = os.path.join("abundance", "fasta_chunks")
            blast_dir = os.path.join("annotations", "blastn_mirbase_mature")
            number_chunks_file = os.path.join("abundance", "estimated_number_of_chunks.txt")
            
            infiles = []
            outfiles = []
            dones = []
            
            #num_chunks = 0
            #if os.path.exists(number_chunks_file) and os.path.getsize(number_chunks_file) > 0:
            #    with open(number_chunks_file, 'r') as file:
            #        num_chunks = file.read().replace('\n', '')
            #        num_chunks = int(num_chunks)
            #else:
            #    raise Exception(str(number_chunks_file) + " file does not exist\nPlease run exonerate step before running array jobs (blastn, diamond-blastp, hmmscan, rpsblast etc.\n")
            
            #for i in range(num_chunks):
            #    infiles.append(os.path.join(chunks_dir, "obs_filtered.fasta_chunk_{:07d}".format(i)))
            #    outfiles.append(os.path.join(blast_dir, "blastn_chunk_{:07d}.tsv".format(i)))
            #    dones.append(os.path.join(blast_dir, "blastn_chunk_{:07d}.done".format(i)))
         
            job = mirna.blastn(
                os.path.join("abundance", "obs_filtered.fasta"),
                os.path.join("annotations", "blastn_mirbase_mature.tsv"),
                config.param("blastn", 'db_mirbase_mature', 1, 'string')
            )
            job.name = "blastn_mirbase_mature"
            job.subname = "blastn"
            jobs.append(job)
          
            # Merge output chunks
            #job = shotgun_metagenomics.merge_chunks(
            #    blast_dir,
            #    os.path.join("annotations", "blastn_mirbase.tsv"),
            #    num_chunks,
            #    "blastn" 
            #)
            #job.name = "blastn_mirbase_mature_merge"
            #job.subname = "blastn_merge"
            #jobs.append(job)

            job = mirna.filter_blast_hits(
                os.path.join("annotations", "blastn_mirbase_mature.tsv"),
                os.path.join("annotations", "blastn_mirbase_mature_filtered.tsv")
                #os.path.join("annotations", "blastn_mirbase_mature_filtered_best.tsv")
            )
            job.name = "filter_blastn_table_mirna_mature"
            job.subname = "filter"
            jobs.append(job)
        
        else:
            sys.stderr.write("*** Step not implemented and necessary for denovo or denovo-sub worflows. ***\n")
        
        return jobs
    
    # TODO: only for denovo or denovo-sub
    def identify_mirbase_hairpin(self):
        """
        Step blastn mirna against DB(s).
        """
        if not os.path.exists(os.path.join("annotations", "blastn_mirbase_hairpin")):
            os.makedirs(os.path.join("annotations", "blastn_mirbase_hairpin"))
        
        jobs = []
        
        # Do blastn on nt for big assembly  
        chunks_dir = os.path.join("abundance", "fasta_chunks")
        blast_dir = os.path.join("annotations", "blastn_mirbase_hairpin")
        number_chunks_file = os.path.join("abundance", "estimated_number_of_chunks.txt")
        
        infiles = []
        outfiles = []
        dones = []
        
        #num_chunks = 0
        #if os.path.exists(number_chunks_file) and os.path.getsize(number_chunks_file) > 0:
        #    with open(number_chunks_file, 'r') as file:
        #        num_chunks = file.read().replace('\n', '')
        #        num_chunks = int(num_chunks)
        #else:
        #    raise Exception(str(number_chunks_file) + " file does not exist\nPlease run exonerate step before running array jobs (blastn, diamond-blastp, hmmscan, rpsblast etc.\n")
        #
        #for i in range(num_chunks):
        #    infiles.append(os.path.join(chunks_dir, "obs_filtered.fasta_chunk_{:07d}".format(i)))
        #    outfiles.append(os.path.join(blast_dir, "blastn_chunk_{:07d}.tsv".format(i)))
        #    dones.append(os.path.join(blast_dir, "blastn_chunk_{:07d}.done".format(i)))
     
        job = mirna.blastn(
            os.path.join("abundance", "obs_filtered.fasta"),
            os.path.join("annotations", "blastn_mirbase_hairpin.tsv"),
            config.param("blastn", 'db_mirbase_hairpin', 1, 'string')
        )
        job.name = "blastn_mirbase_hairpin"
        job.subname = "blastn"
        jobs.append(job)
      
        # Merge output chunks
        #job = shotgun_metagenomics.merge_chunks(
        #    blast_dir,
        #    os.path.join("annotations", "blastn_mirbase.tsv"),
        #    num_chunks,
        #    "blastn" 
        #)
        #job.name = "blastn_mirbase_mature_merge"
        #job.subname = "blastn_merge"
        #jobs.append(job)

        job = mirna.filter_blast_hits(
            os.path.join("annotations", "blastn_mirbase_hairpin.tsv"),
            os.path.join("annotations", "blastn_mirbase_hairpin_filtered.tsv"),
            os.path.join("annotations", "blastn_mirbase_hairpin_filtered_best.tsv")
        )
        job.name = "filter_blastn_table_mirna_hairpin"
        job.subname = "filter"
        jobs.append(job)



        return jobs
    
    def identify_targets(self):
        """
        Step blastn mirna against DB(s).
        """
        jobs = []
        
        if not os.path.exists(os.path.join("annotations", "ssearch", "fwd_strand")):
            os.makedirs(os.path.join("annotations", "ssearch", "fwd_strand"))
        if not os.path.exists(os.path.join("annotations", "ssearch", "rev_strand")):
            os.makedirs(os.path.join("annotations", "ssearch", "rev_strand"))
        #if not os.path.exists(os.path.join("annotations", "blastn")):
        #    os.makedirs(os.path.join("annotations", "blastn"))

        infile_query = ""
        if(self._type == "reference"):
            infile_query_rc = os.path.join("abundance", "miRNA_abundance_rc.fasta")
            infile_query = os.path.join("abundance", "miRNA_abundance.fasta")

        elif(self._type == "denovo" or self._type == "denovo_sub"):
            #infile_query_rc = os.path.join("annotations", "blastn_mirbase_mature_filtered_best_revcomp.fasta")
            infile_query_rc = os.path.join("abundance", "obs_filtered_rc.fasta")



        db_prefix = config.param("ssearch", "ncbi_genomes_prefix", 1, "string")
        for i in range(0, config.param("ssearch", "ncbi_genomes_chunks", 1, 'posint')):
            #sys.stderr.write("loop i: " + str(i) + "\n")
            infile_db = os.path.join(db_prefix + "{:07d}".format(i))

            # For fwd strand. rcomped miRNA must be used as input
            job = mirna.ssearch(
                infile_query_rc,
                infile_db,
                os.path.join("annotations", "ssearch", "fwd_strand", "mirna_mature_targets_" + "{:07d}.tsv".format(i))
            )
            job.name = "ssearch_fwd_" + "{:07d}".format(i)
            job.subname = "ssearch"
            jobs.append(job)
        
        # Then do reverse strand
        for i in range(0, config.param("ssearch", "ncbi_genomes_chunks", 1, 'posint')):
            #sys.stderr.write("loop i: " + str(i) + "\n")
            infile_db_dir = os.path.join(db_prefix + "{:07d}".format(i))

            # For rev strand, normal (non-revcomped)miRNA must be used vs genomes fwd strands
            job = mirna.ssearch(
                infile_query,
                infile_db_dir,
                os.path.join("annotations", "ssearch", "rev_strand", "mirna_mature_targets_" + "{:07d}.tsv".format(i))
            )
            job.name = "ssearch_rev_" + "{:07d}".format(i)
            job.subname = "ssearch"
            jobs.append(job)
        
            #job = mirna.blastn(
            #    infile_query,
            #    os.path.join("annotations", "blastn", "mirna_mature_targets_" + "{:07d}.tsv".format(i)),
            #    infile_db
            #)
            #job.name = "blastn_identify_ncbi_genomes_" + "{:07d}".format(i)
            #job.subname = "blastn"
            #jobs.append(job)

        outfiles_fwd = []
        outfiles_rev = []
        for i in range(0, config.param("ssearch", "ncbi_genomes_chunks", 1, 'posint')):

            job = mirna.parse_mirna_targets(
                os.path.join("annotations", "ssearch", "fwd_strand", "mirna_mature_targets_" + "{:07d}.tsv".format(i)),
                os.path.join("annotations", "ssearch", "fwd_strand", "mirna_mature_targets_" + "{:07d}_raw.tsv".format(i)),
                os.path.join("annotations", "ssearch", "fwd_strand", "mirna_mature_targets_" + "{:07d}_parsed.tsv".format(i)),
                "fwd_strand"
            )
            job.name = "parse_mirna_targets_fwd_" + "{:07d}".format(i)
            job.subname = "parse_mirna_targets"
            jobs.append(job)
            
            outfiles_fwd.append(os.path.join("annotations", "ssearch", "fwd_strand", "mirna_mature_targets_" + "{:07d}_parsed.tsv".format(i)))
        
            job = mirna.parse_mirna_targets(
                os.path.join("annotations", "ssearch", "rev_strand", "mirna_mature_targets_" + "{:07d}.tsv".format(i)),
                os.path.join("annotations", "ssearch", "rev_strand", "mirna_mature_targets_" + "{:07d}_raw.tsv".format(i)),
                os.path.join("annotations", "ssearch", "rev_strand", "mirna_mature_targets_" + "{:07d}_parsed.tsv".format(i)),
                "rev_strand"
            )
            job.name = "parse_mirna_targets_rev_" + "{:07d}".format(i)
            job.subname = "parse_mirna_targets"
            jobs.append(job)

            outfiles_rev.append(os.path.join("annotations", "ssearch", "rev_strand", "mirna_mature_targets_" + "{:07d}_parsed.tsv".format(i)))
        
        # Merge output chunks
        job = mirna.merge_chunks_generic(
            outfiles_fwd,
            os.path.join("annotations", "ssearch", "mirna_mature_targets_fwd_parsed.tsv")
        )
        job.name = "ssearch_merge_fwd"
        job.subname = "ssearch_merge"
        jobs.append(job)
        
        job = mirna.merge_chunks_generic(
            outfiles_rev,
            os.path.join("annotations", "ssearch", "mirna_mature_targets_rev_parsed.tsv")
        )
        job.name = "ssearch_merge_rev"
        job.subname = "ssearch_merge"
        jobs.append(job)

        job = mirna.merge_fwd_rev_targets(
            os.path.join("annotations", "ssearch", "mirna_mature_targets_fwd_parsed.tsv"),
            os.path.join("annotations", "ssearch", "mirna_mature_targets_rev_parsed.tsv"),
            os.path.join("annotations", "ssearch", "mirna_mature_targets_parsed.tsv")
        )
        job.name = "merge_mirna_targets_fwd_and_rev"
        job.subname = "parse_mirna_target"
        jobs.append(job)
        
        #job = mirna.filter_blast_hits_for_targets_finding(
        #    os.path.join("annotations", "ssearch", "mirna_mature_targets.tsv"),
        #    os.path.join("annotations", "ssearch", "mirna_mature_targets_filtered.tsv")
        #)
        #job.name = "filter_blastn_table_mirna_hairpin"
        #job.subname = "filter"
        #jobs.append(job)
           #
            #job = mirna.binding_site_prediction(
            #    infile,
            #    os.path.join("annotations", "mirna_mature_targets.tsv"),
            #    os.path.join("annotations", "mirna_mature_targets.tsv")
            #)
            #job.name = "viennaRNA_RNAup"
            #job.subname = "rnaup"
            #jobs.append(job)

        return jobs

    def finalize(self):
        jobs = []

        if(self._type == "denovo" or self._type == "denovo_sub"):

            job = mirna.merge_tables(
                os.path.join("abundance", "obs_filtered.tsv"),
                #os.path.join("annotations", "blastn_mirbase_hairpin_filtered_best.tsv"),
                os.path.join("annotations", "blastn_mirbase_mature_filtered.tsv"),
                os.path.join("annotations", "mirna_table.tsv")
            )
            job.name = "merge_mirna_tables"
            job.subname = "merge"
            jobs.append(job)
            
            job = mirna.normalize_mirna_table(
                os.path.join("annotations", "mirna_table.tsv"),
                os.path.join("abundance", "obs_filtered.fasta"),
                os.path.join("annotations", "mirna_table_rpkm.tsv")
            )
            job.name = "normalize_miRNA_matrix"
            job.name = "normalize"
            jobs.append(job)

        return jobs
    
    def alpha_diversity(self):
        jobs = []
        
        # Rarefaction curves using all data - get an idea if sequencing depth is appropriate for each sample.
        # Note if really applicable to this type of data though...
        job = mirna.rtk(
            os.path.join("abundance", "miRNA_abundance.tsv"),
            os.path.join("alpha_div"),
            "cluster",
            remove_last_col = True
        )
        job.name = "alpha_diversity_saturation"
        job.subname = "alpha_diversity"
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
            self.trim,
            self.remove_contam,
            self.filter_by_size,
            self.align,
            self.abundance,
            #self.exonerate,
            self.identify_mirbase_mature,
            #self.identify_mirbase_hairpin,
            self.identify_targets,
            self.finalize,
            self.alpha_diversity
            #self.cleanup
        ]

    def set_local_variables(self):
        self._parser_local = self.argparser

        # barcodes
        self._args_local = self._parser_local.parse_args()
        config.parse_files(self._args_local.config)
        self._config = self._args_local.config[0]
        self._config = os.path.abspath(self._config.name)
        self._type = self._args_local.type
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
        mkdir_p(os.path.join(root_dir, "qced_reads"))

        #prefixes = ["5S", "8S", "16S", "18S", "23S", "28S"]
    
    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_readset_file(self.args.readsets.name)
        return self._readsets


    def __init__(self):
        version = open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "..", "VERSION"), 'r').read().split('\n')[0]
        pipeline_string = """
###############################################################################
                             _ _____  _   _          
                            (_)  __ \| \ | |   /\    
                   _ __ ___  _| |__) |  \| |  /  \   
                  | '_ ` _ \| |  _  /| . ` | / /\ \  
                  | | | | | | | | \ \| |\  |/ ____ \ 
                  |_| |_| |_|_|_|  \_\_| \_/_/    \_\ 
                                     
               Support: jtremblay514@gmail.com
             Home page: jtremblay.github.io/pipelines.html
               Version: """ + version + """

###############################################################################"""
        sys.stderr.write(pipeline_string + '\n')
        time.sleep(1)
        # Add pipeline specific arguments
        self.argparser.add_argument("-r", "--readsets", help="readset file",  type=argparse.FileType('r'), required=False)
        self.argparser.add_argument("-t", "--type", help="denovo, denovo_sub or reference", type=None, required=True)
        self.set_local_variables()
        sys.stderr.write('miRNA pipeline\n')
        super(MiRNA, self).__init__()
                
MiRNA().submit_jobs()
