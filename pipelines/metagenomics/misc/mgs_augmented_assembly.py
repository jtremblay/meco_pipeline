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

class MGS_augmented_assembly(illumina.Illumina):

    def trim(self):
        jobs = []
        # Merge all demultiplexed fastq files in one file. One file for reads1 and one file for reads2 of Illumina paired-end.
        # If library is Illumina single end or 454 or PacBio, only one file will be generated.
        outdir = self._root_dir
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
            elif readset.run_type == "SINGLE_END":
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (SINGLE_END not implemented yet)!")
            
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
     
    def first_assembly(self):
        root = self._root_dir
        jobs = []
        fastq_list = []
        contigs_list = []
        contigs_out_list = []
        
        for readset in self.readsets:
            infile = os.path.join(root, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
            outdir = os.path.join(root, "first_assembly", readset.sample.name) 
            curr_contigs = os.path.join(root, "first_assembly", readset.sample.name, "Contigs.fasta") 
            curr_contigs_out = os.path.join(root, "first_assembly", readset.sample.name, "Contigs_renamed.fasta") 
            fastq_list.append(infile)
            contigs_list.append(curr_contigs)
            contigs_out_list.append(curr_contigs)

            job = shotgun_metagenomics.ray(
                infile,
                outdir,
                "pe"
            )
            job.name = "ray_" + readset.sample.name
            job.subname = "ray"
            jobs.append(job)
             
        job = shotgun_metagenomics.ray_big(
            fastq_list,
            os.path.join(root, "first_assembly", "big_assembly"),
            "pe"
        )
        job.name = "ray_big_assembly"
        job.subname = "ray_big"
        jobs.append(job)

        #job = shotgun_metagenomics.rename_contigs(
        #    contigs_list,
        #    contigs_out_list
        #)
        #job.name = "rename_contigs"
        #job.subname = "rename_contigs"
        #jobs.append(job)

        return jobs

    def gene_prediction(self):
        root = self._root_dir
        jobs = []
        gff_list = []
        fna_list = []
        
        for readset in self.readsets:
            infile = os.path.join(root, "first_assembly", readset.sample.name, "Contigs_renamed.fasta")
            outdir = os.path.join(root, "gene_prediction", readset.sample.name)
            gff = os.path.join(root, "gene_prediction", readset.sample.name, readset.sample.name + ".gff")
            fna = os.path.join(root, "gene_prediction", readset.sample.name, readset.sample.name + ".fna")
        
            gff_list.append(gff)
            fna_list.append(fna)

            if not os.path.exists(outdir):
                os.makedirs(os.path.join(outdir))
        
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
         
        return jobs

    def gene_clustering(self):
        jobs = []
        
        infiles = []
        for readset in self.readsets:
            infiles.append(os.path.join(self._root_dir, "gene_prediction", readset.sample.name, readset.sample.name + "_renamed.fna"))

        merged_fasta = os.path.join(self._root_dir, "gene_clustering", "genes.fasta")
        barcodes = os.path.join(self._root_dir, "gene_clustering", "barcodes.fasta")

        job = shotgun_metagenomics.preprocess_for_clustering(
            infiles,
            merged_fasta,
            barcodes
        )
        job.name = "preprocess_for_clustering"
        job.subname = "preprocess_for_clustering"
        jobs.append(job)
        
        job = shotgun_metagenomics.clustering(
            merged_fasta,
            barcodes,
            os.path.join(self._root_dir, "gene_clustering")
        )   
        job.name = "generate_clusters"
        job.subname = "clustering"
        jobs.append(job)
        
        job = shotgun_metagenomics.filter_clusters(
            os.path.join(self._root_dir, "gene_clustering", "cdhit.fasta"),
            os.path.join(self._root_dir, "gene_clustering", "cdhit_filtered.fasta")
        )   
        job.name = "filter_clusters"
        job.subname = "filter_clusters"
        jobs.append(job)
        
        return jobs

    ## Gene abundance and canopy clustering
    def gene_abundance(self):
        jobs = []
        cov_list = []

        # Make index and bed file.
        job = shotgun_metagenomics.make_index(
            os.path.join(self._root_dir, "gene_clustering", "cdhit_filtered.fasta"),
        )
        job.name = "make_index"
        job.subname = "make_index"
        jobs.append(job)

        for readset in self.readsets:
            infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
            out1 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.sample.name + ".ncontam_paired_R1.fastq.gz")
            out2 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.sample.name + ".ncontam_paired_R2.fastq.gz")
            bam = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.sample.name + ".bam")
            cov = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.sample.name + ".cov")
            cov_list.append(cov)

            outdir = self._root_dir + "gene_abundance/" + readset.sample.name
            if not os.path.exists(outdir):
                os.makedirs(os.path.join(outdir))
            
            job = shotgun_metagenomics.split_pairs(
                infile,
                out1,
                out2
            )
            job.name = "split_pairs_" + readset.sample.name
            job.subname = "split_pairs"
            jobs.append(job)

            job = shotgun_metagenomics.bwa_mem(
                os.path.join(self._root_dir, "gene_clustering", "cdhit_filtered.fasta"),
                out1,
                out2,
                bam,
                os.path.join(self._root_dir, "gene_clustering", "cdhit_filtered.fasta.bwt"),
            )
            job.name = "bwa_mem-" + readset.sample.name
            job.subname = "bwa"
            jobs.append(job)
            
            job = shotgun_metagenomics.coverage_bed(
                bam,
                self._root_dir + "gene_clustering/cdhit_filtered.bed",
                cov
            )
            job.name = "bedtoolsCov-" + readset.sample.name
            job.subname = "bedtools"
            jobs.append(job)
            
        # Once all coverage has been computed, merge all tables.
        job = shotgun_metagenomics.merge_counts(
            cov_list,
            self._root_dir + "gene_abundance/merged_gene_abundance.tsv"
        )
        job.name = "merge_gene_abundance"
        job.subname = "merge_gene_abundance"
        jobs.append(job)
         
        return jobs
    
    def gene_binning(self):
        jobs = []
        
        job = shotgun_metagenomics.canopy(
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            "canopy",
            os.path.join(self._root_dir, "gene_binning", "clusters.tsv"),
            os.path.join(self._root_dir, "gene_binning", "profiles.tsv")
        )
        job.name = "canopy"
        job.subname = "canopy"
        jobs.append(job)

        contigs = []
        for readset in self.readsets:
            curr_contigs = os.path.join(self._root_dir, "first_assembly", readset.sample.name, "Contigs_renamed.fasta")
            contigs.append(curr_contigs)
         
        if not os.path.exists(os.path.join(self._root_dir + "gene_binning", "canopy_contigs")):
            os.makedirs(os.path.join(self._root_dir, "gene_binning", "canopy_contigs"))
       
        # Generate one fasta file per canopy contigs and make bwa index.
        job = shotgun_metagenomics.extract_canopy_contigs(
            contigs,
            os.path.join(self._root_dir, "gene_binning", "canopy_contigs"),
            os.path.join(self._root_dir, "gene_binning", "clusters.tsv")
        )
        job.name = "extract_canopy_contigs"
        job.subname = "extract_canopy_contigs"
        jobs.append(job)
        
        job = shotgun_metagenomics.parse_canopies(
            os.path.join(self._root_dir, "gene_binning", "clusters.tsv"),
            os.path.join(self._root_dir, "gene_binning", "clusters_parsed.tsv")
        )
        job.name = "parse_canopies"
        job.subname = "parse_canopies"
        jobs.append(job)

        bams = []
        for readset in self.readsets:
            reads1 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
            reads2 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
            bam = os.path.join(self._root_dir, "gene_binning", "canopy_contigs", "canopy." + readset.sample.name + ".bam")
            bams.append(bam)

            job = shotgun_metagenomics.bwa_mem_canopies(
                os.path.join(self._root_dir, "gene_binning", "canopy_contigs", "canopy_contigs.fasta"),
                reads1,
                reads2,
                bam
            )
            job.name = "map_raw_reads_on_canopy_contigs" + readset.name
            job.subname = "bwa"
            jobs.append(job)
        
        job = shotgun_metagenomics.extract_reads1(
            bams,
            os.path.join(self._root_dir, "gene_binning", "reads"),
            os.path.join(self._root_dir, "gene_binning", "clusters_parsed.tsv")
        )
        job.name = "extract_reads_1"
        job.subname = "extract_reads"
        jobs.append(job)
        
        job = shotgun_metagenomics.extract_reads2(
            bams,
            os.path.join(self._root_dir, "gene_binning", "reads"),
            os.path.join(self._root_dir, "gene_binning", "clusters_parsed.tsv")
        )
        job.name = "extract_reads_2"
        job.subname = "extract_reads"
        jobs.append(job)

        return jobs
    
    def mgs_augmented_assembly(self):
        jobs = []

        range = config.param("ray", "range", 1, "string")
        range_list = range.split(":")

        for curr_range in range_list:
            indir = os.path.join(self._root_dir, "gene_binning", "reads")
            outdir = os.path.join(self._root_dir, "mgs_augmented_assembly")
            dummy_infile1 = os.path.join(self._root_dir, "gene_binning", "extract_reads1.done")
            dummy_infile2 = os.path.join(self._root_dir, "gene_binning", "extract_reads2.done")
       
            job = shotgun_metagenomics.ray_wrapper(
                dummy_infile1,
                dummy_infile2,
                indir,
                outdir,
                curr_range
            )
            job.name = "ray_wrapper_" + curr_range
            job.subname = "ray"
            jobs.append(job)
        
        return jobs
    
    def taxonomic_annotation(self):
        jobs = []
        #job = shotgun_metagenomics.mymethod(
        #)
        #job.name = "myjobname"
        #jobs.append(job)
        sys.stderr.write('not implemented yet\n')
        return jobs
    
    def gene_annotation(self):
        jobs = []
        #job = shotgun_metagenomics.mymethod(
        #)
        #job.name = "myjobname"
        #jobs.append(job)
        sys.stderr.write('not implemented yet\n')
        return jobs 
     
    def rescaffold(self):
        jobs = []
        #module load mugqic_dev/SSPACE/3.0 && SSPACE_Standard_v3.0.pl -l libraries.txt -s ./assembly/Ray_k_31/Scaffolds.fasta -x 0 -m 32 -o 20 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 1 -S 0 -b standard_out_k31
        return jobs

    def statistics(self):
        jobs = []
        #job = shotgun_metagenomics.mymethod(
        #)
        #job.name = "myjobname"
        #jobs.append(job)
        sys.stderr.write('not implemented yet\n')
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
            self.trim,
            self.duk,
            self.first_assembly,
            self.gene_prediction,
            self.gene_clustering,
            self.gene_abundance,
            self.gene_binning,
            self.mgs_augmented_assembly,
            self.taxonomic_annotation,
            self.gene_annotation,
            #self.gene_abundance_clustering,
            self.statistics,
            self.generate_report,
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
        if not os.path.exists(root_dir + "/qced_reads"):
            os.makedirs(os.path.join(root_dir, "qced_reads"))
        if not os.path.exists(root_dir + "/first_assembly"):
            os.makedirs(os.path.join(root_dir, "first_assembly"))
        if not os.path.exists(root_dir + "/gene_prediction"):
            os.makedirs(os.path.join(root_dir, "gene_prediction"))
        if not os.path.exists(root_dir + "/gene_clustering"):
            os.makedirs(os.path.join(root_dir, "gene_clustering"))
        if not os.path.exists(root_dir + "/gene_abundance"):
            os.makedirs(os.path.join(root_dir, "gene_abundance"))
        if not os.path.exists(root_dir + "/taxonomic_annotation"):
            os.makedirs(os.path.join(root_dir, "taxonomic_annotation"))
        if not os.path.exists(root_dir + "/gene_annotation"):
            os.makedirs(os.path.join(root_dir, "gene_annotation"))
        if not os.path.exists(root_dir + "/gene_abundance_clustering"):
            os.makedirs(os.path.join(root_dir, "gene_abundance_clustering"))
        if not os.path.exists(root_dir + "/mgs_augmented_assembly"):
            os.makedirs(os.path.join(root_dir, "mgs_augmented_assembly"))
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
        sys.stderr.write('Running MGS (Metagenomic Genomic Specie) augmented pipeline\n')
        super(MGS_augmented_assembly, self).__init__()
                
MGS_augmented_assembly().submit_jobs()
