#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys
#from os.path import basename

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

class Blast(illumina.Illumina):
 
    def exonerate(self):
        jobs = []

        # Split for big assembly
        #infile_fna = os.path.join(self._root_dir, "gene_prediction", "big_assembly", "big_assembly.fna")
        infile_fna = self._infile
        chunks_dir = os.path.join(self._root_dir, "blastn_nt", "fasta_chunks")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
        base = os.path.basename(self._infile)
        prefix = os.path.splitext(base)[0]
        
        #FNA
        job = shotgun_metagenomics.exonerate(
            infile_fna,
            chunks_dir,
            num_chunks,
            base 
        )
        job.name = "exonerate_blast"
        job.subname = "exonerate"
        jobs.append(job)
        
        return jobs    

    def blastn_nt(self):
        jobs = []
        
        # Do blastn on nt for big assembly  
        chunks_dir = os.path.join(self._root_dir, "blastn_nt", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "blastn_nt")
        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint')
        base = os.path.basename(self._infile)
        prefix = os.path.splitext(base)[0]

        for i in range(num_chunks):
            job = shotgun_metagenomics.blastn(
                os.path.join(chunks_dir, base + "_chunk_{:07d}".format(i)),
                os.path.join(blast_dir, "blastn_chunk_{:07d}.tsv".format(i)),
                blast_dir 
            )
            job.name = "blastn_nt"
            job.subname = "blastn"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join(self._root_dir, "blastn_nt", "blastn_nt.tsv"),
            num_chunks,
            "blastn" 
        )
        job.name = "blastn_nt"
        job.subname = "blastn"
        jobs.append(job)
        
        job = shotgun_metagenomics.keep_blast_best_hit(
            os.path.join(self._root_dir, "blastn_nt", "blastn_nt.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "big_assembly", "blastn_nt_besthit.tsv")
        )
        job.name = "blastn_nt_best_hit"
        job.subname = "blastn_best_hit"
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
            self.exonerate,#5
            self.blastn_nt,#6
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
        self._parser_local.add_argument("-i", "--infile-fasta", help="fasta file file", type=file, required=True)

        # barcodes
        self._args_local = self._parser_local.parse_args()
        config.parse_files(self._args_local.config)
        self._config = self._args_local.config[0]
        self._config = os.path.abspath(self._config.name)
        self._infile = self._args_local.infile_fasta
        self._infile = os.path.abspath(self._infile.name)
        #self._infile = self._args_local.infile_fasta
        sys.stderr.write('self._infile: ' + self._infile + '\n')
        
        self._root_dir = self._args_local.output_dir
        if not os.path.exists(self._root_dir):
            os.makedirs(self._root_dir)
        
        # Make directories
        self.make_directories(self._root_dir)
  
    # Define and make directories. Also desing initial infile.
    def make_directories(self, root_dir):
         
        if not os.path.exists(root_dir):
            os.makedirs(root_dir)
        if not os.path.exists(root_dir + "/blastn_nt"):
            os.makedirs(os.path.join(root_dir, "blastn_nt"))


    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            #self._readsets = parse_nanuq_readset_file(self.args.readsets.name)
            self._readsets = parse_readset_file(self.args.readsets.name)
        return self._readsets


    def __init__(self):
        # Add pipeline specific arguments
        self.set_local_variables()
        self.argparser.add_argument("-i", "--infile-fasta", help="infile in fasta format", type=file)
        sys.stderr.write('Running blast pipeline\n')
        super(Blast, self).__init__()
                
Blast().submit_jobs()
