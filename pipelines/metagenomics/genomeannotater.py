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
from bio import genome_assembly
from bio import microbial_ecology

from pipelines import common

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# Global scope variables
log = logging.getLogger(__name__)

class GenomeAnnotater(common.CAFPipeline):
    
    def gene_prediction(self):
        
        """
        Step gene_prediction(): Use Prodigal to predict ORFs of bacterial genomes.
        """
        jobs = []
        
        for bin in self.magsets:   
            #log.info("MAG name: " + bin.name)
            #log.info("MAG file: " + bin.fasta_file)
            if not os.path.exists(os.path.join(bin.name)):
                os.makedirs(os.path.join(bin.name))
            if not os.path.exists(os.path.join(bin.name, "gene_prediction")):
                os.makedirs(os.path.join(bin.name, "gene_prediction"))

            infile = bin.fasta_file
            outdir = os.path.join(bin.name, "gene_prediction")
            
            if not os.path.exists(outdir):
                os.makedirs(os.path.join(outdir))
            
            job = genome_assembly.prodigal(
                infile,
                os.path.join(outdir, "Contigs.gff"),
                os.path.join(outdir, "Contigs.fna"),
                os.path.join(outdir, "Contigs.faa"),
                os.path.join(outdir, "Contigs_renamed.gff"),
                os.path.join(outdir, "Contigs_renamed.fna"),
                os.path.join(outdir, "Contigs_renamed.faa")
            )
            job.name = bin.name + "_" + "prodigal"
            job.subname = "prodigal"
            jobs.append(job)
     
        return jobs

    def exonerate(self):
        
        """
        Step exonerate(): Split co-assembly and predicted genes fasta files in smaller chunks to perform downstream annotations.
        """
        jobs = []
        
        for bin in self.magsets:   
            #log.info("MAG name: " + bin.name)
            #log.info("MAG file: " + bin.fasta_file)
            if not os.path.exists(os.path.join(bin.name)):
                os.makedirs(os.path.join(bin.name))
            if not os.path.exists(os.path.join(bin.name, "assembly", "fasta_chunks")):
                os.makedirs(os.path.join(bin.name, "assembly", "fasta_chunks"))
            if not os.path.exists(os.path.join(bin.name, "gene_prediction", "fasta_chunks")):
                os.makedirs(os.path.join(bin.name, "gene_prediction", "fasta_chunks"))
                
            if not os.path.lexists(os.path.join(bin.name, "assembly", "Contigs.fasta")):
                os.symlink(os.path.join("..", "..", bin.fasta_file), os.path.join(bin.name, "assembly", "Contigs.fasta"))

            infile_contigs = os.path.join(bin.name, "assembly", "Contigs.fasta")

            chunks_dir_contigs = os.path.join(bin.name, "assembly", "fasta_chunks")
            number_of_chunks_file_contigs = os.path.join(bin.name, "assembly", "estimated_number_of_chunks_contigs.txt")
            infile_faa = os.path.join(bin.name, "gene_prediction", "Contigs_renamed.faa")
            chunks_dir = os.path.join(bin.name, "gene_prediction", "fasta_chunks")
            number_of_chunks_file_genes = os.path.join(bin.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
            
            #FNA contigs assembly (i.e. contigs)
            job = shotgun_metagenomics.estimate_number_of_chunks(
                infile_contigs,
                number_of_chunks_file_contigs,
                str(config.param('exonerate', 'targeted_chunk_file_size_contigs', 1, 'int'))
            )
            job.name = bin.name + "_" + "estimate_chunks_file_size"
            jobs.append(job)

            job = shotgun_metagenomics.exonerate(
                infile_contigs,
                chunks_dir_contigs,
                number_of_chunks_file_contigs,
                "Contigs.fasta" 
            )
            job.name = bin.name + "_" + "exonerate_contigs_fna"
            job.subname = "exonerate"
            jobs.append(job)
            
            #FAA genes (i.e. predicted genes)
            job = shotgun_metagenomics.estimate_number_of_chunks(
                infile_faa,
                number_of_chunks_file_genes,
                str(config.param('exonerate', 'targeted_chunk_file_size_genes', 1, 'int'))
            )
            job.name = bin.name + "_" + "estimate_chunks_file_size"
            jobs.append(job)
            
            job = shotgun_metagenomics.exonerate(
                infile_faa,
                chunks_dir,
                number_of_chunks_file_genes,
                "Contigs_renamed.faa" 
            )
            job.name = bin.name + "_" + "exonerate_genes_fna"
            job.subname = "exonerate"
            jobs.append(job)
        
        return jobs    
    
    def diamond_blastp_kegg(self):
        """
        Step diamond_blastp_kegg(): DIAMOND BLASTp of predicted genes vs KEGG db.
        """

        jobs = []
        
        for bin in self.magsets:   
            #log.info("MAG name: " + bin.name)
            #log.info("MAG file: " + bin.fasta_file)
            if not os.path.exists(os.path.join(bin.name)):
                os.makedirs(os.path.join(bin.name))
            if not os.path.exists(os.path.join(bin.name, "annotations", "blastp_kegg")):
                os.makedirs(os.path.join(bin.name, "annotations", "blastp_kegg"))
            
            # Do blastp on KEGG for co-assembly 
            chunks_dir = os.path.join(bin.name, "gene_prediction", "fasta_chunks")
            blast_dir = os.path.join(bin.name, "annotations", "blastp_kegg")
            number_chunks_file = os.path.join(bin.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
            
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
            job.name = bin.name + "_" + "blastp_array_kegg"
            job.subname = "blastp_kegg"
            job.job_array_num_task = num_chunks
            jobs.append(job)
            
            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks(
                blast_dir,
                os.path.join(bin.name, "annotations", "blastp_kegg.tsv"),
                num_chunks,
                "blastp"
            )
            job.name = bin.name + "_" + "blastp_kegg_merge"
            job.subname = "merge"
            jobs.append(job)
            
            # Generate a clean table of module/pathways.
            job = shotgun_metagenomics.parse_kegg(
                os.path.join(bin.name, "annotations", "blastp_kegg.tsv"),
                os.path.join(bin.name, "annotations", "blastp_kegg_parsed.tsv")
            )
            job.name = bin.name + "_" + "parse_kegg"
            job.subname = "merge"
            jobs.append(job)

        return jobs 
    
    def ncrna(self):
        """
        Step ncrna(): scan for non-coding RNA in contigs: rRNA and tRNA sequence.
        """

        jobs = []
        
        for bin in self.magsets:   
            #log.info("MAG name: " + bin.name)
            #log.info("MAG file: " + bin.fasta_file)
            if not os.path.exists(os.path.join(bin.name)):
                os.makedirs(os.path.join(bin.name))
            if not os.path.exists(os.path.join(bin.name, "annotations", "trna")):
                os.makedirs(os.path.join(bin.name, "annotations", "trna"))
            if not os.path.exists(os.path.join(bin.name, "annotations", "rrna", "barrnap")):
                os.makedirs(os.path.join(bin.name, "annotations", "rrna", "barrnap"))
                    
            # Do tRNA of co-assembly 
            chunks_dir = os.path.join(bin.name, "assembly", "fasta_chunks")
            outdir = os.path.join(bin.name, "annotations", "trna")
            number_chunks_file = os.path.join(bin.name, "assembly", "estimated_number_of_chunks_contigs.txt")
            
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
         
            job = shotgun_metagenomics.trnascanse_array_job(
                os.path.join(chunks_dir),
                "Contigs.fasta_chunk_",
                os.path.join(outdir),
                "trna_chunk_",
                "trnascanse",
                infiles,
                outfiles,
                dones,
                self._curr_scheduler
            )
            job.name = bin.name + "_" + "trna_array_job"
            job.subname = "trnascanse"
            job.job_array_num_task = num_chunks
            jobs.append(job)
            
            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks(
                outdir,
                os.path.join(bin.name, "annotations", "trna.tsv"),
                num_chunks,
                "trna"
            )
            job.name = bin.name + "_" + "merge_tRNAScan-SE"
            job.subname = "merge"
            jobs.append(job)

            # Do rRNA prediction from co-assembly with barrnap
            prefixes = ["bac", "arc", "euk", "mito"]
            for prefix in prefixes:
                job = shotgun_metagenomics.barrnap(
                    os.path.join(bin.name, "assembly", "Contigs.fasta"),
                    prefix,
                    os.path.join(bin.name, "annotations", "rrna", "barrnap", prefix + ".fna"),
                    os.path.join(bin.name, "annotations", "rrna", "barrnap", prefix + ".tsv")
                )
                job.name = bin.name + "_" + "barrnap_" + prefix
                job.subname = "barrnap"
                jobs.append(job)
           
                job = shotgun_metagenomics.split_barrnap_fna(
                    os.path.join(bin.name, "annotations", "rrna", "barrnap", prefix + ".fna"),
                    prefix,
                    os.path.join(bin.name, "annotations", "rrna", "barrnap", prefix + "_5S.fna"), 
                    os.path.join(bin.name, "annotations", "rrna", "barrnap", prefix + "_SSU.fna"), 
                    os.path.join(bin.name, "annotations", "rrna", "barrnap", prefix + "_LSU.fna"), 
                    os.path.join(bin.name, "annotations", "rrna", "barrnap", prefix + "_12S.fna"), 
                    os.path.join(bin.name, "annotations", "rrna", "barrnap", prefix + "_5-8S.fna") 
                )
                job.name = bin.name + "_" + "split_barrnap_fna_" + prefix
                jobs.append(job)
        
            prefixes = ["arc_SSU", "arc_LSU", "bac_SSU", "bac_LSU", "euk_SSU", "euk_LSU"]
            for prefix in prefixes:
                # create bed files for abundance - we have to only get the rea
                # Just take taxonomy file which contains coord for each rRNA genes, then generate bed file.
                job = shotgun_metagenomics.barrnap_fna_to_bed(
                    os.path.join(bin.name, "annotations", "rrna", "barrnap", prefix + ".fna"),
                    os.path.join(bin.name, "annotations", "rrna", "barrnap", prefix + ".bed")
                )
                job.name = bin.name + "_" + "barrnap_to_bed_" + prefix
                job.subname = "barrnap_to_bed"
                jobs.append(job)

        return jobs
    
    def rpsblast_cog(self):
        """
        Step rpsblast_cog(): RPSBLAST of predicted genes vs COG db
        """
        
        jobs = []
        
        for bin in self.magsets:   
            #log.info("MAG name: " + bin.name)
            #log.info("MAG file: " + bin.fasta_file)
            if not os.path.exists(os.path.join(bin.name)):
                os.makedirs(os.path.join(bin.name))
            if not os.path.exists(os.path.join(bin.name, "annotations", "rpsblast_cog")):
                os.makedirs(os.path.join(bin.name, "annotations", "rpsblast_cog"))
        
        # Do rpsblast on COG for big assembly  
            if(config.param('DEFAULT', 'skip_cog', 1, 'string') == 'yes'):
                fname = os.path.join(bin.name, "annotations", "rpsblast_cog.tsv")
                open(fname, 'a').close()
                os.utime(fname, None)
            else:
                number_chunks_file = os.path.join(bin.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
            
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
                
                chunks_dir = os.path.join(bin.name, "gene_prediction", "fasta_chunks")
                blast_dir = os.path.join(bin.name, "annotations", "rpsblast_cog")
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
                job.name = bin.name + "_" + "rpsblast_cog"
                job.subname = "rpsblast"
                job.job_array_num_task = num_chunks
                jobs.append(job)
        
                # Merge output chunks
                job = shotgun_metagenomics.merge_chunks(
                    blast_dir,
                    os.path.join(bin.name, "annotations", "rpsblast_cog.tsv"),
                    num_chunks,
                    "rpsblast"
                )
                job.name = bin.name + "_" + "rpsblast_cog"
                job.subname = "merge"
                jobs.append(job)
        
        return jobs 
    
    def rpsblast_kog(self):
        """
        Step rpsblast_kog(): RPSBLAST of predicted genes vs KOG db
        """
        
        jobs = []
        
        for bin in self.magsets:   
            #log.info("MAG name: " + bin.name)
            #log.info("MAG file: " + bin.fasta_file)
            if not os.path.exists(os.path.join(bin.name)):
                os.makedirs(os.path.join(bin.name))
            if not os.path.exists(os.path.join(bin.name, "annotations", "rpsblast_kog")):
                os.makedirs(os.path.join(bin.name, "annotations", "rpsblast_kog"))
        
        # Do rpsblast on KOG for big assembly  
            if(config.param('DEFAULT', 'skip_kog', 1, 'string') == 'yes'):
                fname = os.path.join(bin.name, "annotations", "rpsblast_kog.tsv")
                open(fname, 'a').close()
                os.utime(fname, None)
            else:
                number_chunks_file = os.path.join(bin.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
            
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
                
                chunks_dir = os.path.join(bin.name, "gene_prediction", "fasta_chunks")
                blast_dir = os.path.join(bin.name, "annotations", "rpsblast_kog")
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
                job.name = bin.name + "_" + "rpsblast_kog"
                job.subname = "rpsblast"
                job.job_array_num_task = num_chunks
                jobs.append(job)
        
                # Merge output chunks
                job = shotgun_metagenomics.merge_chunks(
                    blast_dir,
                    os.path.join(bin.name, "annotations", "rpsblast_kog.tsv"),
                    num_chunks,
                    "rpsblast"
                )
                job.name = bin.name + "_" + "rpsblast_kog"
                job.subname = "merge"
                jobs.append(job)
        
        return jobs 
    
    def hmmscan_pfam(self):
        
        """
        Step hmmscan_pfam(): HMMSCAN of predicted genes vs PFAM-A DB.
        """
        
        jobs = []
        
        for bin in self.magsets:   
            #log.info("MAG name: " + bin.name)
            #log.info("MAG file: " + bin.fasta_file)
            if not os.path.exists(os.path.join(bin.name)):
                os.makedirs(os.path.join(bin.name))
            if not os.path.exists(os.path.join(bin.name, "annotations", "rpsblast_kog")):
                os.makedirs(os.path.join(bin.name, "annotations", "rpsblast_kog"))
        
            chunks_dir = os.path.join(bin.name, "gene_prediction", "fasta_chunks")
            hmmscan_out_dir = os.path.join(bin.name, "annotations", "hmmscan_pfam")
            number_chunks_file = os.path.join(bin.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
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
                tblouts.append(os.path.join(hmmscan_out_dir, "hmmscan_chunk_{:07d}.tblout".format(i)))
                domtblouts.append(os.path.join(hmmscan_out_dir, "hmmscan_chunk_{:07d}.domtblout".format(i)))
                pfamtblouts.append(os.path.join(hmmscan_out_dir, "hmmscan_chunk_{:07d}.pfamtblout".format(i)))
                dones.append(os.path.join(hmmscan_out_dir, "hmmscan_chunk_{:07d}.done".format(i)))

            job = shotgun_metagenomics.hmmscan_array_job(
                os.path.join(chunks_dir),
                "Contigs_renamed.faa_chunk_",
                os.path.join(hmmscan_out_dir),
                "hmmscan_chunk_",
                infiles,
                tblouts, domtblouts, pfamtblouts,
                dones,
                config.param('pfam', 'db', required=True),
                self._curr_scheduler
            )
            job.name = bin.name + "_" + "hmmscan_pfam"
            job.subname = "hmmscan"
            job.job_array_num_task = num_chunks
            jobs.append(job)

            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks_hmms(
                hmmscan_out_dir,
                os.path.join(bin.name, "annotations"),
                num_chunks,
                "hmmscan",
                "hmmscan_pfam"
            )
            job.name = bin.name + "_" + "hmmscan_pfam_merge"
            job.subname = "merge"      
            jobs.append(job)
            
        return jobs

    def diamond_blastp_nr(self):
        
        """
        Step diamond_blastp_nr(): DIAMOND BLASTp of predicted genes vs NCBI nr
        """
        
        jobs = []
        
        for bin in self.magsets:   
            #log.info("MAG name: " + bin.name)
            #log.info("MAG file: " + bin.fasta_file)
            if not os.path.exists(os.path.join(bin.name)):
                os.makedirs(os.path.join(bin.name))
            if not os.path.exists(os.path.join(bin.name, "annotations", "blastp_nr")):
                os.makedirs(os.path.join(bin.name, "annotations", "blastp_nr"))
       
            fname = os.path.join(bin.name, "annotations", "blastp_nr_annotated.tsv")
            open(fname, 'a').close()
            os.utime(fname, None)
            
            # Do blastp on KEGG for big assembly  
            chunks_dir = os.path.join(bin.name, "gene_prediction", "fasta_chunks")
            blast_dir = os.path.join(bin.name, "annotations", "blastp_nr")
            number_chunks_file = os.path.join(bin.name, "gene_prediction", "estimated_number_of_chunks_genes.txt")
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
            job.name = bin.name + "_" + "diamond_blastp_array_nr"
            job.subname = "blastp_nr"
            job.job_array_num_task = num_chunks
            jobs.append(job)
            
            # Merge output chunks
            job = shotgun_metagenomics.merge_chunks(
                blast_dir,
                os.path.join(bin.name, "annotations", "blastp_nr.tsv"),
                num_chunks,
                "blastp"
            )
            job.name = bin.name + "_" + "merge_diamond_blastp_nr"
            job.subname = "merge"
            jobs.append(job)
        
        return jobs 

    def taxonomic_annotation(self):
        """
        Step taxonomic_annotation(): Using DIAMOND blastp results generated in the diamond_blastp vs nr step, each contig
                                     will be assigned a taxonimic lineage using the CAT package. Taxonomic summaries and 
                                     figures are then generated.
        """
        jobs = []
        for bin in self.magsets:   
            
            #log.info("MAG name: " + bin.name)
            #log.info("MAG file: " + bin.fasta_file)
            if not os.path.exists(os.path.join(bin.name)):
                os.makedirs(os.path.join(bin.name))
            if not os.path.exists(os.path.join(bin.name, "annotations", "taxonomy_consensus")):
                os.makedirs(os.path.join(bin.name, "annotations", "taxonomy_consensus"))
        
            prefixes = ["consensus"]
            
            job = shotgun_metagenomics.convert_orf_ids_for_cat(
               os.path.join(bin.name, "gene_prediction", "Contigs_renamed.gff"),
               os.path.join(bin.name, "annotations", "blastp_nr.tsv"),
               os.path.join(bin.name, "annotations", "blastp_nr_alt_orf_ids.tsv")
            )
            job.name = bin.name + "_" + "convert_ord_ids_for_CAT"
            job.subname = "convert_ord_ids_for_CAT"
            jobs.append(job)
            
            job = shotgun_metagenomics.CAT(
               os.path.join(bin.name, "assembly", "Contigs.fasta"),
               os.path.join(bin.name, "gene_prediction", "Contigs.faa"),
               os.path.join(bin.name, "annotations", "blastp_nr_alt_orf_ids.tsv"),
               os.path.join(bin.name, "annotations", "taxonomy_consensus", "out")
            )
            job.name = bin.name + "_" + "CAT"
            job.subname = "CAT"
            jobs.append(job)
        
            job = genome_assembly.generate_taxonomy_table_consensus(
               os.path.join(bin.name, "annotations", "taxonomy_consensus", "out.contig2classification_with_names.tsv"),
               os.path.join(bin.name, "annotations", "taxonomy_consensus", "taxonomy.tsv")
            )
            job.name = "generate_taxonomy_table_consensus"
            job.subname = "generate_taxonomy_table"
            jobs.append(job)
            
        return jobs


    # Here generate final GFF (for viewing data in a genome browser). 
    # And generate final DDA sheets. with logFC, gene_name and actual normalized values.
    def finalize(self):
        jobs = []
        
        for bin in self.magsets:   
            
            #log.info("MAG name: " + bin.name)
            #log.info("MAG file: " + bin.fasta_file)
            if not os.path.exists(os.path.join(bin.name)):
                os.makedirs(os.path.join(bin.name))
        
            job = shotgun_metagenomics.generate_gff(
                # infiles
                os.path.join(bin.name, "gene_prediction", "Contigs_renamed.gff"),
                os.path.join(bin.name, "assembly", "Contigs.fasta"),
                os.path.join(bin.name, "annotations", "blastp_kegg_parsed.tsv"),
                os.path.join(bin.name, "annotations", "hmmscan_pfam_tblout.tsv"),
                os.path.join(bin.name, "annotations", "rpsblast_cog.tsv"),
                os.path.join(bin.name, "annotations", "rpsblast_kog.tsv"),
                os.path.join(bin.name, "annotations", "taxonomy_consensus", "taxonomy.tsv"),
                os.path.join(bin.name, "annotations", "blastp_nr_annotated.tsv"),
                # outfiles
                os.path.join(bin.name, "annotations", "Contigs.gff"),
                os.path.join(bin.name, "annotations", "Contigs.fasta"),
                os.path.join(bin.name, "annotations", "annotations.tsv")
            )
            job.name = bin.name + "_" + "generate_gff"
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
    

    @property
    def steps(self):
        
        return [
            # Core steps.
            self.gene_prediction,
            self.exonerate,
            self.diamond_blastp_kegg,
            self.rpsblast_cog,
            self.rpsblast_kog,
            self.hmmscan_pfam,
            self.diamond_blastp_nr,
            self.ncrna,
            self.taxonomic_annotation,
            self.finalize
        ]

    def set_local_variables(self):
        self._parser_local = self.argparser

        
        # barcodes
        self._args_local = self._parser_local.parse_args()
        config.parse_files(self._args_local.config)
        self._config = self._args_local.config[0]
        self._config = os.path.abspath(self._config.name)
        self._curr_scheduler = self._args_local.job_scheduler
        
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
#################################################################################################
   _____                                                            _        _   _
  / ____|                                   /\                     | |      | | (_)
 | |  __  ___ _ __   ___  _ __ ___   ___   /  \   _ __  _ __   ___ | |_ __ _| |_ _  ___  _ __
 | | |_ |/ _ \ '_ \ / _ \| '_ ` _ \ / _ \ / /\ \ | '_ \| '_ \ / _ \| __/ _` | __| |/ _ \| '_ \\
 | |__| |  __/ | | | (_) | | | | | |  __// ____ \| | | | | | | (_) | || (_| | |_| | (_) | | | |
  \_____|\___|_| |_|\___/|_| |_| |_|\___/_/    \_\_| |_|_| |_|\___/ \__\__,_|\__|_|\___/|_| |_|

               Support: jtremblay514@gmail.com
             Home page: jtremblay.github.io/pipelines.html

#################################################################################################"""
        sys.stderr.write(metagenomics_string + '\n')
        time.sleep(1)
        # Add pipeline specific arguments
        self.argparser.add_argument("-m", "--magsets", help="magset file", type=argparse.FileType('r'), required=True)
        #self.argparser.add_argument("-r", "--readsets", help="readset file", type=file, required=True)
        self.set_local_variables()
        sys.stderr.write('Genome Annotation pipeline\n')
        super(GenomeAnnotater, self).__init__()
                
GenomeAnnotater().submit_jobs()
