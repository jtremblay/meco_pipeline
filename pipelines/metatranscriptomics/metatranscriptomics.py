#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys
import errno

# Append caf_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# CAF Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.design import *
from bio.readset import *
from bio.utils import *

from bio import shotgun_metagenomics
from bio import shotgun_metatranscriptomics
from bio import microbial_ecology
from bio import rnaseq
from bio import utils

from pipelines.illumina import illumina

# Global scope variables
log = logging.getLogger(__name__)

class Metatranscriptomics(illumina.Illumina):

    def trim(self):
        jobs = []
        # Merge all demultiplexed fastq files in one file. One file for reads1 and one file for reads2 of Illumina paired-end.
        # If library is Illumina single end or 454 or PacBio, only one file will be generated.
        outdir = self._root_dir
        trim_stats = []
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

                trim_stats.append(trim_file_prefix + "stats.csv")

            elif readset.run_type == "SINGLE_END":
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (SINGLE_END not implemented yet)!")
            
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
            job.name = "trimmomatic." + readset.name
            job.subname = "trim"
            jobs.append(job)
        
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
        
        return jobs
            
    def duk(self):
        jobs=[]
        outdir = self._root_dir
        logs = []
        readset_ids = []

        for readset in self.readsets:
            trim_file_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".trim.")
            outfile_prefix = os.path.join(outdir, "qced_reads", readset.sample.name, readset.name + ".")
            out1up = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_unpaired_R1.fastq.gz")
            out2up = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_unpaired_R2.fastq.gz")
            out1p = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
            out2p = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
            log = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".duk_contam_interleaved_log.txt")

            job = shotgun_metagenomics.duk_gz(
                trim_file_prefix + "interleaved.fastq.gz",
                outfile_prefix + "contam.fastq.gz",
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
     
    def abundance(self):
        jobs = []

        cov_list_contigs = []
        cov_list_genes = []
        flagstats_contigs_list = []
        trimmomatic_list = []
        #htseq_list_contigs = []
        
        # Will make index for bwa. and also bed file for computing reads spanning later.
        reference_contigs = os.path.join(self._root_dir, "assembly", "Contigs.fasta")
        bed_contigs =       os.path.join(self._root_dir, "assembly", "Contigs.bed")
        bwt_contigs =       os.path.join(self._root_dir, "assembly", "Contigs.fasta.bwt")
        bed_genes =         os.path.join(self._root_dir, "assembly", "Contigs_genes.bed")
        
        # Will make index for bwa. and also bed file for computing reads spanning later.
        for readset in self.readsets:
            # Trimmomatic reads.
            trimmomatic = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".trim.stats.csv")
            trimmomatic_list.append(trimmomatic)
            # BAM contigs
            bam_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".bam")
            bam_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".bam")
            flagstats_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".flagstats")
            flagstats_contigs_list.append(flagstats_contigs)
            # Cov contigs
            cov_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".cov")
            cov_list_contigs.append(cov_contigs)
            # Cov genes
            cov_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".cov")
            cov_list_genes.append(cov_genes)
            # Outdir
            outdir_genes = os.path.join(self._root_dir, "gene_abundance", readset.sample.name)
            outdir_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name)

            if not os.path.exists(outdir_genes):
                os.makedirs(os.path.join(outdir_genes))
            if not os.path.exists(outdir_contigs):
                os.makedirs(os.path.join(outdir_contigs))
       
            if( readset.run_type == "PAIRED_END"):
                infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
                infile_R1 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
                infile_R2 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
                flag = "0x2"
                
                # map against contigs
                if config.param('star_align', 'do_star', 1, 'string') == "yes":
                    job = rnaseq.star_align(
                        infile_R1,
                        infile_R2,
                        outdir_genes,
                        config.param('star_align', 'genome_index_folder', 1, 'string'),
                        bam_genes
                    )
                    job.name = "star_align" + readset.sample.name
                    job.subname = "star_align"
                    jobs.append(job)

                    curr_bam = bam_genes

                else:
                    job = shotgun_metagenomics.bwa_mem_samtools(
                        reference_contigs,
                        bwt_contigs,
                        infile_R1,
                        infile_R2,
                        bam_contigs
                    )
                    job.name = "bwa_mem-contigs" + readset.sample.name
                    job.subname = "bwa"
                    jobs.append(job)
                    
                    curr_bam = bam_contigs
                
                job = shotgun_metagenomics.flagstats(
                    curr_bam,
                    flagstats_contigs
                )
                job.name = "flagstats-contigs-" + readset.sample.name
                job.subname = "flagstats"
                jobs.append(job)
                
            #elif(readset.run_type == "SINGLE_END"):
            #    infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_R1.fastq.gz")
            #    
            #    job = shotgun_metagenomics.bwa_mem_se(
            #        reference,
            #        infile,
            #        bam
            #    )
            #    job.name = "bwa_mem-" + readset.sample.name
            #    job.subname = "bwa_short"
            #    jobs.append(job)

            if config.param('star_align', 'do_star', 1, 'string') == "yes":
                # raw counts with htseq. For Htseq, use contigs abundance AND gene models (GFF).
                job = shotgun_metatranscriptomics.htseq_count(
                    bam_genes,
                    config.param('star_align', 'gff', 1, 'filepath'),
                    cov_genes
                )
                job.name = "htseq-genes" + readset.sample.name
                job.subname = "htseq"
                jobs.append(job)

            else:
                # raw counts with bedtools. Bedtools used for genes abundance only
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
            
            ## raw counts with htseq. For Htseq, use contigs abundance AND gene models (GFF).
            #job = shotgun_metatranscriptomics.htseq_count(
            #    bam_contigs,
            #    config.param('DEFAULT', 'gff', 1, 'filepath'),
            #    htseq_contigs
            #)
            #job.name = "htseq-contigs" + readset.sample.name
            #job.subname = "htseq"
            #jobs.append(job)

        # merge alignment rates for contigs
        job = shotgun_metagenomics.merge_flagstats(
            flagstats_contigs_list,
            trimmomatic_list,
            os.path.join(self._root_dir, "contigs_abundance", "qc_mapping_stats_contigs.tsv")
        )
        job.name = "flagstats_merge_contigs"
        job.subname = "flagstats"
        jobs.append(job)
        
        if config.param('star_align', 'do_star', 1, 'string') == "yes":
            #sys.stderr.write("Not implemented yet...\n")
            job = shotgun_metatranscriptomics.merge_htseq_counts(
                cov_list_genes,
                os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
                os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance_cpm.tsv")
            )
            job.name = "merge_gene_abundance_genes"
            job.subname = "merge_gene_abundance"
            jobs.append(job)

        else:
            # Count reads with bedtools and htseq
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
        
        return jobs 
    
    def DEG(self):
        
        jobs = []

        #job = shotgun_metagenomics.draw_pcoa(
        #    os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance_cpm.tsv"),
        #    config.param('DEFAULT', 'mappingFile'),
        #    os.path.join(self._root_dir, "gene_abundance", "plots")
        #)
        #job.name = "drawPcoa_bedtools"
        #job.subname = "drawPcoa"
        #jobs.append(job)
        
        # edgeR for DEG - bedtools
        if(config.param('DEG', 'do_deg_pairwise', 1, 'string') == 'yes'):
            design_files = config.param('DEFAULT', 'designFile', 1, 'string')
            designs = design_files.split(":")
            for design_file in designs:
                job = shotgun_metatranscriptomics.edger(
                    os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
                    os.path.join(self._root_dir, design_file),
                    os.path.join(self._root_dir, "DEG_bedtools")
                )
                job.name = "DEG_bedtools"
                job.subname = "DEG"
                jobs.append(job)
         
        # Do not necessarily do GLMs.
        if config.param('DEG', 'do_deg_glm', 1, 'string') == 'yes':
            mapping_files = config.param('DEG', 'mapping_files', 1, 'string')
            mappings = mapping_files.split(":")
            for map in mappings:
                job = shotgun_metatranscriptomics.edger_glm(
                    os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
                    map,
                    os.path.join(self._root_dir, "DEG_GLM")
                )
                job.name = "DEG_GLM"
                job.subname = "DEG"
                jobs.append(job)

        # if rnaseq, 
        if config.param('star_align', 'do_star', 1, 'string') == "yes":
            #sys.stderr.write("Not implemented yet...\n")
            job = shotgun_rnaseq.annotate_edger(
                cov_list_genes,
                os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
                os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance_cpm.tsv")
            )
            job.name = "merge_gene_abundance_genes"
            job.subname = "merge_gene_abundance"
            jobs.append(job)
        # For metagenomic binning centric analyses.
        #if config.param('DEG', 'do_deg_bins', 1, 'string') == 'yes':
        #    job = shotgun_metagenomics.deg_bins(
        #        os.path.join(self._root_dir, "DEG_bedtools", "DEG_normalized_significant.tsv"),
        #        os.path.join(self._root_dir, "assembly", "Contigs_genes.bed"),
        #        #os.path.join(self._root_dir, "DEG_bedtools")
        #    )
        #    job.name = "DEG_GLM"
        #    job.subname = "DEG"
        #    jobs.append(job)   
    
        return jobs 
        
    def betadiv(self):
        jobs = []
        # Beta diversity.
        def mkdir_p(path):
            try:
                os.makedirs(path)
            except OSError as exc: # Python >2.5
                if exc.errno == errno.EEXIST and os.path.isdir(path):
                    pass
                else: raise
        mkdir_p(os.path.join(self._root_dir, "betadiv"))
        fname = os.path.join(self._root_dir, "betadiv", "tree.fasttree")
        open(fname, 'a').close()
        os.utime(fname, None)

        job = microbial_ecology.beta_diversity(
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "otu_table_normalized.biom"),
            "bray_curtis",
            self._root_dir + "betadiv/",
            self._root_dir + "betadiv/tree.fasttree",
            "bray_curtis_otu_table_normalized"
        )
        job.name = "beta_diversity_bray_curtis_all"
        job.subname = "beta_diversity"
        jobs.append(job)
        
        job = microbial_ecology.beta_diversity(
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "otu_table_normalized_bacteriaArchaea.biom"),
            "bray_curtis",
            self._root_dir + "betadiv/",
            self._root_dir + "betadiv/tree.fasttree",
            "bray_curtis_otu_table_normalized_bacteriaArchaea"
        )
        job.name = "beta_diversity_bray_curtis_bactarch"
        job.subname = "beta_diversity"
        jobs.append(job)
        
        job = microbial_ecology.principal_coordinates(
            self._root_dir + "betadiv/bray_curtis_" + "otu_table_normalized" + ".txt",
            self._root_dir + "betadiv/bray_curtis_" + "otu_table_normalized" + "_coords.tsv"
        )
        job.name = "beta_diversity_bray_curtis_PC_all"
        job.subname = "pc"
        jobs.append(job)
        
        job = microbial_ecology.principal_coordinates(
            self._root_dir + "betadiv/bray_curtis_" + "otu_table_normalized_bacteriaArchaea" + ".txt",
            self._root_dir + "betadiv/bray_curtis_" + "otu_table_normalized_bacteriaArchaea" + "_coords.tsv"
        )
        job.name = "beta_diversity_bray_curtis_PC_bactarch"
        job.subname = "pc"
        jobs.append(job)

        return jobs 
     
    def taxonomy(self):
        jobs = []
            
        # We already have taxonomy from metagenomics (assembly) pipelne.
        job = shotgun_metagenomics.generate_otu_table(
            os.path.join(self._root_dir, "assembly", "taxonomy.tsv"),
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "otu_table.tsv")
        )
        job.name = "generate_otu_table"
        job.subname = "generate_otu_table"
        jobs.append(job)
        
        prefixes = [
            "genes"
        ]
        infiles = [
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes", "otu_table.tsv")
        ]
        directories = [
            os.path.join(self._root_dir, "gene_annotation", "taxonomy_genes")
        ]
        abundances = [
            os.path.join(self._root_dir, "gene_abundance", "merged_gene_abundance_cpm.tsv")
        ] 
        
        for i in range(0, len(infiles)):
            job = microbial_ecology.convert_otu_to_biom_hdf5(
                infiles[i],
                os.path.join(directories[i], "otu_table.biom")
            )
            job.name = "convert_otu_to_biom_" + prefixes[i]
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
                "bactArch"
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
   
           
            normalizations = ["normalized"]
            types = ["absolute", "relative"]
            organisms = ["", "_others", "_bacteriaArchaea"]
            organisms2 = ["all", "others", "bacteriaArchaea"]

            for n in range(0, len(normalizations)):
                for j in range(0, len(types)):
                    for k in range(0, len(organisms)):
                        for m in range(1, config.param('summarize_taxonomy', 'taxonomyDepth', 1, 'int')):
                            if types[j] == "absolute":
                                job = microbial_ecology.summarize_taxonomy_absolute(
                                    os.path.join(directories[i], "otu_table_" + normalizations[n] + organisms[k] + ".biom"),
                                    m,
                                    os.path.join(directories[i], organisms2[k], types[j]) 
                                )
                                job.name = "tax_summary_" + normalizations[n] + "_" + types[j] + "_" + organisms2[k]  + "_L" + str(m) + "_" + prefixes[i]
                                job.subname = "summarize_taxonomy"
                                jobs.append(job)
                            elif types[j] == "relative":
                                job = microbial_ecology.summarize_taxonomy_relative(
                                    os.path.join(directories[i], "otu_table_" + normalizations[n] + organisms[k] + ".biom"),
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
                                os.path.join(directories[i], organisms2[k], types[j], "otu_table_" + normalizations[n] + organisms[k] + "_L" + str(m) + ".txt"),
                                os.path.join(directories[i], organisms2[k], types[j], "plots"),
                                "taxonomy_L" + str(m) + "_" + types[j]
                            )
                            job.name = "plot_taxa_single_" + normalizations[n] + "_" + types[j]  + "_L" + str(m) + "_" + prefixes[i] + "_" + organisms2[k]
                            job.subname = "plot_taxa"
                            jobs.append(job)
                    
        return jobs

    def overrepresentation (self):
        jobs = []

        types = ["pathways", "modules", "K"]
        for type in types:
            job = shotgun_metagenomics.kegg_overrep(
                os.path.join(self._root_dir, "DEG_bedtools"),
                os.path.join(self._root_dir, "gene_annotation", "blastp_kegg_parsed.tsv"),
                type,
                os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
                "DEG"
            )
            job.name = "getKegg_" + type
            job.subname = "getKegg"
            jobs.append(job)
         
        # COG and KOG
        job = shotgun_metagenomics.cog_overrep(
            os.path.join(self._root_dir, "DEG_bedtools"),
            os.path.join(self._root_dir, "gene_annotation", "rpsblast_cog.tsv"),
            os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
            "DEG"
        )
        job.name = "getCOG"
        job.subname = "getCOG"
        jobs.append(job)
        
        job = shotgun_metagenomics.kog_overrep(
            os.path.join(self._root_dir, "DEG_bedtools"),
            os.path.join(self._root_dir, "gene_annotation", "rpsblast_kog.tsv"),
            os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
            "DEG"
        )
        job.name = "getKOG"
        job.subname = "getKOG"
        jobs.append(job)
        
        # Pfam and tigrfam
#        job = shotgun_metagenomics.pfam_overrep(
#            os.path.join(self._root_dir, "DEG_bedtools"),
#            os.path.join(self._root_dir, "gene_annotation", "hmmscan_pfam_tblout.tsv"),
#            os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
#            "DEG"
#        )
#        job.name = "getPfam"
#        job.subname = "getPfam"
#        jobs.append(job)
#         
#        job = shotgun_metagenomics.tigrfam_overrep(
#            os.path.join(self._root_dir, "DEG_bedtools"),
#            os.path.join(self._root_dir, "gene_annotation", "hmmscan_tigrfam_tblout.tsv"),
#            os.path.join(self._root_dir, "DEG_bedtools", "DEG.done"),
#            "DEG"
#        )
#        job.name = "getTigrfam"
#        job.subname = "getTigrfam"
#        jobs.append(job)
#        
#        # heatmaps
#        prefixes = ["DEG_normalized_significant_keggK.tsv", "DEG_normalized_significant_keggModules.tsv", "DEG_normalized_significant_keggPathways.tsv", "DEG_normalized_significant_COG.tsv", "DEG_normalized_significant_KOG.tsv", "DEG_normalized_significant_PFAM.tsv", "DEG_normalized_significant_TIGRFAM.tsv"]
#        for prefix in prefixes:
#            job = shotgun_metagenomics.draw_heatmap(
#                os.path.join(self._root_dir, "DEG_bedtools"),
#                config.param('DEFAULT', 'mappingFile'),
#                prefix,
#                os.path.join(self._root_dir, "DEG_bedtools", "DEG.done")
#            )
#            job.name = "draw_heatmap"
#            job.subname = "draw_heatmap_" +  prefix
#            jobs.append(job)
         
        return jobs

    # Here generate final GFF (for viewing data in a genome browser). 
    # And generate final DEG sheets. with logFC, gene_name and actual normalized values.
    def finalize(self):
        jobs = []
        
        DEGs = ["DEG","DEG_bedtools"]

        for DEG in DEGs:
            job = shotgun_metagenomics.annotate_DDA(
                os.path.join(self._root_dir, DEG),
                config.param('DEFAULT', 'annotations', 1, 'filepath')
            )
            job.name = "annotate_" + DEG
            job.subname = "annotate_DEG"
            jobs.append(job)

        return jobs

    def generate_report(self):
        jobs = []
        
        i = datetime.datetime.now()
        year = "%04d" %i.year
        day = "%02d" %i.day
        month = "%02d" %i.month
        title = config.param('report', 'report.title')
        date = year + "_" + month + "_" + day + "_metatranscriptomics_" + title
        
        job = shotgun_metagenomics.concatenate_reports(
            self._config,
            self._root_dir,
            "Metatranscriptomics",
            self._root_dir + "/" + date
        )
        job.name = "deliverables"
        job.subname = "report"
        jobs.append(job)
        return jobs
    
    def exploratory(self):
        jobs = []
        
        i = datetime.datetime.now()
        year = "%04d" %i.year
        day = "%02d" %i.day
        month = "%02d" %i.month
        title = config.param('report', 'report.title')
        date = year + "_" + month + "_" + day

        log_fc_string = config.param('exploratory', 'cutoff_logFC', 1, 'string')
        log_fcs = log_fc_string.split(":")
                
        gene_list_files = config.param('exploratory', 'gene_list', 1, 'string')
        gene_lists = gene_list_files.split(":")
        
        for gene_list in gene_lists:
            gene_list_prefix = os.path.splitext(gene_list)[0]
            gene_list_prefix = gene_list_prefix.replace("./","")
            for log_fc in log_fcs:
                if not os.path.exists( os.path.join(self._root_dir, "DEG_bedtools", "exploratory", "log_fc-"+log_fc) ):
                    os.makedirs( os.path.join(self._root_dir, "DEG_bedtools", "exploratory", "log_fc-"+log_fc) )
                
                job = shotgun_metatranscriptomics.explore_heatmaps(
                    log_fc,
                    os.path.join(self._root_dir, "DEG_bedtools", "DEG_normalized_significant.tsv"),
                    os.path.abspath(os.path.join(self._root_dir, "DEG_bedtools")),
                    os.path.abspath(os.path.join(self._root_dir, "DEG_bedtools", "exploratory", "log_fc-" + log_fc)),
                    date,
                    gene_list
                )
                job.name = "exploratory_" + gene_list_prefix + "_" + log_fc
                job.subname = "exploratory"
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
            self.trim,#1
            self.duk,#2
            self.abundance,#3
            self.DEG,#4
            self.overrepresentation,#5
            self.finalize,#6
            self.taxonomy,#7
            self.generate_report,#8
            self.exploratory,#9
            self.cleanup
        ]

    def set_local_variables(self):
        self._parser_local = argparse.ArgumentParser(description='Process options.')
        self._parser_local.add_argument("-c", "--config", help="config INI-style file", nargs="+", type=file, required=True)
        self._parser_local.add_argument("-r", "--readsets", help="readset file", type=file, required=False)
        self._parser_local.add_argument("-s", "--steps", help="step range e.g. '1-5', '3,6,7', '2,4-8'", required=True)
        self._parser_local.add_argument("-o", "--output-dir", help="output directory (default: current)", default=os.getcwd())
        self._parser_local.add_argument("-j", "--job-scheduler", help="job scheduler type (default: slurm)", choices=["torque", "batch", "daemon", "slurm"], default="slurm")
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
        if not os.path.exists(root_dir + "/DEG"):
            os.makedirs(os.path.join(root_dir, "DEG"))
        if not os.path.exists(root_dir + "/qced_reads"):
            os.makedirs(os.path.join(root_dir, "qced_reads"))
        if not os.path.exists(root_dir + "/abundance"):
            os.makedirs(os.path.join(root_dir, "abundance"))
        if not os.path.exists(root_dir + "/taxonomy"):
            os.makedirs(os.path.join(root_dir, "taxonomy"))
        if not os.path.exists(root_dir + "gene_annotation/taxonomy_genes/all/absolute/plots"):
            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/all/absolute/plots"))
            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/all/relative/plots"))
        if not os.path.exists(root_dir + "/gene_annotation/taxonomy_genes/others/absolute/plots"):
            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/others/absolute/plots"))
            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/others/relative/plots"))
        if not os.path.exists(root_dir + "gene_annotation/taxonomy_genes/bacteriaArchaea/absolute/plots"):
            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/bacteriaArchaea/absolute/plots"))
            os.makedirs(os.path.join(root_dir, "gene_annotation/taxonomy_genes/bacteriaArchaea/relative/plots"))
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
        sys.stderr.write('Metatranscriptomics pipeline\n')
        super(Metatranscriptomics, self).__init__()
                
Metatranscriptomics().submit_jobs()
