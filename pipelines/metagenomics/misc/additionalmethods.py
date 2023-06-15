
    def generate_report(self):
        jobs = []
        
        i = datetime.datetime.now()
        year = "%04d" %i.year
        day = "%02d" %i.day
        month = "%02d" %i.month
        title = config.param('report', 'report.title')
        date = year + "_" + month + "_" + day + "_metagenomics_" + title
        
        job = shotgun_metagenomics.generate_report(
            self._config,
            self._root_dir,
            "Metagenomics",
            self._root_dir + "/" + date
        )
        job.name = "deliverables"
        job.subname = "report"
        jobs.append(job)
        return jobs
    
    # perform classification based on jgi itagger paper approach for MG data.
    def emirge(self):
        jobs = []
        
        emirge_fastas = []
        emirge_rdps = []
        emirge_taxonomies = []
        names = []

        for readset in self.readsets:
            in1p = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq")
            in2p = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq")
            emirge_fasta = os.path.join(self._root_dir, "gene_annotation", "emirge", readset.sample.name, readset.name + ".fasta")
            emirge_rdp = os.path.join(self._root_dir, "gene_annotation", "emirge", readset.sample.name, readset.name + "_rdp.tsv")
            emirge_taxonomy = os.path.join(self._root_dir, "gene_annotation", "emirge", readset.sample.name, readset.name + "_taxonomy.tsv")
            outdir = os.path.join(self._root_dir, "gene_annotation", "emirge", readset.sample.name)
            emirge_fastas.append(emirge_fasta)
            emirge_rdps.append(emirge_rdp)
            emirge_taxonomies.append(emirge_taxonomy)
            names.append(readset.sample.name)

            # For each of these potential rRNA reads, perform rdp classification
            job = shotgun_metagenomics.emirge(
                outdir,
                in1p,
                in2p,
                readset.sample.name,
                os.path.join(self._root_dir, "contigs_abundance", "lib_stats.tsv")
            )
            job.name = "emirge_" + readset.sample.name
            job.subname = "emirge"
            jobs.append(job)
       
            # BLASTN NCBI 
            job = shotgun_metagenomics.blastn(
                emirge_fasta,
                os.path.join(self._root_dir, "gene_annotation", "emirge", readset.sample.name, readset.name + "_blastn_nt.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "emirge", readset.sample.name, readset.name),
                "blastn"
            )
            job.name = "blastn_nt_emirge_" + readset.sample.name
            job.subname = "blastn"
            jobs.append(job)

            # extract BLASTN Taxonomy
            job = shotgun_metagenomics.extract_taxonomy(
                os.path.join(self._root_dir, "gene_annotation", "emirge", readset.sample.name, readset.name + "_blastn_nt.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "emirge", readset.sample.name, readset.name + "_taxonomy.tsv")
            )
            job.name = "extract_blastn_taxonomy_emirge_" + readset.sample.name
            job.subname = "ncbi_tax"
            jobs.append(job)
            
            # For each of these potential rRNA reads, perform rdp classification
            job = microbial_ecology.rdp_wrapper(
                emirge_fasta,
                os.path.join(self._root_dir, "gene_annotation", "emirge", readset.sample.name, readset.name + "_rdp.tsv")
            )
            job.name = "classify_emirge" + readset.sample.name
            job.subname = "RDP"
            jobs.append(job)
            
        # with rdp results in hand, generate otu table
        job = shotgun_metagenomics.emirge_to_otu_table(
            emirge_fastas, #fastas contains abundance values here...
            emirge_taxonomies,
            os.path.join(self._root_dir, "gene_annotation", "emirge", "otu_table_blastn.tsv"),
            names,
            "blastn"
        )
        job.name = "emirge_otu_table_blastn"
        job.subname = "emirge_otu_table"
        jobs.append(job)
        
        # with rdp results in hand, generate otu table
        job = shotgun_metagenomics.emirge_to_otu_table(
            emirge_fastas, #fastas contains abundance values here...
            emirge_rdps,
            os.path.join(self._root_dir, "gene_annotation", "emirge", "otu_table_rdp.tsv"),
            names,
            "rdp"
        )
        job.name = "emirge_otu_table_rdp"
        job.subname = "emirge_otu_table"
        jobs.append(job)

        return jobs


    def plasmids(self):
        jobs = []
        
        # Do blastn on nt for big assembly 
        chunks_dir = os.path.join(self._root_dir, "assembly", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "plasmids")
        num_chunks = config.param('exonerate_contigs', 'num_fasta_chunks', type='posint')

        for i in range(num_chunks):
            job = shotgun_metagenomics.blastn(
                os.path.join(chunks_dir, "Contigs.fasta_chunk_{:07d}".format(i)),
                os.path.join(blast_dir, "blastn_chunk_{:07d}.tsv".format(i)),
                blast_dir,
                "blastn_plasmids"
            )
            job.name = "blastn_plasmids"
            job.subname = "blastn_plasmids"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_chunks(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "blastn_plasmids.tsv"),
            num_chunks,
            "blastn" 
        )
        job.name = "blastn_plasmids_merge"
        job.subname = "blastn_merge"
        jobs.append(job)
        
        job = shotgun_metagenomics.keep_blast_best_hit(
            os.path.join(self._root_dir, "gene_annotation", "blastn_plasmids.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "blastn_plasmids_besthit.tsv")
        )
        job.name = "blastn_best_hit_plasmids"
        job.subname = "blastn_best_hit"
        jobs.append(job)
        
        job = shotgun_metagenomics.extract_plasmids(
            os.path.join(self._root_dir, "gene_annotation", "blastn_plasmids_besthit.tsv"),
            os.path.join(self._root_dir, "assembly", "Contigs.fasta"),
            os.path.join(self._root_dir, "gene_prediction", "Contigs_renamed.faa"),
            os.path.join(self._root_dir, "gene_annotation", "hmmscan_pfam_tblout.tsv"),
            os.path.join(self._root_dir, "gene_prediction", "Contigs_renamed.gff"),
            os.path.join(self._root_dir, "gene_annotation", "plasmids.fna"),
            os.path.join(self._root_dir, "gene_annotation", "plasmids.faa"),
            os.path.join(self._root_dir, "gene_annotation", "plasmids.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "plasmids_gene_list.tsv")
        )
        job.name = "extract_plasmids"
        job.subname = "extract_plasmids"
        jobs.append(job)
        
        return jobs
#    def hmmscan_tigrfam(self):
#        jobs = [] # Do rpsblast on COG for big assembly
# 
#        chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
#        blast_dir = os.path.join(self._root_dir, "gene_annotation", "hmmscan_tigrfam")
#        num_chunks = config.param('exonerate', 'num_fasta_chunks', type='posint', required=True)
#
#        for i in range(num_chunks):
#            job = shotgun_metagenomics.hmmscan(
#                os.path.join(chunks_dir, "Contigs_renamed.faa_chunk_{:07d}".format(i)),
#                os.path.join(blast_dir, "hmmscan_chunk_{:07d}".format(i)),
#                blast_dir,
#                config.param('tigrfam', 'db', required=True) 
#            )
#            job.name = "hmmscan_tigrfam"
#            job.subname = "hmmscan"
#            jobs.append(job)
#      
#        # Merge output chunks
#        job = shotgun_metagenomics.merge_chunks_hmms(
#            blast_dir,
#            os.path.join(self._root_dir, "gene_annotation"),
#            num_chunks,
#            "hmmscan",
#            "hmmscan_tigrfam"
#        )
#        job.name = "hmmscan_tigrfam_merge"
#        job.subname = "merge"      
#        jobs.append(job)
#
#        return jobs 


    def ublastp_kegg(self):

        jobs = []

        # Do ublast on nr 
        #chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
        blast_dir = os.path.join(self._root_dir, "gene_annotation", "ublastp_kegg")
        num_chunks = config.param('ublastp', 'num_chunks_db_kegg', 1, 'posint')

        # Here loop through all 30 nr_\d+.fasta dbs instead of looping through queries.
        for i in range(1,num_chunks):
            job = shotgun_metagenomics.ublastp(
                os.path.join(self._root_dir, "gene_prediction", "Contigs_renamed.faa".format(i)),
                os.path.join(blast_dir, "ublastp_chunk_{:03d}.tsv".format(i)),
                blast_dir,
                os.path.join(config.param('ublastp', 'db_kegg', 1, 'string') + "_{:03d}.udb".format(i))
            )
            job.name = "ublastp_kegg"
            job.subname = "ublastp"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.keep_ublast_best_hit(
            blast_dir,
            os.path.join(self._root_dir, "gene_annotation", "ublastp_kegg.tsv"),
            num_chunks,
            "ublastp"
        )
        job.name = "ublastp_kegg_keep_best_hit"
        job.subname = "keep_best_hit"
        jobs.append(job)
        
        # Generate a clean table of module/pathways.
        job = shotgun_metagenomics.parse_kegg(
            os.path.join(self._root_dir, "gene_annotation", "ublastp_kegg.tsv"),
            os.path.join(self._root_dir, "gene_annotation", "blastp_kegg_parsed.tsv") # Should be ublastp_kegg_parsed.tsv, but to facilitate with compatibility
        )
        job.name = "parse_kegg"
        job.subname = "merge"
        jobs.append(job)
        
        return jobs 
    
    def ublastp_nr(self):

        jobs = []
        if(config.param('ublastp', 'skip_ublastp_nr', 1, 'string') == 'yes'):
           sys.stderr.write("Skipping ublastp_nr step...\n")
           fname = os.path.join(self._root_dir, "gene_annotation", "ublastp_nr_annotated.tsv")
           open(fname, 'a').close()
           os.utime(fname, None)
        else:
            # Do ublast on nr 
            #chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
            blast_dir = os.path.join(self._root_dir, "gene_annotation", "ublastp_nr")
            num_chunks = config.param('ublastp', 'num_chunks_db_nr', 1, 'posint')

            # Here loop through all 30 nr_\d+.fasta dbs instead of looping through queries.
            for i in range(1,num_chunks):
                job = shotgun_metagenomics.ublastp(
                    os.path.join(self._root_dir, "gene_prediction", "Contigs_renamed.faa".format(i)),
                    os.path.join(blast_dir, "ublastp_chunk_{:03d}.tsv".format(i)),
                    blast_dir,
                    os.path.join(config.param('ublastp', 'db_nr', 1, 'string') + "_{:03d}.udb".format(i))
                )
                job.name = "ublastp_nr"
                job.subname = "ublastp"
                jobs.append(job)
          
            # Merge output chunks
            job = shotgun_metagenomics.keep_ublast_best_hit(
                blast_dir,
                os.path.join(self._root_dir, "gene_annotation", "ublastp_nr.tsv"),
                num_chunks,
                "ublastp"
            )
            job.name = "ublastp_nr_keep_best_hit"
            job.subname = "ublastp_keep_best_hit"
            jobs.append(job)
            
            # Annotate ublast output.
            job = shotgun_metagenomics.annotate_ublast(
                os.path.join(self._root_dir, "gene_annotation", "ublastp_nr.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "ublastp_nr_annotated.tsv"),
                annotations =  config.param('ublastp', 'annotations_nr', 1, 'filepath')
            )
            job.name = "annotate_ublast"
            job.subname = "annotate_ublast"
            jobs.append(job)
        
        return jobs 
    
    def kmergenie(self):
        jobs = []
        root = self._root_dir
        
        for readset in self.readsets:
            outfile = os.path.join(root, "assembly", "kmergenie", readset.name + ".kmergenie.txt")
            
            if(readset.run_type == "SINGLE_END"):
                infile = os.path.join(root, "qced_reads", readset.sample.name, readset.name + ".ncontam_R1.fastq.gz")

            elif(readset.run_type == "PAIRED_END"):
                infile = os.path.join(root, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
             
            job = shotgun_metagenomics.kmergenie(
                infile,
                outfile
            )
            job.name = "kmergenie_" + readset.sample.name
            job.subname = "kmergenie"
            jobs.append(job)
        
        return jobs


    def rescaffold(self):
        jobs = []
        
        for readset in self.readsets:
            infile = root + "assembly/" + readset.sample.name + "/Scaffolds.fasta"
            
            job = shotgun_metagenomics.sspace(
                infile,
                outdir
            )
            job.name = "ray_" + readset.sample.name
            job.subname = "ray"
            jobs.append(job)
        
        return jobs
    
    def groopm(self):
        jobs = []
        
        bams = []
        for readset in self.readsets:
            bam_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".bam")
            bams.append(bam_contigs)
        
        job = shotgun_metagenomics.groopm_parse(
            os.path.join(self._root_dir, "assembly", "Contigs.fasta"),
            os.path.join(self._root_dir, "binning", "groopm", "db.gm"),
            os.path.join(self._root_dir, "binning", "groopm", "groopm_parse.done"),
            bams
        )
        job.name = "groopm_parse"
        job.subname = "groopm_parse"
        jobs.append(job)
        
        job = shotgun_metagenomics.groopm_core(
            os.path.join(self._root_dir, "binning", "groopm", "db.gm"),
            os.path.join(self._root_dir, "binning", "groopm", "groopm_parse.done"),
            os.path.join(self._root_dir, "binning", "groopm", "groopm_core.done")
        )
        job.name = "groopm_core"
        job.subname = "groopm_core"
        jobs.append(job)
        
        job = shotgun_metagenomics.groopm_recruit(
            os.path.join(self._root_dir, "binning", "groopm", "db.gm"),
            os.path.join(self._root_dir, "binning", "groopm", "groopm_core.done"),
            os.path.join(self._root_dir, "binning", "groopm", "groopm_recruit.done")
        )
        job.name = "groopm_core"
        job.subname = "groopm_core"
        jobs.append(job)

        return jobs

    def maxbin(self):
        jobs = []
        bams = []
        abundance_files = []

        for readset in self.readsets:
            bam_contigs = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + ".bam")
            bams.append(bam_contigs)
            abundance_file = os.path.join(self._root_dir, "contigs_abundance", readset.sample.name, readset.name + "_abundance.txt")
            abundance_files.append(abundance_file)
        
        job = shotgun_metagenomics.get_abundance_for_maxbin(
            bams,
            os.path.join(self._root_dir, "binning", "maxbin", "coverages.txt"),
            os.path.join(self._root_dir, "binning", "maxbin", "coverages"),
            os.path.join(self._root_dir, "binning", "maxbin", "abundance_fofn.txt")
        )
        job.name = "bamm"
        job.subname = "bamm"
        jobs.append(job)
        
        job = shotgun_metagenomics.maxbin(
            os.path.join(self._root_dir, "assembly", "Contigs.fasta"),
            os.path.join(self._root_dir, "binning", "maxbin", "abundance_fofn.txt"),
            os.path.join(self._root_dir, "binning", "maxbin", "out") 
        )
        job.name = "maxbin"
        job.subname = "maxbin"
        jobs.append(job)

        return jobs


    # Blast contigs against ncbi complete bacterial and archaeal genomes.
#    def ublastn_ncbi_genomes(self):
#
#        jobs = []
#
#        # Do ublastn on ncbi genomes
#        #chunks_dir = os.path.join(self._root_dir, "gene_prediction", "fasta_chunks")
#        blast_dir = os.path.join(self._root_dir, "gene_annotation", "ublastn_genomes")
#        num_chunks = config.param('ublastn', 'num_chunks_db_genomes', 1, 'posint')
#
#        # Here loop through all 30 nt_\d+.fasta dbs instead of looping through queries.
#        for i in range(1,num_chunks):
#            job = shotgun_metagenomics.ublastp(
#                os.path.join(self._root_dir, "assembly", "Contigs.fasta".format(i)),
#                os.path.join(blast_dir, "ublastn_chunk_{:03d}.tsv".format(i)),
#                blast_dir,
#                os.path.join(config.param('ublastn', 'db_genomes', 1, 'string') + "_{:03d}.udb".format(i))
#            )
#            job.name = "ublastn_genomes"
#            job.subname = "ublastn"
#            jobs.append(job)
#      
#        # Merge output chunks
#        job = shotgun_metagenomics.keep_ublast_best_hit(
#            blast_dir,
#            os.path.join(self._root_dir, "gene_annotation", "ublastn_genomes.tsv"),
#            num_chunks,
#            "ublastn"
#        )
#        job.name = "ublastn_genomes_keep_best_hit"
#        job.subname = "ublastn_genomes_keep_best_hit"
#        jobs.append(job)
#        
#        # Annotate ublastn output.
#        #job = shotgun_metagenomics.annotate_ublastn(
#        #    os.path.join(self._root_dir, "gene_annotation", "ublastn_genomes.tsv"),
#        #    os.path.join(self._root_dir, "gene_annotation", "ublastn_genomes_annotated.tsv"),
#        #    annotations =  config.param('ublastn', 'annotations', 1, 'filepath')
#        #)
#        #job.name = "annotate_ublastn_genomes"
#        #job.subname = "annotate_ublastn_genomes"
#        #jobs.append(job)
#        
#        return jobs 


    def rnammer(self):
        jobs = []
        
        chunks_dir = os.path.join(self._root_dir, "assembly", "fasta_chunks")
        rnammer_dir = os.path.join(self._root_dir, "gene_annotation", "rnammer", "chunks")
        num_chunks = config.param('exonerate_contigs', 'num_fasta_chunks', type='posint', required=True)

        for i in range(num_chunks):
            job = shotgun_metagenomics.rnammer(
                os.path.join(chunks_dir, "Contigs.fasta_chunk_{:07d}".format(i)),
                os.path.join(rnammer_dir, "rnammer_chunk_{:07d}.fasta".format(i)),
                rnammer_dir,
                i
            )
            job.name = "rnammer"
            job.subname = "rnammer"
            jobs.append(job)
      
        # Merge output chunks
        job = shotgun_metagenomics.merge_rnammer_chunks(
            rnammer_dir,
            os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer.fasta"),
            num_chunks,
            "rnammer"
        )
        job.name = "rnammer_merge"
        job.subname = "merge"
        jobs.append(job)
        
        # Split rnammer fasta output between subunits.
        job = shotgun_metagenomics.split_rnammer(
            os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_5S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_8S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_16S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_18S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_23S.fasta"),
            os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_28S.fasta")
        )
        job.name = "rnammer_split"
        job.subname = "rnammer_split"
        jobs.append(job)
        
        #prefixes = ["rnammer_5S", "rnammer_8S", "rnammer_16S", "rnammer_18S", "rnammer_23S", "rnammer_28S"]
        #prefixes = ["rnammer_5S", "rnammer_16S", "rnammer_23S"]
        prefixes = ["rnammer_16S"]
            
        # get ncbi id with blastn against nt AND perform RDP classifier.
        # Note that blastn otu tables will be generated in a later step: taxonomy_annotation.
        for prefix in prefixes:
            # create bed files for abundance - we have to only get the rea
            # Just take taxonomy file which contains coord for each rRNA genes, then generate bed file.
            job = shotgun_metagenomics.rnammer_to_bed(
                os.path.join(self._root_dir, "gene_annotation", "rnammer", prefix + ".fasta"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", prefix + ".bed")
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
                    os.path.join(self._root_dir, "gene_annotation", "rnammer", prefix + ".bed"),
                    os.path.join(self._root_dir, "gene_annotation", "rnammer", "abundance", prefix + "_" + readset.name + ".cov"),
                    flag
                )
                job.name = "bedtoolsCov-contigs-rnammer_" + readset.sample.name + "_" + prefix
                job.subname = "bedtools"
                jobs.append(job)
                cov_list.append(os.path.join(self._root_dir, "gene_annotation", "rnammer", "abundance", prefix + "_" + readset.name + ".cov"))
        
            job = shotgun_metagenomics.merge_counts(
                cov_list,
                os.path.join(self._root_dir, "gene_annotation", "rnammer", prefix + "_merged_abundance.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", prefix + "_merged_abundance_RPKM.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", prefix + "_merged_abundance_cpm.tsv"),
                "genes"
            )
            job.name = "merge_gene_abundance_contigs_" + prefix
            job.subname = "merge_gene_abundance"
            jobs.append(job)
            
            # RDP classifier with rnammer results.
            job = microbial_ecology.rdp_wrapper(
                os.path.join(self._root_dir, "gene_annotation", "rnammer", prefix + ".fasta"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_rdp", prefix + "_rdp.tsv")
            )
            job.name = "classify_" + prefix
            job.subname = "RDP"
            jobs.append(job)
           
            # Convert RDP table to in-house taxonomy format.
            job = shotgun_metagenomics.rdp_to_taxonomy(
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_rdp", prefix + "_rdp.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_rdp", prefix + "_rdp_taxonomy.tsv")
            )
            job.name = "rdp_to_taxonomy_rnammer_" + prefix
            job.subname = "rdp_to_taxonomy"
            jobs.append(job)

            # Then with rdp output, generate otu table
            job = shotgun_metagenomics.generate_otu_table(
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_rdp", prefix + "_rdp_taxonomy.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", prefix + "_merged_abundance.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_rdp", prefix +  "_otu_table.tsv")
            )
            job.name = "generate_otu_table"
            job.subname = "generate_otu_table"
            jobs.append(job)
            
            # Then, once abundance is done, BLASTN on nt NCBI. OTU table generation will be done later in taxonomic_annotation step!
            job = shotgun_metagenomics.blastn(
                os.path.join(self._root_dir, "gene_annotation", "rnammer", prefix + ".fasta"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_blastn", "blastn_nt_" + prefix + ".tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_blastn"),
                "blastn"
            )
            job.name = "blastn_nt_" + prefix
            job.subname = "blastn"
            jobs.append(job)
            
            job = shotgun_metagenomics.extract_taxonomy(
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_blastn", "blastn_nt_" + prefix + ".tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_blastn", "blastn_nt_taxonomy_" + prefix + ".tsv")
            )
            job.name = "extract_taxonomy_rnammer_blastn_" + prefix
            job.subname = "ncbi_tax"
            jobs.append(job)
   
            job = shotgun_metagenomics.generate_otu_table(
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_blastn", "blastn_nt_taxonomy_" + prefix + ".tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", prefix + "_merged_abundance.tsv"),
                os.path.join(self._root_dir, "gene_annotation", "rnammer", "rnammer_blastn", prefix + "_otu_table.tsv")
            )
            job.name = "generate_otu_table_rnammer_blastn_" + prefix
            job.subname = "generate_otu_table"
            jobs.append(job)

        return jobs
