
    def gene_abundance(self):
        jobs = []
        cov_list = []

        # if doing assembly by samples + big assembly.
        if(self.big_assembly_only == False):
	        for readset in self.readsets:
	            bam = os.path.join(self._root_dir, "gene_abundance", readset.sample.name, readset.name + ".bam")
	            cov = os.path.join(self._root_dir, "assembly", readset.sample.name, "Scaffolds.cov")
	            bed = os.path.join(self._root_dir, "assembly", readset.sample.name, "Scaffolds.bed")
	            reference = os.path.join(self._root_dir, "assembly", readset.sample.name, "Scaffolds.fasta")
	            #cov_list.append(cov)
	    
	            sys.stderr.write("OUTDIR IN GENEABUN: " + reference + "\n")
	
	            outdir = self._root_dir + "gene_abundance/" + readset.sample.name
	            if not os.path.exists(outdir):
	                os.makedirs(os.path.join(outdir))
	       
	            # Will make index for bwa. and also bed file for computing reads spanning later.
	            job = shotgun_metagenomics.make_index(
	                reference
	            )
	            job.name = "make_index"
	            job.subname = "make_index"
	            jobs.append(job)
	            
	            if( readset.run_type == "PAIRED_END"):
	                infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired.fastq.gz")
	                out1 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R1.fastq.gz")
	                out2 = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_paired_R2.fastq.gz")
	                
	                job = shotgun_metagenomics.split_pairs(
	                    infile,
	                    out1,
	                    out2
	                )
	                job.name = "split_pairs_" + readset.sample.name
	                job.subname = "split_pairs"
	                jobs.append(job)
	
	                job = shotgun_metagenomics.bwa_mem(
	                    reference,
	                    out1,
	                    out2,
	                    bam
	                )
	                job.name = "bwa_mem-" + readset.sample.name
	                job.subname = "bwa"
	                jobs.append(job)
	
	            elif(readset.run_type == "SINGLE_END"):
	                infile = os.path.join(self._root_dir, "qced_reads", readset.sample.name, readset.name + ".ncontam_R1.fastq.gz")
	                
	                job = shotgun_metagenomics.bwa_mem_se(
	                    reference,
	                    infile,
	                    bam
	                )
	                job.name = "bwa_mem-" + readset.sample.name
	                job.subname = "bwa"
	                jobs.append(job)
	
	            job = shotgun_metagenomics.coverage_bed(
	                bam,
	                #self._root_dir + "gene_clustering/Scaffolds.bed",
	                bed,
	                cov
	            )
	            job.name = "bedtoolsCov-" + readset.sample.name
	            job.subname = "bedtools"
	            jobs.append(job)
            

        # Once all coverage has been computed, merge all tables.
        #sys.stderr.write('gene abundance: ' + ','.join([str(x) for x in cov_list] ) + '\n')
        job = shotgun_metagenomics.merge_counts(
            cov_list,
            self._root_dir + "gene_abundance/merged_gene_abundance.tsv"
        )
        job.name = "merge_gene_abundance"
        job.subname = "merge_gene_abundance"
        jobs.append(job)
         
        return jobs 
