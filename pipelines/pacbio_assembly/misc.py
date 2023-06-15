
    def pacbio_tools_assembly_stats(self):

        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)
                preassembly_directory = os.path.join(coverage_directory, "preassembly", "data")

                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list'):
                    mer_size_text = "merSize" + mer_size
                    sample_cutoff_mer_size = "_".join([sample.name, cutoff_x, mer_size_text])
                    mer_size_directory = os.path.join(coverage_directory, mer_size_text)
                    blast_directory = os.path.join(mer_size_directory, "blast")
                    mummer_file_prefix = os.path.join(mer_size_directory, "mummer", sample.name + ".")
                    report_directory = os.path.join(mer_size_directory, "report")

                    polishing_rounds = config.param('DEFAULT', 'polishing_rounds', type='posint')
                    fasta_consensus = os.path.join(mer_size_directory, "polishing" + str(polishing_rounds), "data", "consensus.fasta")
                    # Number of unique run-smartcells per sample
                    smartcells = len(set([(readset.run, readset.smartcell) for readset in sample.readsets]))

                    jobs.append(concat_jobs([
                        Job(command="mkdir -p " + report_directory),
                        # Generate table(s) and figures
                        pacbio_tools.assembly_stats(
                            os.path.join(preassembly_directory, "filtered_shortreads.fa"),
                            os.path.join(preassembly_directory, "filtered_longreads.fa"),
                            os.path.join(preassembly_directory, "corrected.fasta"),
                            os.path.join(sample.name, "filtering", "data", "filtered_summary.csv"),
                            fasta_consensus,
                            sample.name,
                            cutoff_x + "_" + mer_size,
                            sample.readsets[0].estimated_genome_size,
                            smartcells,
                            report_directory
                        )
                    ], name="pacbio_tools_assembly_stats." + sample_cutoff_mer_size))

                    report_file = os.path.join(report_directory, "PacBioAssembly.pacbio_tools_assembly_stats.md")
                    jobs.append(
                        Job(
                            [
                                os.path.join(report_directory, "summaryTableReads.tsv"),
                                os.path.join(report_directory, "summaryTableReads2.tsv"),
                                os.path.join(report_directory, "summaryTableAssembly.tsv"),
                                os.path.join(report_directory, "pacBioGraph_readLengthScore.pdf"),
                                os.path.join(report_directory, "pacBioGraph_readLengthScore.jpeg"),
                                os.path.join(report_directory, "pacBioGraph_histoReadLength.pdf"),
                                os.path.join(report_directory, "pacBioGraph_histoReadLength.jpeg"),
                                fasta_consensus + ".gz"
                            ],
                            [report_file],
                            [['pacbio_tools_assembly_stats', 'module_pandoc']],
                            command="""\
cp {fasta_consensus}.gz {report_directory}/ && \\
total_subreads=`grep -P '^"Total subreads"\t"' {report_directory}/summaryTableReads.tsv | cut -f2 | sed 's/"//g'` && \\
average_subreads_length=`grep -P '^"Average subread length"\t"' {report_directory}/summaryTableReads.tsv | cut -f2 | sed 's/"//g'` && \\
summary_table_reads2=`LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, sprintf("%\\47d", $2), sprintf("%\\47d", $3), sprintf("%\\47d", $4), sprintf("%\\47d", $5), sprintf("%\\47d", $6), sprintf("%\\47d", $7), sprintf("%\\47d", $8), sprintf("%\\47d", $9), sprintf("%\\47d", $10), sprintf("%\\47d", $11)}}}}' {report_directory}/summaryTableReads2.tsv | sed 's/^#//'` && \\
summary_table_assembly=`awk -F"\t" '{{OFS="\t"; if (NR==1) {{print; gsub(/[^\t]/, "-")}} print}}' {report_directory}/summaryTableAssembly.tsv | sed 's/"//g' | sed 's/\t/|/g'` && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable total_subreads="$total_subreads" \\
  --variable average_subreads_length="$average_subreads_length" \\
  --variable summary_table_reads2="$summary_table_reads2" \\
  --variable summary_table_assembly="$summary_table_assembly" \\
  --variable smartcells="{smartcells}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                                fasta_consensus=fasta_consensus,
                                report_directory=report_directory,
                                smartcells=smartcells,
                                report_template_dir=self.report_template_dir,
                                basename_report_file=os.path.basename(report_file),
                                report_file=report_file
                            ),
                            report_files=[os.path.relpath(report_file, mer_size_directory)],
                            name="pacbio_tools_assembly_stats_report." + sample_cutoff_mer_size)
                    )

        return jobs
