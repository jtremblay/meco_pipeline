
def clustering1(infile, barcodes, outdir):
        
    job = Job(
        [infile],
        ["outdir/obs_filtered.fasta", "outdir/obs_filtered.tsv"],
        [
            ['usearch', 'moduleVersion.usearch'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl']
        ]
    )
    job.command="""
clustering1.pl \\
--infile_fastq {infile} \\
--ref_db {ref_db} \\
--barcodes {barcodes} \\
--outdir {outdir} \\
--num_threads {num_threads}.""".format(
    infile = infile,
    barcodes = barcodes,
    ref_db =  config.param( 'DB', 'chimeras', 'filepath'),
    outdir = outdir,
    num_threads = config.param( 'clustering', 'num_threads', 1, 'int')
    )

    return job

def clustering2(infile, barcodes, outdir):
        
    job = Job(
        [infile],
        ["outdir/obs_filtered.fasta", "outdir/obs_filtered.tsv"],
        [
            ['usearch', 'moduleVersion.usearch'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl']
        ]
    )
    job.command="""
clustering2.pl \\
--infile_fastq {infile} \\
--ref_db {ref_db} \\
--barcodes {barcodes} \\
--outdir {outdir} \\
--num_threads {num_threads}""".format(
        infile = infile,
        ref_db = config.param( 'DB', 'chimeras', 1, 'path'),
        barcodes = barcodes,
        outdir = outdir,
        num_threads = config.param( 'clustering', 'num_threads', 'int')
    )
    return job


def kmernator(infiles, outfiles, outdir):
        
    #outfiles = []
    #for infile in infiles:
    #    basename =  os.path.splitext(os.path.basename(infile))[0]
    #    outfile = os.path.join(outdir, "FilteredData" + basename + ".fastq")
    #    outfiles.append(outfile)

    job = Job(
        infiles,
        outfiles,
        [
            ['memtime', 'moduleVersion.memtime'],
            ['meco_tools', 'moduleVersion.meco_tools'],
            ['kmernator', 'moduleVersion.kmernator'],
            ['gcc', 'moduleVersion.gcc'],
            ['openmpi', 'moduleVersion.openmpi']
        ]
    )

    job.command="""
memtime mpirun -n {num_threads} FilterReads-P \\
  --output-file {outdir}/FilteredData \\
  --max-kmer-output-depth {max_depth} \\
  --min-depth {min_depth} \\
  --skip-artifact-filter=1 \\
  {kmer_size} \\
  {infiles}""".format(
    num_threads = config.param("kmernator", "num_threads", 1, "posint"),
    infiles = " ".join(infiles),
    max_depth = config.param("kmernator", "max_depth", 1, "posint"),
    min_depth = config.param("kmernator", "min_depth", 1, "posint"),
    kmer_size = config.param("kmernator", "kmer_size", 1, "posint"),
    outdir = outdir
     )

    return job

def ray_big(infiles, outdir, type):
   
    infile_string = ""
    if(type == "pe"):
        for infile in infiles:
            infile_string = infile_string + " -i " + infile

    if(type == "se"):
        for infile in infiles:
            infile_string = infile_string + " -s " + infile;
    
    job = Job(
        infiles,
        [outdir + "/Contigs.fasta"],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['gcc', 'moduleVersion.gcc'],
            ['ray', 'moduleVersion.ray'],
            ['openmpi', 'moduleVersion.openmpi']
        ]
    )
    job.command="""
rm {outdir} -rf && \\
memtime mpiexec -n {num_threads} Ray -k {kmer} \\
  -minimum-contig-length {min_contig_length} \\
  {infiles} \\
  -o {outdir}""".format(
    min_contig_length = config.param("ray_big", "minContigLength", 1, "posint"),
    num_threads = config.param("ray_big", "num_threads", 1, "posint"),
    kmer = config.param("ray_big", "kmer", 1, "posint"),
    #routing_graph_degree = config.param("ray", "routingGraphDegree", 1, "posint"),
    #connection_type = config.param("ray", "connectionType", 1, "string"),
    infiles = infile_string,
    outdir = outdir
  )
  #-write-kmers \\
  #-route-messages \\
  #-connection-type {connection_type} \\
  #-routing-graph degree {routing_graph_degree} \\

    return job



#def merge_trimstats(infiles, outfile):
#        
#    job = Job(
#        infiles,
#        [outfile],
#        [
#            ['memtime', 'moduleVersion.memtime'],
#            ['meco_tools', 'moduleVersion.meco_tools']
#        ]
#    )
#    job.command="""
#memtime samtools flagstats {infile} > {outfile}""".format(
#    infile = " ".join(infiles),
#    outfile = outfile
#    )


#def bwa_mem(reference, infile1, infile2, outfile):
#    job = Job(
#        [infile1, infile2, reference],
#        [outfile],
#        [
#            ['memtime', 'moduleVersion.memtime'],
#            ['java', 'moduleVersion.java'],
#            ['picard', 'moduleVersion.picard'],
#            ['samtools', 'moduleVersion.samtools'],
#            ['bwa', 'moduleVersion.bwa']
#        ]
#    )
#    job.command="""
#memtime bwa mem -M \\
#  -t {num_threads} \\
#  {reference} \\
#  {infile1} \\
#  {infile2} \\
#  | java \\
#  -Djava.io.tmpdir=/localscratch/ \\
#  -XX:ParallelGCThreads=1 \\
#  -Dsamjdk.use_async_io=true \\
#  -Dsamjdk.buffer_size=4194304 \\
#  -Xmx15G -jar \${{PICARD_HOME}}/SortSam.jar \\
#  INPUT=/dev/stdin \\
#  VALIDATION_STRINGENCY=SILENT \\
#  OUTPUT={outfile} \\
#  SORT_ORDER=coordinate \\
#  CREATE_INDEX=true \\
#  CREATE_MD5_FILE=true \\
#  MAX_RECORDS_IN_RAM={max_records_in_ram} \\
#  && samtools view -F 0x100 -bh {outfile} > {outfile}.tmp && mv {outfile}.tmp {outfile} \\
#  && memtime samtools index {outfile}""".format(
#    num_threads = config.param('bwa', 'num_threads', type='int', required=True), 
#    max_records_in_ram = config.param('bwa', 'max_records_in_ram', type='int', required=True),
#    infile1 = infile1,
#    infile2 = infile2,
#    reference = reference,
#    outfile = outfile
#    )
#  #SORT_ORDER=coordinate \\
#
#    return job

#def bwa_mem(reference, bwt, infile1, infile2):
#    out_sam = None
#    job = Job(
#        [infile1, infile2, reference, bwt],
#        [out_sam],
#        [
#            ['memtime', 'moduleVersion.memtime'],
#            ['java', 'moduleVersion.java'],
#            ['picard', 'moduleVersion.picard'],
#            ['samtools', 'moduleVersion.samtools'],
#            ['bwa', 'moduleVersion.bwa']
#        ]
#    )
#    job.command="""
#memtime bwa mem -M \\
#  -t {num_threads} \\
#  {reference} \\
#  {infile1} \\
#  {infile2} \\
#  {out_sam}""".format(
#    num_threads = config.param('bwa', 'num_threads', type='int', required=True), 
#    mem_per_thread = config.param('samtools', 'mem_per_thread', required=True),
#    infile1 = infile1,
#    infile2 = infile2,
#    reference = reference,
#    out_sam=" \\  > " + out_sam if out_sam else ""
#    )
#
#    return job
#
#def samtools(input, output):
#    job = Job(
#        [input],
#        [output],
#        [
#            ['memtime', 'moduleVersion.memtime'],
#            ['samtools', 'moduleVersion.samtools']
#        ]
#    )
#    job.command="""
#  samtools view -Sbh -F 0x100 {input} > {output}.tmp && \\
#  memtime samtools sort -@ {num_threads} -m {mem_per_thread} {output}.tmp {output}.tmp.sorted && \\
#  memtime mv {output}.tmp.sorted.bam {output} && \\
#  memtime rm {output}.tmp && \\
#  memtime samtools index {output}""".format(
#    num_threads = config.param('bwa', 'num_threads', type='int', required=True), 
#    mem_per_thread = config.param('samtools', 'mem_per_thread', required=True),
#    input = input,
#    output = output
#    )
#    
#    return job

#def bwa_mem_canopies(reference, infile1, infile2, outfile):
#    job = Job(
#        [infile1, infile2],
#        [outfile],
#        [
#            ['memtime', 'moduleVersion.memtime'],
#            ['java', 'moduleVersion.java'],
#            ['picard', 'moduleVersion.picard'],
#            ['samtools', 'moduleVersion.samtools'],
#            ['bwa', 'moduleVersion.bwa']
#        ]
#    )
#    job.command="""
#memtime bwa mem -M \\
#  -t {num_threads} \\
#  {reference} \\
#  {infile1} \\
#  {infile2} \\
#  | java \\
#  -Djava.io.tmpdir=/localscratch/ \\
#  -XX:ParallelGCThreads=1 \\
#  -Dsamjdk.use_async_io=true \\
#  -Dsamjdk.buffer_size=4194304 \\
#  -Xmx15G -jar \${{PICARD_HOME}}/SortSam.jar \\
#  INPUT=/dev/stdin \\
#  VALIDATION_STRINGENCY=SILENT \\
#  OUTPUT={outfile} \\
#  SORT_ORDER=queryname \\
#  CREATE_INDEX=true \\
#  CREATE_MD5_FILE=true \\
#  MAX_RECORDS_IN_RAM={max_records_in_ram} \\
#  && samtools view -F 0x100 -bh {outfile} > {outfile}.tmp && mv {outfile}.tmp {outfile}""".format(
#    num_threads = config.param('bwa', 'num_threads', type='int', required=True), 
#    max_records_in_ram = config.param('bwa', 'max_records_in_ram', type='int', required=True),
#    infile1 = infile1,
#    infile2 = infile2,
#    reference = reference,
#    outfile = outfile
#    )
#  #SORT_ORDER=coordinate \\
#
#    return job


def canopy(infile, prefix, out_clusters, out_profiles):
    job = Job(
        [infile],
        [out_clusters, out_profiles],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.meco_tools'],
            ['canopy', 'moduleVersion.canopy']
        ]
    )
    job.command="""
memtime cc_x86.bin \\
  -i {infile} \\
  -c {out_profiles} \\
  -o {out_clusters} \\
  -p {prefix} \\
  -n {num_threads} \\
  --max_canopy_dist {max_canopy_dist} \\
  --max_merge_dist {max_merge_dist}""".format(
    infile = infile,
    out_clusters = out_clusters,
    out_profiles = out_profiles,
    prefix = prefix,
    max_canopy_dist = config.param('canopy', 'max_canopy_dist', required=True, type='float'),
    max_merge_dist = config.param('canopy', 'max_merge_dist', type='float', required=True),
    num_threads = config.param('canopy', 'num_threads', type='int', required=True) 
    )

    return job

def extract_canopy_contigs(contigs, outdir, canopies):
    #dummy_outfile = os.path.join(outdir, "extract_canopy_contigs.done")
    outfile1 = os.path.join(outdir, "canopy_contigs.fasta")
    outfile2 = os.path.join(outdir, "canopy_contigs.fasta.bwt")

    job = Job(
        [canopies],
        [outfile1, outfile2],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.meco_tools'],
            ['bwa', 'moduleVersion.bwa']
        ]
    )
    job.command="""
memtime contigsToCanopy.pl \\
  --infile_canopy {canopies} \\
  --contigs {contigs} \\
  --outdir {outdir} \\
  --n {n}""".format(
    canopies = canopies,
    contigs = ",".join(contigs),
    outdir = outdir,
    n = config.param('extract_canopy_contigs', 'n', 'int')
    )

    return job

def parse_canopies(infile, outfile):
    job = Job(
        [infile],
        [outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['meco_tools', 'moduleVersion.meco_tools']
        ]
    )

    job.command="""
memtime parseCanopies.pl \\
  --infile {infile} \\
  --n {n} \\
  > {outfile}""".format(
    infile = infile,
    outfile = outfile,
    n = config.param('extract_canopy_contigs', 'n', 1, 'int')
    )

    return job

#def filter_canopy(in_clusters, in_profiles, out_clusters, out_profiles, out_abundance, stats):
#    job = Job(
#        [in_clusters, in_profiles],
#        [out_clusters, out_profiles],
#        [
#            ['memtime', 'moduleVersion.memtime'],
#            ['tools', 'moduleVersion.meco_tools']
#        ]
#    )
#    job.command="""
#memtime filterCanopy.pl \\
#  --infile_clusters {in_clusters} \\
#  --infile_profiles {in_profiles} \\
#  --outfile_clusters {out_clusters} \\
#  --outfile_profiles {out_profiles} \\
#  --outfile_canopy_abundance {out_abundance} \\
#  --min_number_of_genes {min_number_of_genes} \\
#  --high_abundance_perc {high_abundance_perc} \\
#  --high_abundance_freq {high_abundance_freq} \\
#   > {stats}""".format(
#    in_clusters = in_clusters,
#    in_profiles = in_profiles,
#    out_clusters = out_clusters,
#    out_profiles = out_profiles,
#    out_abundance = out_abundance,
#    stats = stats,
#    min_number_of_genes =  config.param('filter_canopy', 'min_number_of_genes', 'int'),
#    high_abundance_perc =  config.param('filter_canopy', 'high_abundance_perc', 'int'),
#    high_abundance_freq =  config.param('filter_canopy', 'high_abundance_freq', 'int')
#    )
#
#    return job


def extract_reads1(bams, outdir, clusters):
    dummy_outfile = outdir + "/extract_reads1.done"
    #sys.stderr.write(",".join(bams) +'\n')

    infiles = list(bams)
    infiles.append(clusters)
        
    job = Job(
        infiles, 
        [dummy_outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['samtools', 'moduleVersion.samtools'],
            ['meco_tools', 'moduleVersion.meco_tools']
        ]
    )
    
    job.command="""
memtime extractReads.pl \\
  --infiles {bams} \\
  --outdir {outdir} \\
  --clusters {clusters} \\
  --flag 0x040 \\
  && touch {dummy_outfile}""".format(
    bams = ",".join(bams),
    outdir = outdir,
    clusters = clusters,
    dummy_outfile = dummy_outfile
    #num_threads = config.param('extract_reads', 'num_threads', 'int')
    )

    return job

def extract_reads2(bams, outdir, clusters):
    dummy_outfile = outdir + "/extract_reads2.done"

    infiles = list(bams)
    infiles.append(clusters)
        
    job = Job(
        infiles, 
        [dummy_outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['samtools', 'moduleVersion.samtools'],
            ['meco_tools', 'moduleVersion.meco_tools']
        ]
    )
    
    job.command="""
memtime extractReads.pl \\
  --infiles {bam_infiles} \\
  --outdir {outdir} \\
  --clusters {clusters} \\
  --flag 0x080 \\
  && touch {dummy_outfile}""".format(
    bam_infiles = ",".join(bams),
    outdir = outdir,
    clusters = clusters,
    dummy_outfile = dummy_outfile
    #num_threads = config.param('extract_reads', 'num_threads', 'int')
    )

    return job


def kegg_overrep_DDA(indir, blastp_table, type, DDA_done, prefix):
   

    script = ""
    dummy_outfile = ""
    if type == "pathways":
        script = "getKeggPathwaysLoop.py"
        dummy_outfile = os.path.join(indir, "kegg_pathways_overrep.done")
    if type == "modules":
        script = "getKeggModulesLoop.py"
        dummy_outfile = os.path.join(indir, "kegg_modules_overrep.done")
    if type == "K":
        script = "getKeggKLoop.py"
        dummy_outfile = os.path.join(indir, "kegg_K_overrep.done")

    job = Job(
        [blastp_table, DDA_done],
        [dummy_outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['lapack', 'moduleVersion.lapack'],
            ['qiime-dependencies', 'moduleVersion.qiime-dependencies']
        ]
    )
  #--infile-gene-abundance {gene_abundance} \\
    job.command="""
memtime {script} \\
  --infile-blastp {blastp_table} \\
  --indir {indir} \\
  --num-threads {num_threads} \\
  --prefix {prefix} \\
  && touch {dummy_outfile}""".format(
  blastp_table = blastp_table,
  indir = indir,
  script = script,
  dummy_outfile = dummy_outfile,
  num_threads = config.param('getKegg', 'num_threads', 1, 'int'),
  prefix = prefix
  )

    return job



def cog_overrep(indir, rpsblast_table, DDA_done, prefix):
   
    dummy_outfile = os.path.join(indir, "cog_overrep.done")

    job = Job(
        [rpsblast_table, DDA_done],
        [dummy_outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['lapack', 'moduleVersion.lapack'],
            ['qiime-dependencies', 'moduleVersion.qiime-dependencies']
        ]
    )
  #--infile-gene-abundance {gene_abundance} \\
    job.command="""
memtime getCOGLoop.py \\
  --infile-cog {rpsblast_table} \\
  --indir {indir} \\
  --num-threads {num_threads} \\
  --prefix {prefix} \\
  && touch {dummy_outfile}""".format(
  rpsblast_table = rpsblast_table,
  indir = indir,
  dummy_outfile = dummy_outfile,
  num_threads = config.param('getCOG', 'num_threads', 1, 'int'),
  prefix = prefix
  )
    return job

def kog_overrep(indir, rpsblast_table, DDA_done, prefix):
   
    dummy_outfile = os.path.join(indir, "kog_overrep.done")

    job = Job(
        [rpsblast_table, DDA_done],
        [dummy_outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['lapack', 'moduleVersion.lapack'],
            ['qiime-dependencies', 'moduleVersion.qiime-dependencies']
        ]
    )
  #--infile-gene-abundance {gene_abundance} \\
    job.command="""
memtime getKOGLoop.py \\
  --infile-kog {rpsblast_table} \\
  --indir {indir} \\
  --num-threads {num_threads} \\
  --prefix {prefix} \\
  && touch {dummy_outfile}""".format(
  rpsblast_table = rpsblast_table,
  indir = indir,
  dummy_outfile = dummy_outfile,
  num_threads = config.param('getKOG', 'num_threads', 1, 'int'),
  prefix = prefix
  )
    return job

def pfam_overrep(indir, hmmscan_output, DDA_done, prefix):
   
    dummy_outfile = os.path.join(indir, "pfam_overrep.done")

    job = Job(
        [hmmscan_output, DDA_done],
        [dummy_outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['lapack', 'moduleVersion.lapack'],
            ['qiime-dependencies', 'moduleVersion.qiime-dependencies']
        ]
    )
  #--infile-gene-abundance {gene_abundance} \\
    job.command="""
memtime getPfamLoop.py \\
  --infile-pfam {hmmscan_output} \\
  --indir {indir} \\
  --num-threads {num_threads} \\
  --prefix {prefix} \\
  && touch {dummy_outfile}""".format(
  hmmscan_output = hmmscan_output,
  indir = indir,
  dummy_outfile = dummy_outfile,
  num_threads = config.param('getKOG', 'num_threads', 1, 'int'),
  prefix = prefix
  )
    return job

def tigrfam_overrep(indir, hmmscan_output, DDA_done, prefix):
   
    dummy_outfile = os.path.join(indir, "tigrfam_overrep.done")

    job = Job(
        [hmmscan_output, DDA_done],
        [dummy_outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['lapack', 'moduleVersion.lapack'],
            ['qiime-dependencies', 'moduleVersion.qiime-dependencies']
        ]
    )
  #--infile-gene-abundance {gene_abundance} \\
    job.command="""
memtime getTigrfamLoop.py \\
  --infile-pfam {hmmscan_output} \\
  --indir {indir} \\
  --num-threads {num_threads} \\
  --prefix {prefix} \\
  && touch {dummy_outfile}""".format(
  hmmscan_output = hmmscan_output,
  indir = indir,
  dummy_outfile = dummy_outfile,
  num_threads = config.param('getKOG', 'num_threads', 1, 'int'),
  prefix = prefix
  )
    return job



def upgma_clustering(infile, outfile):
    #outdir = os.path.abspath(infile)
    outdir = os.path.dirname(os.path.relpath(infile))
    outdir = os.path.join(outdir, "upgma")
    outdir = "./" + outdir + "/"
    #prefix = os.path.basename(infile)
    
    base = os.path.basename(outfile)
    prefix = os.path.splitext(base)[0]
    
    job = Job(
        [infile],
        [outfile],
        [
            ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
            ['qiime', 'moduleVersion.qiime'],
            ['python', 'moduleVersion.python'],
            ['R', 'moduleVersion.R'],
            ['gcc', 'moduleVersion.gcc'],
            ['tools', 'moduleVersion.tools'],
            ['lapack', 'moduleVersion.lapack']
        ]
    )

    job.command="""
export OMPI_MCA_mtl=^psm && \\
upgma_cluster.py \\
  -i {infile} \\
  -o {outfile} \\
  && plotTree.R \\
  -t {outfile} \\
  -o {outdir} \\
  -p {prefix} \\
  && plotTreeWithMappingFile.R \\
  -t {outfile} \\
  -o {outdir} \\
  -p {prefix} \\
  -m {mapping_file} \\
  && unset OMPI_MCA_mtl""".format(
    infile = infile,
    outfile = outfile,
    prefix = prefix,
    outdir = outdir,
    mapping_file = config.param('default', 'mappingFile', 1, 'filepath')
    )

    return job
