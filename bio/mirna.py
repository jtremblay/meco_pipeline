#!/usr/bin/env python

#LICENSE AND COPYRIGHT

#Copyright (C) 2022 INRS - Centre Armand-Frappier

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

# CAF Modules
from core.config import *
from core.job import *

def trimmomatic_se(input1, output1, quality_offset, trim_log, trim_stats):

    job = Job(
        [input1], 
        [output1, trim_log, trim_stats],
        [
            ['java', 'module_java'], 
            ['trimmomatic', 'module_trimmomatic']
        ]
    )
    threads = config.param('trim', 'threads', type='posint')
    adapter_file = config.param('trim', 'adapter_fasta', required=False, type='filepath')
    illumina_clip_settings = config.param('trim', 'illumina_clip_settings')
    trailing_min_quality = config.param('trim', 'trailing_min_quality', type='int')
    min_length = config.param('trim', 'min_length', type='posint')
    headcrop = config.param('trim', 'headcrop', required=False, type='int')
    average_quality = config.param('trim', 'average_quality', required=False, type='posint')

    job.command = """
java -XX:ParallelGCThreads={threads} -Xmx2G -jar \$TRIMMOMATIC_JAR {mode} \\
  -threads {threads} \\
  -phred{quality_offset} \\
  {input1} {output1} \\
  ILLUMINACLIP:{adapter_file}{illumina_clip_settings} \\
  TRAILING:{trailing_min_quality} \\
  MINLEN:{min_length} \\
  HEADCROP:{headcrop} \\
  LEADING:{leading} \\
  AVGQUAL:{average_quality}""".format(
        mode = "SE",
        threads = threads,
        quality_offset = quality_offset,
        input1 = input1,
        output1 = output1,
        illumina_clip_settings=illumina_clip_settings,
        trailing_min_quality=trailing_min_quality,
        adapter_file = config.param('trim', 'adapter_fasta', type='filepath'),
        min_length = min_length,
        headcrop = str(headcrop),
        average_quality = average_quality,
        leading = config.param('trim', 'leading_min_quality', 1, type='posint')
    )

    job.command += " \\\n  2> " + trim_log

    # Compute statistics
    job.command += " && \\\ngrep ^Input " + trim_log + " | perl -pe 's/^Input Read Pairs: (\\d+).*Both Surviving: (\\d+).*Forward Only Surviving: (\\d+).*$/Raw Fragments,\\1\\nFragment Surviving,\\2\\nSingle Surviving,\\3/' > " + trim_stats

    return job

def size_select(infile_fastq, outfile_passed, outfile_failed):

    basename =  os.path.splitext(os.path.basename(outfile_failed))[0]
    dirname = os.path.dirname(outfile_failed)
   

    job = Job(
        [infile_fastq],
        [outfile_passed, outfile_failed], 
        [
            ['caf_tools', 'module_tools']
        ]
    )
    
  #--gi_to_tax {gi_to_tax} \\
    job.command="""
filterFastqByLengthGz.pl \\
  --infile {infile_fastq} \\
  --min_length {min_length} \\
  --max_length {max_length} \\
  --remove_Ns \\
  --outfile_gtcutoff {dirname}/{basename} | gzip -f > {outfile_passed} && \\
  gzip -f {dirname}/{basename}""".format(
        infile_fastq = infile_fastq,
        basename = basename,
        dirname = dirname,
        outfile_passed = outfile_passed,
        min_length = config.param('size_select', 'min', required=True, type='posint'),
        max_length = config.param('size_select', 'max', required=True, type='posint')
    )

    return job


def bwa_aln_samtools(reference, bwt, infile, outfile):
    
    basename =  os.path.splitext(os.path.basename(outfile))[0]
    dirname = os.path.dirname(outfile)
    
    job = Job(
        [infile, reference, bwt],
        [outfile],
        [
            ['samtools', 'module_samtools'],
            ['bwa', 'module_bwa']
        ]
    )
    job.command="""
bwa aln -n {mismatch} -l {seed} \\
  -t {num_threads} \\
  {reference} \\
  {infile} \\
  > {dirname}/{basename}.sai && \\
bwa samse {reference} {dirname}/{basename}.sai {infile} \\
  | samtools view -Sbh -F 0x100 - > {outfile}.tmp && \\
  samtools sort -@ {num_threads} -m {mem_per_thread} {outfile}.tmp -o {outfile}.tmp.sorted.bam && \\
  mv {outfile}.tmp.sorted.bam {outfile} && \\
  rm {outfile}.tmp && \\
  samtools index {outfile}""".format(
    num_threads = config.param('bwa', 'num_threads', type='int', required=True), 
    mismatch = config.param('bwa', 'mismatch', type='posint', required=True), 
    seed = config.param('bwa', 'seed', type='posint', required=True), 
    mem_per_thread = config.param('samtools', 'mem_per_thread', required=True),
    infile = infile,
    reference = reference,
    outfile = outfile,
    basename = basename,
    dirname = dirname
    )

    return job

def make_bowtie2_index(fasta, bwt):
    #basename =  os.path.splitext(os.path.basename(fasta))[0]
    #basename =  os.path.splitext(fasta)[0]

    job = Job(
        [fasta],
        #[basename + ".fasta.bwt", basename + ".bed"],
        #[basename + ".bed", bwt],
        [bwt],
        [
            
            ['samtools', 'module_samtools'],
            ['caf_tools', 'module_tools'],
            ['bwa', 'module_bwa']
        ]
    )
    job.command="""
bowtie2-build {fasta}""".format(
    fasta = fasta,
    bwt = bwt
    )

    return job

def bowtie2_samtools(bti, infile, outfile_bam, unaligned_fastq):
    bti_file = os.path.join(bti + ".1.bt2")

    job = Job(
        [infile, bti_file],
        [outfile_bam, unaligned_fastq],
        [
            ['samtools', 'module_samtools'],
            ['bowtie2', 'module_bowtie2']
        ]
    )
    job.command="""
bowtie2 -p {num_threads} --local -D 20 -R 3 -N 1 -L 5 -p 8 --gbar 1 -mp 3 --un-gz {unaligned_fastq} \\
  -x {bti} \\
  -U {infile} \\
  | samtools view -Sbh -F 0x100 - > {outfile_bam}.tmp && \\
  samtools sort -@ {num_threads} -m {mem_per_thread} {outfile_bam}.tmp -o {outfile_bam}.tmp.sorted.bam && \\
  mv {outfile_bam}.tmp.sorted.bam {outfile_bam} && \\
  rm {outfile_bam}.tmp && \\
  samtools index {outfile_bam}""".format(
    num_threads = config.param('bowtie2', 'num_threads', type='int', required=True), 
    mem_per_thread = config.param('samtools', 'mem_per_thread', required=True),
    infile = infile,
    outfile_bam = outfile_bam,
    unaligned_fastq = unaligned_fastq,
    bti = bti
    )

    return job


def extract_unmapped(infile, outfile_bam, outfile_fastq):
    
    basename =  os.path.splitext(os.path.basename(outfile_fastq))[0]
    dirname = os.path.dirname(outfile_fastq)
    
    job = Job(
        [infile],
        [outfile_bam, os.path.join(outfile_fastq + ".gz")],
        [
            ['samtools', 'module_samtools'],
            ['bedtools', 'module_bedtools']
        ]
    )
    job.command="""
samtools view -Sbh -f 4 {infile} > {outfile_bam} && \\
  bedtools bamtofastq -i {outfile_bam} -fq {dirname}/{basename}.fastq && \\
  gzip -f {dirname}/{basename}.fastq""".format(
    infile = infile,
    outfile_bam = outfile_bam,
    basename = basename,
    dirname = dirname
    #outfile_fastq = outfile_fastq
    )

    return job

def extract_mapped(infile, outfile_bam, outfile_fastq):
    
    basename =  os.path.splitext(os.path.basename(outfile_fastq))[0]
    dirname = os.path.dirname(outfile_fastq)
    
    job = Job(
        [infile],
        [outfile_bam, os.path.join(outfile_fastq + ".gz")],
        [
            ['samtools', 'module_samtools'],
            ['bedtools', 'module_bedtools']
        ]
    )
    job.command="""
samtools view -Sbh -F 4 {infile} > {outfile_bam} && \\
  bedtools bamtofastq -i {outfile_bam} -fq {dirname}/{basename}.fastq && \\
  gzip -f {dirname}/{basename}.fastq""".format(
    infile = infile,
    outfile_bam = outfile_bam,
    basename = basename,
    dirname = dirname
    #outfile_fastq = outfile_fastq
    )

    return job

def dereplicate(infiles, outfile):
    job = Job(
        infiles,
        [outfile],
        [
            ['caf_tools', 'module_tools']
        ]
    )
    job.command="""
dereplicateGz.pl --fastq {infiles} \\
        > {outfile}""".format(
    infiles = ",".join(infiles),
    outfile = outfile
    )

    return job

def generate_matrix(infile, outdir):
    job = Job(
        [infile],
        [os.path.join(outdir, "obs.tsv"), os.path.join(outdir, "obs.fasta")],
        [
            ['caf_tools', 'module_tools']
        ]
    )
    job.command="""
dereplicateWriteTable.pl \\
    --fasta {infile} \\
    --outdir {outdir}""".format(
    infile = infile,
    outdir = outdir
    )

    return job


def filter_by_abundance(infile_tsv, infile_fasta, outfile_tsv, outfile_fasta, outfile_fasta_rc):
    job = Job(
        [infile_tsv, infile_fasta],
        [outfile_tsv, outfile_fasta, outfile_fasta_rc],
        [
            ['caf_tools', 'module_tools']
        ]
    )
    job.command="""
filterMiRNAByAbundance.pl \\
 --infile_tsv {infile_tsv} \\
 --infile_fasta {infile_fasta} \\
 --outfile_tsv {outfile_tsv} \\
 --outfile_fasta {outfile_fasta} \\
 --min_reads {min_reads} && \\
revComp.pl --infile {outfile_fasta} --rna > {outfile_fasta_rc}""".format(
    infile_tsv = infile_tsv,
    infile_fasta = infile_fasta,
    outfile_fasta_rc = outfile_fasta_rc,
    outfile_tsv = outfile_tsv,
    outfile_fasta = outfile_fasta,
    min_reads = config.param("filter", "min_reads", 1, "posint")
    )

    return job

def normalize(infile, fasta, outfile_rpkm):
    
    dirname = os.path.dirname(fasta)
    
    job = Job(
        [infile, fasta],
        [outfile_rpkm, os.path.join(dirname, "obs_length.tsv")],
        [
            ['caf_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    job.command="""
cat {fasta} | awk '\$0 ~ \\">\\" {{if (NR > 1) {{print c;}} c=0;printf substr(\$0,2,100) \\"\\\\t\\"; }} \$0 !~ \\">\\" {{c+=length(\$0);}} END {{ print c; }}' > {dirname}/obs_length.tsv && \\
normalizeMiRNA.R \\
 -i {infile} \\
 -l {dirname}/obs_length.tsv \\
 -o {outfile_rpkm}""".format(
    infile = infile,
    fasta = fasta,
    outfile_rpkm = outfile_rpkm,
    dirname = dirname
    )

    return job

def normalize_mirna_table(infile, fasta, outfile_rpkm):
    
    dirname = os.path.dirname(fasta)
    
    job = Job(
        [infile, fasta],
        [outfile_rpkm, os.path.join(dirname, "obs_length.tsv")],
        [
            ['caf_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    job.command="""
cat {fasta} | awk '\$0 ~ \\">\\" {{if (NR > 1) {{print c;}} c=0;printf substr(\$0,2,100) \\"\\\\t\\"; }} \$0 !~ \\">\\" {{c+=length(\$0);}} END {{ print c; }}' > {dirname}/obs_length.tsv && \\
normalizeMiRNA2.R \\
 -i {infile} \\
 -l {dirname}/obs_length.tsv \\
 -o {outfile_rpkm}""".format(
    infile = infile,
    fasta = fasta,
    outfile_rpkm = outfile_rpkm,
    dirname = dirname
    )

    return job

# Because array job : touch done instead of relying on the standard done mechanism.
def blastn_array_job(indir, inprefix, outdir, outprefix, prefix, infiles, outfiles, dones, db):

    infile = os.path.join(indir, inprefix + "\$SLURM_ARRAY_TASK_ID2")
    outfile = os.path.join(outdir, outprefix + "\$SLURM_ARRAY_TASK_ID2.tsv")
    done = os.path.join(outdir, outprefix + "\$SLURM_ARRAY_TASK_ID2.done")

    all_outfiles = []
    for i in range(0, len(dones)):
        all_outfiles.append(dones[i])
    for i in range(0, len(outfiles)):
        all_outfiles.append(outfiles[i])
    

    job = Job(
        infiles,
        all_outfiles,
        [
            ['blast', 'module_blast'],
            ['caf_tools', 'module_tools']
        ]
    )

    job.command="""
SLURM_ARRAY_TASK_ID2=\$(printf "%07d" \$SLURM_ARRAY_TASK_ID)
mkdir -p {blast_dir} && \\
if [[ -s {infile} ]] ; then
    blastn \\
     -db {db} \\
     -query {infile} \\
     -out {outfile} \\
     -word_size {word_size} \\
     -outfmt \\"{outfmt}\\" \\
     -max_target_seqs 10 \\
     -num_threads {num_threads} \\
     -evalue {evalue} && \\
     touch {done} 
else
    touch {outfile}
fi ;""".format(
    infile = infile,
    outfile = outfile,
    #db = config.param(prefix, 'db', 1, 'string'),
    db = db,
    infiles = infiles,
    outfiles = outfiles,
    word_size = str(config.param(prefix, 'word_size', 1, 'int')),
    outfmt = config.param(prefix, 'outfmt', required=True),
    num_threads = config.param(prefix, 'num_threads', 1, 'int'),
    #other_params = config.param(prefix, 'other_params', 1, 'string'),
    blast_dir = outdir,
    evalue = config.param(prefix, "evalue", 1, "string"),
    SLURM_ARRAY_TASK_ID = "SLURM_ARRAY_TASK_ID",
    done = done
    #SLURM_ARRAY_TASK_ID2 = "SLURM_ARRAY_TASK_ID2"
    )

    return job

def blastn(infile, outfile, db):

    job = Job(
        [infile],
        [outfile],
        [
            ['blast', 'module_blast'],
            ['caf_tools', 'module_tools']
        ]
    )

    job.command="""blastn \\
     -db {db} \\
     -query {infile} \\
     -out {outfile} \\
     -word_size {word_size} \\
     -outfmt \\"{outfmt}\\" \\
     -max_target_seqs 5 \\
     -num_threads {num_threads} \\
     -evalue {evalue}""".format(
    infile = infile,
    outfile = outfile,
    #db = config.param("blastn", 'db', 1, 'string'),
    db = db,
    word_size = str(config.param("blastn", 'word_size', 1, 'int')),
    outfmt = config.param("blastn", 'outfmt', required=True),
    num_threads = config.param("blastn", 'num_threads', 1, 'int'),
    #other_params = config.param("blastn", 'other_params', 1, 'string'),
    evalue = config.param("blastn", "evalue", 1, "string")
    )

    return job

def filter_blast_hits(infile, outfile_filtered):

    job = Job(
        [infile],
        [outfile_filtered],
        [
            ['caf_tools', 'module_tools']
        ]
    )

    job.command="""
keepBlastBestHitsMiRNATargets.pl \\
    --infile {infile} \\
    --e {evalue} \\
    --length {length} \\
    --max_mismatch {mismatch} > {outfile_filtered}""".format(
    infile = infile,
    outfile_filtered = outfile_filtered,
    #outfile_filtered_best = outfile_filtered_best,
    evalue = config.param("blastn", "evalue", 1, "string"),
    length = config.param("blastn", "length", 1, "string"),
    mismatch = config.param("blastn", "mismatch", 1, "string"),
    perc = config.param("blastn", "perc", 1, "string"),
    n = config.param("blastn", "n", 1, "string")
    )

    return job

def filter_blast_hits_for_targets_finding(infile, outfile_filtered, outfile_filtered_best):

    job = Job(
        [infile],
        [outfile_filtered, outfile_filtered_best],
        [
            ['caf_tools', 'module_tools']
        ]
    )

    job.command="""
keepBlastBestHitsMiRNA.pl \\
    --infile {infile} \\
    --e {evalue} \\
    --length {length} \\
    --mismatch {mismatch} \\
    --perc {perc} \\
    > {outfile_filtered} && \\
keepNBestBlastHit.pl \\
    --infile {outfile_filtered} \\
    --n {n} \\
    > {outfile_filtered_best}""".format(
    infile = infile,
    outfile_filtered = outfile_filtered,
    outfile_filtered_best = outfile_filtered_best,
    evalue = config.param("blastn_targets", "evalue", 1, "string"),
    length = config.param("blastn_targets", "length", 1, "string"),
    mismatch = config.param("blastn_targets", "mismatch", 1, "string"),
    perc = config.param("blastn_targets", "perc", 1, "string"),
    n = config.param("blastn_targets", "n", 1, "string")
    )

    return job

def merge_tables(infile_abundance, infile_mature, outfile):

    job = Job(
        [infile_abundance, infile_mature],
        [outfile],
        [
            ['caf_tools', 'module_tools']
        ]
    )

    job.command="""
mergeMiRNATables.R \\
    -a {infile_abundance} \\
    -m {infile_mature} \\
    -o {outfile}""".format(
    infile_abundance = infile_abundance,
    #infile_hairpin = infile_hairpin,
    infile_mature = infile_mature,
    outfile = outfile
    )

    return job

def coverage_bed(infile, gff, outfile):
    job = Job(
        [infile, gff],
        [outfile],
        [
            ['bedtools', 'module_bedtools']
        ]
    )
    
    if isinstance( config.param('bedtools', 'other_args', 1, 'string'), str ):
        if config.param('bedtools', 'other_args', 1, 'string') != "":
            job.command="""
bedtools coverage {other_args} \\""".format(
                 other_args = config.param('bedtools', 'other_args', 1, 'string')
                )
        else:
            job.command="""
bedtools coverage \\"""

    job.command+="""
  -b {infile} \\
  -a {gff} \\
  -counts \\
  > {outfile}""".format(
    infile = infile,
    gff = gff,
    outfile = outfile
    )
    
    return job

def merge_counts(infiles, outfile_raw, outfile_cpm, outfile_rpkm):

    job = Job(
        infiles,
        [outfile_raw, outfile_cpm, outfile_rpkm],
        [
            ['tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
  
  #-n {outfile_normalized} \\
    job.command="""
mergeGeneAbundanceMiRNA.R \\
  -i {infile} \\
  -o {outfile_raw} \\
  -c {outfile_cpm} \\
  -r {outfile_rpkm}""".format(
        infile = ",".join(infiles),
        outfile_raw = outfile_raw,
        outfile_cpm = outfile_cpm,
        outfile_rpkm = outfile_rpkm
    )

    return job

def generate_fasta_file(infile, outfile, outfile_revcomp, ref):

    job = Job(
        [infile, ref],
        [outfile, outfile_revcomp],
        [
            ['tools', 'module_tools']
        ]
    )
  
  #-n {outfile_normalized} \\
    job.command="""
cat {infile} | awk '{{print \\$1}}' > {infile}.txt && \\
getFastaSeqsByListOfGenes.pl \\
  --infile {ref} \\
  --gene_list {infile}.txt \\
  > {outfile} && \\
revComp.pl --infile {outfile} --rna > {outfile_revcomp}""".format(
        infile = infile,
        outfile = outfile,
        outfile_revcomp = outfile_revcomp,
        ref = ref
    )

    return job

def ssearch(infile_fasta, infile_reference_dir, outfile):

    job = Job(
        [infile_fasta],
        [outfile],
        [
            ['caf_tools', 'module_tools'],
            ['fasta36', 'module_fasta36']
        ]
    )
#ssearch36 -f 2 -g 0.5 -m 8CC -U -E 1000 -D -T 8 -b'=5' -d'=5'
#-m {m} \\
    job.command="""
echo '' > {outfile} && \\
for file in {infile_reference_dir}/*.fna ; do
   ssearch36 \\
    -f {f} \\
    -g {g} \\
    -E {E} \\
    -T {t} \\
    -b {b} \\
    {flags} \\
    {infile_fasta} \\$file >> {outfile}
done""".format(
    f = config.param("ssearch", "gap_open", 1, "int"),
    g = config.param("ssearch", "gap_extension", 1, "float"),
    #m = config.param("ssearch", "out_format", 1, "string"),
    E = config.param("ssearch", "expectation", 1, "posint"),
    t = config.param("ssearch", "num_threads", 1, "posint"),
    b = config.param("ssearch", "high_scores_reported", 1, "posint"),
    #d = config.param("ssearch", "number_of_alignments_reported", 1, "posint"),
    flags = config.param("ssearch", "other_flags", 1, "string"),
    infile_fasta = infile_fasta,
    infile_reference_dir = infile_reference_dir,
    outfile = outfile
    )

    return job

def merge_chunks(indir, outfile, number_of_chunks, prefix):

    infiles = []
    for i in range(number_of_chunks):
        infiles.append(os.path.join(indir, prefix + "{:07d}.tsv".format(i)))

    job = Job(
        infiles, 
        [outfile],
        [
            ['caf_tools', 'module_tools']
        ]
    )
    
    job.command="""
cat {infiles} \\
  > {outfile}""".format(
    outfile = outfile,
    infiles = " ".join(infiles)
    )

    return job

def merge_chunks_generic(infiles, outfile):

    job = Job(
        infiles, 
        [outfile],
        [
            ['caf_tools', 'module_tools']
        ]
    )
    
    job.command="""
cat {infiles} \\
  > {outfile}""".format(
    outfile = outfile,
    infiles = " ".join(infiles)
    )

    return job

def rtk(infile, outdir, header_first_element, remove_last_col=False):
        
    job = Job(
        [infile],
        [os.path.join(outdir, "rtk.done")],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools'],
            ['rtk', 'module_rtk']
        ]
    )
    
    job.command="""
mkdir -p {outdir} && \\""".format(outdir = outdir)
    
#sed -i '1 s/^\\\\t/{header_first_element}\\\\t/' {infile} && \\
    job.command+="""
generateDepthPointsForRarefaction.pl \\"""

    if remove_last_col == True:
        job.command+="""
  --remove_last_col true \\"""
    else:
        job.command+="""\
  --remove_last_col false \\"""

    job.command+="""
  --infile {infile} \\
  --number_of_points {number_of_points} \\
  > {outdir}/depth_list.txt && \\
sed '1 s/^\\\S\\\+\\\\t/\\\\t/' {infile} > {outdir}/tmp.tsv && \\
rtk {mode} \\
  -i {outdir}/tmp.tsv \\
  -o {outdir}/ \\
  -r {perm} \\
  -d \\`cat {outdir}/depth_list.txt\\` && \\
rm {outdir}/tmp.tsv && \\
touch {outdir}/rtk.done""".format(
        infile = infile,
        outdir = outdir,
        header_first_element = header_first_element,
        number_of_points = config.param('alphadiv', 'number_of_points', 1, 'posint'),
        perm = config.param('alphadiv', 'perm', 1, 'posint'),
        mode = config.param('alphadiv', 'mode', 1, 'string')
    )

    return job

def parse_mirna_targets(infile, outfile_raw, outfile_parsed, strand):

    job = Job(
        [infile],
        [outfile_raw, outfile_parsed],
        [
            ['tools', 'module_tools']
        ]
    )
  
  #-n {outfile_normalized} \\
    if strand == "fwd_strand":
        job.command="""
parseSsearch.pl \\
  --infile {infile} \\
  > {outfile_raw} && \\
parseMiRNATargets.pl \\
  --infile {outfile_raw} \\
  > {outfile_parsed}""".format(
            infile = infile,
            outfile_raw = outfile_raw,
            outfile_parsed = outfile_parsed
        )
    else:
        job.command="""
parseSsearch.pl \\
  --infile {infile} \\
  --rev \\
  > {outfile_raw} && \\
parseMiRNATargets.pl \\
  --infile {outfile_raw} \\
  > {outfile_parsed}""".format(
            infile = infile,
            outfile_raw = outfile_raw,
            outfile_parsed = outfile_parsed
        )

    return job

def merge_fwd_rev_targets(infile_fwd, infile_rev, outfile):

    job = Job(
        [infile_fwd, infile_rev],
        [outfile],
        [
            ['tools', 'module_tools']
        ]
    )
  
  #-n {outfile_normalized} \\
    job.command="""
mergeMiRNATargets.pl \\
  --infile_fwd {infile_fwd} \\
  --infile_rev {infile_rev} \\
  > {outfile}""".format(
        infile_fwd = infile_fwd,
        infile_rev = infile_rev,
        outfile = outfile
    )

    return job

def merge_flagstats(infiles_fs, infiles_qc, outfile):
        
    job = Job(
        infiles_fs + infiles_qc,
        [outfile],
        [
            ['caf_tools', 'module_tools']
        ]
    )

    job.command="""
mergeFlagStatsAndQCMiRNA.pl --infilesFlagstat {infile_flagstat} --infilesQC {infile_qc} > {outfile}""".format(
    infile_flagstat = ",".join(infiles_fs),
    infile_qc = ",".join(infiles_qc),
    outfile = outfile
    )

    return job
