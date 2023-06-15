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

# CAF Modules
from core.config import *
from core.job import *


def bbduk(infile, contam, ncontam, log, db, infile_done=False):

    if(infile_done == False):
        infiles = [infile]
    else:
        infiles = [infile, infile_done]

    #ncontam_gz = ncontam + ".gz"
    #contam_gz = contam + ".gz"

    job = Job(
        infiles, 
        [contam, ncontam, log],
        [
            ['bbmap', 'module_bbmap']
        ]
    )
        
    job.command="""
bbduk.sh \\
  in={infile} \\
  stats={log} \\
  out={ncontam} \\
  outm={contam} \\
  k={k} \\
  minkmerhits={c} \\
  ref={db} \\
  overwrite=true \\
  threads=1""".format(
    infile = infile,
    log = log,
    ncontam = ncontam,
    contam = contam,
    k = config.param('bbduk', 'k', 'int'),
    c = config.param('bbduk', 'c', 'int'),
    db = db
    ) 
    return job

def create_interleaved_fastq(reads1, reads2, tmp, outfile):
    job = Job(
        [reads1, reads2],
        [outfile],
        [
            
            ['caf_tools', 'module_tools'],
            ['pigz', 'module_pigz']
        ]
    )
    
    job.command="""
createInterleavedFastq.pl \\
  --reads1 {reads1} \\
  --reads2 {reads2} \\
  > {tmp} && pigz -p {num_threads} -f {tmp}""".format(
    reads1 = reads1,
    reads2 = reads2,
    tmp = tmp,
    num_threads = config.param("interleaved_fastq", "num_threads", 1, "posint") 
    )

    return job


def trimmomatic(input1, input2, paired_output1, unpaired_output1, paired_output2, unpaired_output2, quality_offset, trim_log, trim_stats):

    job = Job(
        [input1, input2], 
        [paired_output1, unpaired_output1, paired_output2, unpaired_output2, trim_log, trim_stats],
        [
            ['java', 'module_java'], 
            ['trimmomatic', 'module_trimmomatic']
        ]
    )

    threads = config.param('trim', 'threads', type='posint')
    adapter_file = config.param('trim', 'adapter_fasta', type='filepath')
    illumina_clip_settings = config.param('trim', 'illumina_clip_settings')
    trailing_min_quality = config.param('trim', 'trailing_min_quality', type='int')
    min_length = config.param('trim', 'min_length', type='posint')
    headcrop = config.param('trim', 'headcrop', required=False, type='int')
    sliding_window1 = config.param('trim', 'sliding_window1', required=False, type='int')
    sliding_window2 = config.param('trim', 'sliding_window2', required=False, type='int')
   
    if not isinstance( sliding_window1, int ):
        if not(sliding_window1 and sliding_window1.strip()):
            sliding_window1 = 4
    
    if not isinstance( sliding_window2, int ):
        if not(sliding_window2 and sliding_window2.strip()):
            sliding_window2 = 15

    job.command = """
java -XX:ParallelGCThreads={threads} -Xmx2G -jar \$TRIMMOMATIC_JAR {mode} \\
  -threads {threads} \\
  -phred{quality_offset} \\
  {input1} {input2} \\
  {paired_output1} {unpaired_output1} {paired_output2} {unpaired_output2} \\
  ILLUMINACLIP:{adapter_file}{illumina_clip_settings} \\
  TRAILING:{trailing_min_quality} \\
  SLIDINGWINDOW:{sliding_window1}:{sliding_window2} \\
  MINLEN:{min_length} \\
  HEADCROP:{headcrop}""".format(
        mode = "PE",
        threads = threads,
        quality_offset = quality_offset,
        input1 = input1,
        input2 = input2,
        paired_output1 = paired_output1,
        paired_output2 = paired_output2,
        unpaired_output1 = unpaired_output1,
        unpaired_output2 = unpaired_output2,
        adapter_file=adapter_file,
        illumina_clip_settings=illumina_clip_settings,
        trailing_min_quality=trailing_min_quality,
        min_length = min_length,
        sliding_window1 = sliding_window1,
        sliding_window2 = sliding_window2,
        headcrop = str(headcrop)
    )

    job.command += " \\\n  2> " + trim_log

    # Compute statistics
    job.command += " && \\\ngrep ^Input " + trim_log + " | perl -pe 's/^Input Read Pairs: (\\d+).*Both Surviving: (\\d+).*Forward Only Surviving: (\\d+).*$/Raw Fragments,\\1\\nFragment Surviving,\\2\\nSingle Surviving,\\3/' > " + trim_stats

    return job

def templateSub(outdir):
        
    job = Job(
        ["undef"],
        ["undef"],
        [
            ['memtime', 'module_memtime']
        ]
    )
    job.command="memtime"

    return job


def megahit(infile, outdir, type):
     
    if(type == "pe"):
        infile_string = "--12 " + ",".join(infile);
    if(type == "se"):
        infile_string = "-r " + ",".join(infile);

    job = Job(
        infile,
        [os.path.join(outdir, "final.contigs.fa"), os.path.join(outdir, "Contigs.fasta")],
        [
            ['megahit', 'module_megahit']
        ]
    )
    job.command="""
rm {outdir} -rf && \\
megahit -t {num_threads} --k-min {kmin} --k-max {kmax} --k-step {kstep} \\
  --min-contig-len {min_contig_length} \\
  {infile} \\
  --out-dir {outdir} && ln -s -f final.contigs.fa Contigs.fasta && mv Contigs.fasta ./{outdir}/""".format(
    min_contig_length = config.param("megahit", "min_contig_length", 1, "posint"),
    num_threads = config.param("megahit", "num_threads", 1, "posint"),
    kmin = config.param("megahit", "kmin", 1, "posint"),
    kmax = config.param("megahit", "kmax", 1, "posint"),
    kstep = config.param("megahit", "kstep", 1, "int"),
    infile = infile_string,
    outdir = outdir
  )

    return job

def spades(infiles, outdir, type):
     
    if(type == "pe"):
        infiles_string = "--12 " + " --12 ".join(infiles);
    if(type == "se"):
        infiles_string = "-r " + ",".join(infiles);

    job = Job(
        infiles,
        [os.path.join(outdir, "contigs.fasta"), os.path.join(outdir, "Contigs.fasta")],
        [
            ['spades', 'module_spades']
        ]
    )
    job.command="""
rm {outdir} -rf && \\
spades.py --careful -t {num_threads} -k {k} \\
   {infiles} \\
   -o {outdir} && ln -s -f contigs.fasta Contigs.fasta && mv Contigs.fasta ./{outdir}/""".format(
    num_threads = config.param("spades", "num_threads", 1, "posint"),
    k = config.param("spades", "k", 1, "string"),
    infiles = infiles_string,
    outdir = outdir
  )

    return job

def spades_merged(infiles_merged, infiles_R1, infiles_R2, outdir):
    
    if(len(infiles_merged) == 1):
        infiles_merged_string = " --merged " + infiles_merged[0];
    else:
        infiles_merged_string = " --merged " + " --merged ".join(infiles_merged);
    
    if(len(infiles_R1) == 1):
        infiles_R1_string = " -1 " + infiles_R1[0];
    else:
        infiles_R1_string = " -1 " + " -1 ".join(infiles_R1);
    
    if(len(infiles_R2) == 1):
        infiles_R2_string = " -2 " + infiles_R2[0];
    else:
        infiles_R2_string = " -2 " + " -2 ".join(infiles_R2);
    
    infiles = infiles_merged + infiles_R1 + infiles_R2
    #sys.stderr.write("files in spades_merged: " + str(infiles) + "\n")   

    job = Job(
        infiles,
        [os.path.join(outdir, "contigs.fasta"), os.path.join(outdir, "Contigs.fasta")],
        [
            ['spades', 'module_spades']
        ]
    )
    job.command="""
rm {outdir} -rf && \\
spades.py --careful -t {num_threads} -k {k} \\
   {infiles_merged} \\
   {infiles_R1} \\
   {infiles_R2} \\
   -o {outdir} && ln -s -f contigs.fasta Contigs.fasta && mv Contigs.fasta ./{outdir}/""".format(
    num_threads = config.param("spades", "num_threads", 1, "posint"),
    k = config.param("spades", "k", 1, "string"),
    infiles_merged = infiles_merged_string,
    infiles_R1 = infiles_R1_string,
    infiles_R2 = infiles_R2_string,
    outdir = outdir
  )

    return job

def remove_unpaired_reads(infile, outfile_paired, unpairedR1, unpairedR2):
    job = Job(
        [infile],
        [outfile_paired + ".gz"],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]                 
    )    
    job.command="""                    
removeUnpaired.pl \\
  --infile {infile} \\
  --outfile_paired {outfile_paired} \\
  --outfile_1 {unpairedR1} \\
  --outfile_2 {unpairedR2} \\
  --num_threads {num_threads} && \
gzip -f {outfile_paired}""".format(
        infile = infile,
        outfile_paired = outfile_paired,
        unpairedR1 = unpairedR1,
        unpairedR2 = unpairedR2,
        num_threads =  config.param('remove_unpaired', 'num_threads', 'int')
    )
                        
    return job

def compile_assembly_stats(infile, outfile):
        
    job = Job(
        [infile],
        [outfile],
        [
            
            ['caf_tools', 'module_tools']
        ]
    )
    job.command="""
compileAssemblyResultsSingle.pl \\
  --infile {infile} > {outfile}""".format(
    infile = infile,
    outfile = outfile
     )

    return job

def flash(infile, outdir, prefix):
    #sys.stderr.write("file in flash: " + os.path.join(outdir, prefix + ".extendedFrags.fastq") + "\n")   
    job = Job(
        [infile],
        [os.path.join(outdir, prefix + ".extendedFrags.fastq"), 
         os.path.join(outdir, prefix + ".notCombined_1.fastq"),
         os.path.join(outdir, prefix + ".notCombined_2.fastq")],
        [
            
            ['flash', 'module_flash']
        ]
    )
    job.command="""
flash -d {outdir} -M 235 -o {prefix} \\
  --interleaved-input {infile}""".format(
    infile = infile,
    outdir = outdir,
    prefix = prefix
     )

    return job

def bwa_mem_samtools(reference, bwt, infile1, infile2, outfile):
    job = Job(
        [infile1, infile2, reference, bwt],
        [outfile],
        [
            
            ['samtools', 'module_samtools'],
            ['bwa', 'module_bwa']
        ]
    )
    job.command="""
bwa mem -M \\
  -t {num_threads} \\
  {reference} \\
  {infile1} \\
  {infile2} \\
  | samtools view -Sbh -F 0x100 - > {outfile}.tmp && \\
  samtools sort -@ {num_threads} -m {mem_per_thread} {outfile}.tmp -o {outfile}.tmp.sorted.bam && \\
  mv {outfile}.tmp.sorted.bam {outfile} && \\
  rm {outfile}.tmp && \\
  samtools index {outfile}""".format(
    num_threads = config.param('bwa', 'num_threads', type='int', required=True), 
    mem_per_thread = config.param('samtools', 'mem_per_thread', required=True),
    infile1 = infile1,
    infile2 = infile2,
    reference = reference,
    outfile = outfile
    )

    return job

def flagstats(infile, outfile):
        
    job = Job(
        [infile],
        [outfile],
        [
            
            ['samtools', 'module_samtools'],
            ['caf_tools', 'module_tools']
        ]
    )
    job.command="""
samtools flagstat {infile} > {outfile}""".format(
    infile = infile,
    outfile = outfile
    )
    return job


def split_pairs(infile, outfile1, outfile2):
        
    job = Job(
        [infile],
        [outfile1, outfile2],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
    job.command="""
splitPairsGz.pl \\
  --infile {infile} \\
  --outfile_1 {outfile1} \\
  --outfile_2 {outfile2}""".format(
    infile = infile,
    outfile1 = outfile1,
    outfile2 = outfile2
    )

    return job

def make_index(fasta, bwt):
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
bwa index {fasta}""".format(
    fasta = fasta,
    bwt = bwt
    )

    return job


def fasta_to_bed(fasta, bed):
        
    job = Job(
        [fasta],
        [bed],
        [
            
            ['samtools', 'module_samtools'],
            ['caf_tools', 'module_tools'],
            ['bwa', 'module_bwa']
        ]
    )
    job.command="""
fastaToBed.pl --fasta {fasta} > {bed}""".format(
    fasta = fasta,
    bed = bed
    )

    return job


def coverage_bed(infile, bed, outfile, flag):
    job = Job(
        [infile, bed],
        [outfile],
        [
            
            ['bedtools', 'module_bedtools'],
            ['samtools', 'module_samtools']
        ]
    )

    job.command="""
samtools view -b -f {flag} {infile} | \\
  coverageBed -abam stdin \\
  -b {bed} \\
  -counts \\
  > {outfile}""".format(
    infile = infile,
    bed = bed,
    outfile = outfile,
    flag = flag
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
mergeFlagStatsAndQC.pl --infilesFlagstat {infile_flagstat} --infilesQC {infile_qc} > {outfile}""".format(
    infile_flagstat = ",".join(infiles_fs),
    infile_qc = ",".join(infiles_qc),
    outfile = outfile
    )

    return job


def get_insert_size(indir, outfile, flagstats): #flagstats just here for dependency purpose.
    
    job = Job(
        [flagstats],
        [outfile],
        [
            
            ['caf_tools', 'module_tools'],
            ['samtools', 'module_samtools']
        ]
    )
    
    job.command="""
prepareBamsForGamngs.pl --indir {indir} > {outfile}""".format(
    indir = indir,
    outfile = outfile
    )

    return job

def merge_bams(bams, outfile):
    
    bams_string = " ".join(bams)

    job = Job(
        bams,
        [outfile],
        [
            ['samtools', 'module_samtools']
        ]
    )
    
    job.command="""
samtools merge {outfile} \\
  {bams_string}""".format(
    bams_string = bams_string,
    outfile = outfile
    )

    return job


def quast(contigs, bams, outdir):
    
    infiles = contigs + bams

    contigs_string = " ".join(contigs) 
    bams_string = ",".join(bams)

    job = Job(
        infiles,
        [os.path.join(outdir, "report.tsv")],
        [
            ['python3', 'module_python3'],
            ['quast', 'module_quast']
        ]
    )
    
    job.command="""
quast.py {contigs_string} \\
  --bam {bams_string} \\
  -o {outdir}""".format(
    contigs_string = contigs_string,
    bams_string = bams_string,
    outdir = outdir
    )

    return job

def convert_orf_ids_for_cat(infile_gff, infile_blastp, outfile_blastp):
    
    job = Job(
        [infile_gff, infile_blastp],
        [outfile_blastp],
        [
            ['tools', 'module_tools']
        ]
    )
    job.command="""
convertDiamondBlastpORFIDForCAT.pl \\
  --infile_gff {infile_gff} \\
  --infile_blastp {infile_blastp} \\
  > {outfile_blastp}""".format(
    infile_gff = infile_gff,
    infile_blastp = infile_blastp,
    outfile_blastp = outfile_blastp
    )

    return job

def CAT(infile_fna, infile_faa, infile_blastp, prefix):
    
    job = Job(
        [infile_fna, infile_faa, infile_blastp],
        [
            os.path.join(prefix + ".ORF2LCA.txt"),
            os.path.join(prefix + ".log"),
            os.path.join(prefix + ".contig2classification.txt"),
            os.path.join(prefix + ".contig2classification_with_names.tsv")
        ],
        [
            ['CAT', 'module_CAT'],
            ['diamond', 'module_diamond'],
            ['python3', 'module_python3'],
            ['prodigal', 'module_prodigal']
        ]
    )
    #./CAT_pack/CAT contigs -c bin1285/out.1285.fna -p bin1285/out.1285.faa 
    #-d ./2019-10-29_CAT_database -t 2019-10-29_taxonomy 
    #-a bin1285/out.1285_renamed_diamondblastpnr_altids.tsv -o out.1285_ext --force
    job.command="""
CAT contigs \\
 -c {infile_fna} \\
 -p {infile_faa} \\
 -a {infile_blastp} \\
 -d {database_folder} \\
 -t {taxonomy_folder} \\
 -o {prefix} --force && \\
CAT add_names \\
 -i {infile_names} \\
 -o {outfile_names} \\
 -t {taxonomy_folder} --force""".format(
    infile_fna = infile_fna,
    infile_faa = infile_faa,
    infile_blastp = infile_blastp,
    prefix = prefix,
    infile_names = os.path.join(prefix + ".contig2classification.txt"),
    outfile_names = os.path.join(prefix + ".contig2classification_with_names.tsv"),
    database_folder = config.param("CAT", "database_folder", 1, "string"),
    taxonomy_folder = config.param("CAT", "taxonomy_folder", 1, "string")
    )

    return job

def generate_otu_table_consensus(infile_taxonomy, infile_abundance, outfile_taxonomy, outfile):

    job = Job(
        [infile_taxonomy, infile_abundance],
        [outfile], 
        [
            ['caf_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
    
    job.command="""
generateOTUTableFromCAT.pl \\
  --infile_taxonomy {infile_taxonomy} \\
  --infile_abundance {infile_abundance} \\
  --outfile_taxonomy {outfile_taxonomy} \\
  > {outfile}""".format(
    infile_taxonomy = infile_taxonomy,
    infile_abundance = infile_abundance,
    outfile_taxonomy = outfile_taxonomy,
    outfile = outfile
    )

    return job

def generate_gff(infile_gff, infile_fasta, kegg, pfam, cog, kog, 
                 taxonomy, ublast_nr, 
                 outfile_gff, outfile_fasta, outfile_annotations):

    job = Job(
        [infile_gff, infile_fasta, kegg, pfam, cog, kog, taxonomy, ublast_nr], 
        [outfile_gff, outfile_fasta, outfile_annotations],
        [
            ['caf_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
generateGFF.pl \\
  --infile_gff {infile_gff} \\
  --infile_fasta {infile_fasta} \\
  --prefix \\"{prefix}\\" \\
  --pfam {pfam} \\
  --cog {cog} \\
  --kog {kog} \\
  --kegg {kegg} \\
  --taxonomy {taxonomy} \\
  --ublast {ublast_nr} \\
  --outfile_gff {outfile_gff} \\
  --outfile_fasta {outfile_fasta} \\
  > {outfile_annotations}""".format(
    infile_gff = infile_gff,
    infile_fasta = infile_fasta,
    prefix = config.param('DEFAULT', 'supp_info', 1, 'string'),
    pfam = pfam,
    cog = cog,
    kog = kog,
    kegg = kegg,
    taxonomy = taxonomy,
    ublast_nr = ublast_nr,
    outfile_gff = outfile_gff,
    outfile_fasta = outfile_fasta,
    outfile_annotations = outfile_annotations
    )

    return job


def prodigal(infile, outfile_gff, outfile_fna, outfile_faa, renamed_gff, renamed_fna, renamed_faa):
    
    job = Job(
        [infile],
        [outfile_gff, outfile_fna, outfile_faa, renamed_fna, renamed_gff, renamed_faa],
        [
            
            ['prodigal', 'module_prodigal'],
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )
    
    job.command="""
prodigal -i {infile} -f gff -p single \\
  -o {outfile_gff} \\
  -a {outfile_faa} \\
  -d {outfile_fna} && \\
convertProdigalNames.pl \\
  --gff {outfile_gff} \\
  --fna {outfile_fna} \\
  --faa {outfile_faa} \\
  --renamed_gff {renamed_gff} \\
  --renamed_faa {renamed_faa} \\
  > {renamed_fna}""".format(
        infile = infile,
        outfile_gff = outfile_gff,
        outfile_fna = outfile_fna,
        outfile_faa = outfile_faa,
        renamed_fna = renamed_fna,
        renamed_gff = renamed_gff,
        renamed_faa = renamed_faa
    )

    return job

def generate_taxonomy_table_consensus(infile_taxonomy, outfile_taxonomy):

    job = Job(
        [infile_taxonomy],
        [outfile_taxonomy], 
        [
            ['caf_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
    
    job.command="""
generateTaxonomyFromCAT.pl \\
  --infile_taxonomy {infile_taxonomy} \\
  > {outfile_taxonomy}""".format(
    infile_taxonomy = infile_taxonomy,
    outfile_taxonomy = outfile_taxonomy
    )

    return job

def generate_external_genomes_file(infile_done, external_genomes, outfile):

    job = Job(
        [infile_done],
        [outfile], 
        [
            ['perl', 'module_perl']
        ]
    )
    
    job.command="""
echo '{external_genome_line}' \\
  > {outfile}""".format(
    infile_done = infile_done,
    external_genome_line = "\n".join(external_genomes),
    outfile = outfile
    )

    return job

def anvio_gen_genome_storage(external_genomes_file, outfile):

    job = Job(
        [external_genomes_file], 
        [outfile],
        [
            ['anvio', 'module_anvio'],
            ['prodigal', 'module_prodigal'],
            ['hmmer', 'module_hmmer'],
            ['blast', 'module_blast'],
            ['bowtie2', 'module_bowtie2'],
            ['bwa', 'module_bwa'],
            ['diamond', 'module_diamond'],
            ['fastani', 'module_fastani'],
            ['fasttree', 'module_fasttree'],
            ['muscle', 'module_muscle'],
            ['samtools', 'module_samtools'],
            ['R', 'module_R'],
            ['mcl', 'module_mcl']
        ]
    )
    
    job.command="""
rm -rf {outfile} && \\
bash -c 'set +u && source \$ANVIO_VENV/activate && set -u && \\
    anvi-gen-genomes-storage \\
        -e {external_genomes_file} \\
        -o {outfile}'""".format(
    outfile = outfile,
    external_genomes_file = external_genomes_file
    )

    return job

def anvio_pangenome(genomes_db, outdir, done_file):

    job = Job(
        [genomes_db], 
        [done_file],
        [
            ['anvio', 'module_anvio'],
            ['prodigal', 'module_prodigal'],
            ['hmmer', 'module_hmmer'],
            ['blast', 'module_blast'],
            ['bowtie2', 'module_bowtie2'],
            ['bwa', 'module_bwa'],
            ['diamond', 'module_diamond'],
            ['fastani', 'module_fastani'],
            ['fasttree', 'module_fasttree'],
            ['muscle', 'module_muscle'],
            ['samtools', 'module_samtools'],
            ['R', 'module_R'],
            ['mcl', 'module_mcl']
        ]
    )
    
    job.command="""
rm -rf {outdir}/* && \\
bash -c 'set +u && source \$ANVIO_VENV/activate && set -u && \\
    anvi-pan-genome -W -g {genomes_db} \\
        --project-name 'Pangenome' \\
        --output-dir {outdir} \\
        --num-threads {num_threads} \\
        --minbit {min_bit} \\
        --mcl-inflation {mcl_inflation} \\
        --enforce-hierarchical-clustering' && \\
touch {done_file}""".format(
    genomes_db = genomes_db,
    outdir = outdir,
    done_file = done_file,
    num_threads = config.param('anvio_pangenome', 'num_threads', 1, 'posint'),
    min_bit = config.param('anvio_pangenome', 'min_bit', 1, 'float'),
    mcl_inflation = config.param('anvio_pangenome', 'mcl_inflation', 1, 'posint')
    )

    return job

def anvio_compute_similarity(infile_done, indir, external_genomes_file, outdir, done_file):

    genomes_db = os.path.join(indir, "Pangenome-PAN.db")

    job = Job(
        [infile_done], 
        [done_file],
        [
            ['anvio', 'module_anvio'],
            ['python3', 'module_python3'],
            ['fastani', 'module_fastani']
        ]
    )

 
    job.command="""
rm -rf {outdir}/* && \\
bash -c 'set +u && source \$ANVIO_VENV/activate && set -u && \\
    anvi-compute-genome-similarity \\
        -e {external_genomes_file} \\
        --program fastANI \\
        --output-dir {outdir} \\
        --num-threads {num_threads} \\
        --pan-db {genomes_db}' \\
touch {done_file}""".format(
    genomes_db = genomes_db,
    outdir = outdir,
    done_file = done_file,
    external_genomes_file = external_genomes_file,
    num_threads = config.param('anvio_pangenome', 'num_threads', 1, 'posint')
    )

    return job

def merge_feature_tables(prefixes, infiles, outfile):

    job = Job(
        infiles,
        [outfile], 
        [
            ['R', 'module_R'],
            ['tools', 'module_tools']
        ]
    )
    
    job.command="""
mergeFeatureTablesIntoSingleFT.R \\
 -i {infiles} \\
 -p {prefixes} \\
 -o {outfile}""".format(
    infiles = ",".join(infiles),
    prefixes = ",".join(prefixes),
    outfile = outfile
    )

    return job

def merge_fasta(prefixes, infiles, outfile):

    job = Job(
        infiles,
        [outfile], 
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )
    
    job.command="""
mergeBarrnapSequences.pl \\
 --prefixes {prefixes} \\
 --infiles {infiles} > {outfile}""".format(
    infiles = ",".join(infiles),
    prefixes = ",".join(prefixes),
    outfile = outfile
    )

    return job
