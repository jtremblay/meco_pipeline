#!/usr/bin/env python

#LICENSE AND COPYRIGHT

#Copyright (C) 2025 Julien Tremblay

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

# MECO Modules
from core.config import *
from core.job import *

# Added by Patrick.G new monitoring job which is designed to send email when all jobs are finished
# Since it depend on every jobs in step, it will send an INVALID_DEPEND mail is something failed
# an END mail will be sented if everything went correctly
# The stepname argument is a placeholder if I want to proceed with command override in scheduler
# shotgunmg.py will have to be modified to include the step name for every occurence of shotgun_metagenomics.monitoring()
def monitoring(stepname=""):
    job = Job(
        [],
        [],
        [],
    )
    # Since we leveraging Slurm mail sending system, no command are needed here
    # Anyway, is necessary, scheduler will override the command to make it send an email on run
    job.command = """sleep 1"""
    return job

# ADDED BY Patrick Gagne on october 28nd 2024 (improvement to first shotmg steps to obtain quality data)
def fastqc_qa_pe(input1,input2,outdir, output_R1, output_R2):
    threads = config.param('fastqc', 'threads', type='posint')
    job = Job(
        [input1,input2],
        [output_R1,output_R2],
        [
            ['fastqc','module_fastqc']
        ]
    )
    job.command = """
fastqc {input1} {input2} \\
  -t {threads} \\
  -o {outdir}""".format(
        threads = threads,
        input1 = input1,
        input2 = input2,
        outdir = outdir
    )
    return job

# ADDED BY Patrick Gagne on october 28nd 2024 (improvement to first shotmg steps to obtain quality data)
def fastqc_qa_se(input1,outdir,output_f):
    threads = config.param('fastqc', 'threads', type='posint')
    job = Job(
        [input1],
        [output_f],
        [
            ['fastqc','module_fastqc']
        ]
    )
    job.command = """
fastqc {input1} {input2} \\
  -t {threads} \\
  -o {outdir}""".format(
        threads = threads,
        input1 = input1,
        input2 = input2,
        outdir = outdir
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
    crop = config.param('trim', 'crop', required=False, type='int')
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
    
    if isinstance( crop, int ):
        job.command += """ CROP:{crop}""".format(crop = str(crop))

    job.command += " \\\n  2> " + trim_log
    

    # Compute statistics
    job.command += " && \\\ngrep ^Input " + trim_log + " | perl -pe 's/^Input Read Pairs: (\\\d+).*Both Surviving: (\\\d+).*Forward Only Surviving: (\\\d+).*$/Raw Fragments,\\\\1\\\\nFragment Surviving,\\\\2\\\\nSingle Surviving,\\\\3/' > " + trim_stats

    return job

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
        average_quality = average_quality
    )

    job.command += " \\\n  2> " + trim_log

    # Compute statistics
    job.command += " && \\\ngrep ^Input " + trim_log + " | perl -pe 's/^Input Read Pairs: (\\d+).*Both Surviving: (\\d+).*Forward Only Surviving: (\\d+).*$/Raw Fragments,\\1\\nFragment Surviving,\\2\\nSingle Surviving,\\3/' > " + trim_stats

    return job

def bbduk(infile, contam, ncontam, log, db, infile_done=False):

    if(infile_done == False):
        infiles = [infile]
    else:
        infiles = [infile, infile_done]

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
  threads=1 {sf}""".format(
    infile = infile,
    log = log,
    ncontam = ncontam,
    contam = contam,
    k = config.param('bbduk', 'k', 'int'),
    c = config.param('bbduk', 'c', 'int'),
    sf = config.param('bbduk','secondary_flags'),
    db = db
    ) 
    return job

def bbduk_paired(infile_R1, infile_R2, 
                 contam_R1, contam_R2, 
                 ncontam_R1, ncontam_R2, 
                 log, db):

    job = Job(
        [infile_R1, infile_R2], 
        [contam_R1, contam_R2, ncontam_R1, ncontam_R2, log],
        [
            ['bbmap', 'module_bbmap']
        ]
    )
        
    job.command="""
bbduk.sh \\
  in={infile_R1} \\
  in2={infile_R2} \\
  out={ncontam_R1} \\
  out2={ncontam_R2} \\
  outm={contam_R1} \\
  outm2={contam_R2} \\
  stats={log} \\
  k={k} \\
  minkmerhits={c} \\
  ref={db} \\
  overwrite=true \\
  threads=1 {sf}""".format(
    infile_R1 = infile_R1,
    infile_R2 = infile_R2,
    ncontam_R1 = ncontam_R1,
    ncontam_R2 = ncontam_R2,
    contam_R1 = contam_R1,
    contam_R2 = contam_R2,
    log = log,
    k = config.param('bbduk', 'k', 'int'),
    c = config.param('bbduk', 'c', 'int'),
    sf = config.param('bbduk','secondary_flags'),
    db = db
    ) 
    return job

def subsample(infile_R1, infile_R2, outfile_R1, outfile_R2):

    job = Job(
        [infile_R1, infile_R2], 
        [outfile_R1, outfile_R2],
        [
            ['bbmap', 'module_bbmap']
        ]
    )
        
    job.command="""
reformat.sh \\
  in={infile_R1} \\
  in2={infile_R2} \\
  out={outfile_R1} \\
  out2={outfile_R2} \\
  samplereadstarget={samplereadstarget} \\
  sampleseed={sampleseed} \\
  overwrite=true \\
  threads={num_threads}""".format(
    infile_R1 = infile_R1,
    infile_R2 = infile_R2,
    outfile_R1 = outfile_R1,
    outfile_R2 = outfile_R2,
    samplereadstarget = config.param('subsample', 'samplereadstarget', type='int', required=True),
    sampleseed = config.param('subsample', 'sampleseed', type='int', required=True),
    num_threads = config.param('subsample', 'num_threads', type='int', required=True)
    )

    return job

def duk_gz_matched_only(infile, outfile, log, db):
               

    job = Job(
        [infile], 
        [outfile, log],
        [
            
            ['meco_tools', 'module_tools'],
            ['perl', 'module_perl'],
            ['duk', 'module_duk']
        ]
    )
        
    job.command=""" 
gunzip -c {infile} | duk \\
-o {log} \\
-n /dev/null \\
-m {outfile} \\
-k {k} \\
-s {s} \\
-c {c} \\
{db}""".format(
    infile = infile,
    log = log,
    outfile = outfile,
    k = config.param('duk', 'k', 'int'),
    s = config.param('duk', 's', 'int'),
    c = config.param('duk', 'c', 'int'),
    db = db
    ) 
    return job

def join_reads(fastq1, fastq2, outfile):

    job = Job(
        [fastq1, fastq2], 
        [outfile],
        [
            
            ['meco_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
         
    job.command="""
joinReads.pl \\
  --infile_R1 {fastq1} \\
  --infile_R2 {fastq2} \\
  --outfile {outfile}""".format(
    fastq1 = fastq1,
    fastq2 = fastq2,
    outfile = outfile
    )

    return job

def merge_duk_logs_interleaved(logs, readset_ids, outfile):
   
    job = Job(
        logs,
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
mergeDukLogs.pl --logs {logs} --ids {readset_ids} > {outfile}""".format(
    logs = ",".join(logs),
    readset_ids = ",".join(readset_ids),
    outfile = outfile
    )
    
    return job

def merge_duk_unmapped_logs_interleaved(logs, readset_ids, outfile):
   
    job = Job(
        logs,
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
mergeDukLogsSub.pl --logs {logs} --ids {readset_ids} > {outfile}""".format(
    logs = ",".join(logs),
    readset_ids = ",".join(readset_ids),
    outfile = outfile
    )
    
    return job

def create_interleaved_fastq(reads1, reads2, tmp, outfile):
    job = Job(
        [reads1, reads2],
        [outfile],
        [
            
            ['meco_tools', 'module_tools'],
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

def remove_unpaired_reads_and_split(infile, unpaired_reads1, unpaired_reads2, paired_reads1, paired_reads2, outfile):
    
    tmp1 = os.path.splitext(unpaired_reads1)[0]
    tmp2 = os.path.splitext(unpaired_reads2)[0]
    tmp12 = os.path.splitext(outfile)[0]
    
    job = Job(
        [infile],
        [outfile, unpaired_reads1, unpaired_reads2, paired_reads1, paired_reads2],
        [
            
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
removeUnpairedReads.pl \\
  --infile {infile} \\
  --unpaired_reads1 {tmp1} \\
  --unpaired_reads2 {tmp2} \\
  > {tmp12} && \\
gzip -f {tmp1} && gzip -f {tmp2} && gzip -f {tmp12} && \\
splitPairsGz.pl \\
  --infile {tmp12}.gz \\
  --outfile_1 {paired_reads1} \\
  --outfile_2 {paired_reads2}""".format(
    tmp1 = tmp1,
    tmp2 = tmp2,
    tmp12 = tmp12,
    infile = infile, 
    paired_reads1 = paired_reads1,
    paired_reads2 = paired_reads2
    )

    return job

def megahit(infile, outdir, type):
     
    if(type == "pe"):
        infile_string = "--12 " + ",".join(infile);
    if(type == "se"):
        infile_string = "-r " + ",".join(infile);

    job = Job(
        infile,
        [os.path.join(outdir, "Contigs.fasta"), os.path.join(outdir, "final.contigs.fa")],
        [
            ['megahit', 'module_megahit']
        ]
    )
    job.command="""
rm {outdir} -rf && \\
megahit -t {num_threads} --k-min {kmin} --k-max {kmax} --k-step {kstep} \\
  --min-contig-len {min_contig_length} \\
  {infile} \\
  --out-dir {outdir} && ln -s -f final.contigs.fa Contigs.fasta && mv Contigs.fasta ./assembly/""".format(
    min_contig_length = config.param("megahit", "min_contig_length", 1, "posint"),
    num_threads = config.param("megahit", "num_threads", 1, "posint"),
    kmin = config.param("megahit", "kmin", 1, "posint"),
    kmax = config.param("megahit", "kmax", 1, "posint"),
    kstep = config.param("megahit", "kstep", 1, "int"),
    infile = infile_string,
    outdir = outdir
  )

    return job

def megahit_R1_R2(infiles_R1, infiles_R2, outdir, type):
    
    memory = config.param('megahit', 'memory', required=False, type='float')
    requested_memory = config.param('megahit', 'requested_memory', required=False, type='int')

    if(type == "pe"):
        #infiles_string_R1 = "-1 " + ",".join(infiles_R1);
        #infiles_string_R2 = "-2 " + ",".join(infiles_R2);
        
        infiles_string_R1 = "-1 " + " -1 ".join(infiles_R1);
        infiles_string_R2 = "-2 " + " -2 ".join(infiles_R2);

    job = Job(
        infiles_R1 + infiles_R2,
        [os.path.join(outdir, "Contigs.fasta"), os.path.join(outdir, "final.contigs.fa")],
        [
            ['megahit', 'module_megahit']
        ]
    )
    job.command="""
rm {outdir} -rf && \\
megahit -t {num_threads} --k-min {kmin} --k-max {kmax} --k-step {kstep} \\
  --min-contig-len {min_contig_length} \\
  {infiles_R1} \\
  {infiles_R2} \\""".format(
    min_contig_length = config.param("megahit", "min_contig_length", 1, "posint"),
    num_threads = config.param("megahit", "num_threads", 1, "posint"),
    kmin = config.param("megahit", "kmin", 1, "posint"),
    kmax = config.param("megahit", "kmax", 1, "posint"),
    kstep = config.param("megahit", "kstep", 1, "int"),
    infiles_R1 = infiles_string_R1,
    infiles_R2 = infiles_string_R2,
    outdir = outdir
  )
    
    if isinstance( memory, float ) and isinstance( requested_memory, int):
        requested_memory_in_bytes = requested_memory * 1000000
        memory_in_bytes = requested_memory_in_bytes * memory
        memory_in_bytes = round(memory_in_bytes, ndigits=0)
        job.command += """
  --memory {memory_in_bytes} \\
  --out-dir {outdir} && \\""".format(memory_in_bytes = str(int(memory_in_bytes)), outdir = outdir)
    else:
        job.command += """
  --out-dir {outdir} && \\""".format(outdir = outdir)

    job.command+="""
sed -i 's/^>\\\(\\\S\\\+\\\) .*/>\\\\1/' {outdir}/final.contigs.fa && \\
ln -s -f final.contigs.fa Contigs.fasta && mv Contigs.fasta {outdir}/""".format(outdir = outdir)

    return job

def megahit_R1_R2_merged(infiles_R1, infiles_R2, infiles_merged, outdir, type):
    
    memory = config.param('megahit', 'memory', required=False, type='float')
    requested_memory = config.param('megahit', 'requested_memory', required=False, type='int')

    if(type == "pe"):
        infiles_string_R1 = "-1 " + ",".join(infiles_R1);
        infiles_string_R2 = "-2 " + ",".join(infiles_R2);
        infiles_string_merged = "-r " + ",".join(infiles_merged);

    job = Job(
        infiles_R1 + infiles_R2 + infiles_merged,
        [os.path.join(outdir, "Contigs.fasta"), os.path.join(outdir, "final.contigs.fa")],
        [
            ['megahit', 'module_megahit']
        ]
    )
    job.command="""
rm {outdir} -rf && \\
megahit -t {num_threads} --k-min {kmin} --k-max {kmax} --k-step {kstep} \\
  --min-contig-len {min_contig_length} \\
  {infiles_R1} \\
  {infiles_R2} \\
  {infiles_merged} \\""".format(
    min_contig_length = config.param("megahit", "min_contig_length", 1, "posint"),
    num_threads = config.param("megahit", "num_threads", 1, "posint"),
    kmin = config.param("megahit", "kmin", 1, "posint"),
    kmax = config.param("megahit", "kmax", 1, "posint"),
    kstep = config.param("megahit", "kstep", 1, "int"),
    infiles_R1 = infiles_string_R1,
    infiles_R2 = infiles_string_R2,
    infiles_merged = infiles_string_merged,
    outdir = outdir
  )
    
    if isinstance( memory, float ) and isinstance( requested_memory, int):
        requested_memory_in_bytes = requested_memory * 1000000
        memory_in_bytes = requested_memory_in_bytes * memory
        memory_in_bytes = round(memory_in_bytes, ndigits=0)
        job.command += """
  --memory {memory_in_bytes} \\
  --out-dir {outdir} && \\""".format(memory_in_bytes = str(int(memory_in_bytes)), outdir = outdir)
    else:
        job.command += """
  --out-dir {outdir} && \\""".format(outdir = outdir)

    job.command+="""
sed -i 's/^>\\\(\\\S\\\+\\\) .*/>\\\\1/' {outdir}/final.contigs.fa && \\
ln -s -f final.contigs.fa Contigs.fasta && mv Contigs.fasta ./assembly/""".format(outdir = outdir)

    return job

def spades_R1_R2(infiles_R1, infiles_R2, outdir):
    
    #temporary fix; change fna to fasta
    #trusted_contigs = os.path.splitext(trusted_contigs)[0]
    #trusted_contigs = trusted_contigs + ".fa"
    min_contig_length = str(config.param("spades", "min_contig_length", 1, "posint"))

    infiles_string_R1 = ""
    infiles_string_R2 = ""
    for i in range(len(infiles_R1)):
        #infiles_string_R1 += " --pe" + str(i+1) + "-1 " + infiles_R1[i] 
        #infiles_string_R2 += " --pe" + str(i+1) + "-2 " + infiles_R2[i] 
        infiles_string_R1 += " -1 " + infiles_R1[i]
        infiles_string_R2 += " -2 " + infiles_R2[i]

    job = Job(
        infiles_R1 + infiles_R2,
        [
            os.path.join(outdir, "scaffolds.fasta"), 
            os.path.join(outdir, "scaffolds_gt" + min_contig_length+ ".fasta"),
            os.path.join(outdir, "Contigs.fasta")
        ],
        [
            ['spades', 'module_spades'],
            ['tools','module_tools']
        ]
    )
  
    job.command="""
rm {outdir} -rf && \\
spades.py --meta -t {num_threads} -k {kmers} \\
  {infiles_string_R1} \\
  {infiles_string_R2} \\
  -o {outdir} -m {memory} && \\
filterFastaByLength.pl \\
  --infile {outdir}/scaffolds.fasta \\
  --length {length} > \\
  {outdir}/scaffolds_gt{length}.fasta && \\
ln -s -f scaffolds_gt{length}.fasta Contigs.fasta && mv Contigs.fasta ./assembly/""".format(
    num_threads = config.param("spades", "num_threads", 1, "posint"),
    kmers = config.param("spades", "kmers", 1, "string"),
    length = config.param("spades", "min_contig_length", 1, "posint"),
    memory = config.param("spades", "memory", 1, "posint"),
    infiles_string_R1 = infiles_string_R1,
    infiles_string_R2 = infiles_string_R2,
    outdir = outdir
  )

    return job

def porechop(infile, outfile):
     
    job = Job(
        [infile],
        [outfile],
        [
            ['python3', 'module_python3']
        ]
    )

    job.command="""
porechop -t {num_threads} -i {infile} -o {outfile}""".format(
    infile = infile,
    outfile = outfile,
    num_threads = config.param("porechop", "num_threads", 1, "posint")
  )

    return job


def cat_fastqs(infiles, outfile):
     
    job = Job(
        infiles,
        [outfile],
        [
            #['miniasm', 'module_miniasm']
        ]
    )

    job.command="""
cat {infiles} > {outfile}""".format(
    infiles = " ".join(infiles),
    outfile = outfile
  )

    return job

def minimap2_all_vs_all(infile, outfile):
     
    job = Job(
        [infile],
        [outfile],
        [
            ['minimap2', 'module_minimap2']
        ]
    )

    job.command="""
minimap2 -x {type} \\
  -t {num_threads} \\
  {infile} \\
  {infile} \\
  | gzip -1 > {outfile}
""".format(
    num_threads = config.param("minimap2_all_vs_all", "num_threads", 1, "posint"),
    type = config.param("minimap2_all_vs_all", "type", 1, "string"),
    outfile = outfile,
    infile = infile
  )

    return job

def miniasm(infile_fastq, infile_alignment, outfile_gfa, outfile_fasta):
     
    job = Job(
        [infile_fastq, infile_alignment],
        [outfile_gfa, outfile_fasta],
        [
            ['miniasm', 'module_miniasm']
        ]
    )
    job.command="""
miniasm -f \\
  {infile_fastq} \\
  {infile_alignment} \\
  > {outfile_gfa} && \\
awk '/^S/{{print \\\">\\\"\$2\\\"\\\\n\\\"\$3}}' {outfile_gfa} > {outfile_fasta}""".format(
    infile_fastq = infile_fastq,
    infile_alignment = infile_alignment,
    outfile_gfa = outfile_gfa,
    outfile_fasta = outfile_fasta
  )

    return job

def canu(infiles, outdir):
     
    job = Job(
        infiles,
        [os.path.join(outdir, "Contigs.fasta")],
        [
            ['canu', 'module_canu'],
            ['java', 'module_java']
        ]
    )

    job.command="""
rm {outdir} -rf && \\
canu -useGrid=false maxMemory={max_memory} maxThreads={num_threads} \\
  -p shotgun_MG -d {outdir}  \\
  corMinCoverage=0 \\
  corOutCoverage=all \\
  corMhapSensitivity=high \\
  correctedErrorRate=0.105 \\
  genomeSize=5m \\
  corMaxEvidenceCoverageLocal=10 \\
  corMaxEvidenceCoverageGlobal=10 \\
  -nanopore-raw {infiles}""".format(
    num_threads = config.param("canu", "num_threads", 1, "posint"),
    max_memory = config.param("canu", "max_memory", 1, "posint"),
    kstep = config.param("megahit", "kstep", 1, "int"),
    infiles = " ".join(infiles),
    outdir = outdir
  )

    return job

def ropebwt2(infile, outfile_bwt, outfile_npy):
        
    job = Job(
        [infile],
        [outfile_bwt, outfile_npy],
        [
            ['fmlrc', 'module_fmlrc'],
            ['ropebwt2', 'module_ropebwt2']
        ]
    )
    job.command="""
ropebwt2 \\
  -o {outfile_bwt} \\
  {infile} && \\
fmlrc-convert -i {outfile_bwt} {outfile_npy}""".format(
    infile = infile,
    outfile_bwt = outfile_bwt,
    outfile_npy = outfile_npy
  )

    return job

def fmlrc(infile_npy, long_fasta_reads, outfile_fasta):
        
    job = Job(
        [infile_npy, long_fasta_reads],
        [outfile_fasta],
        [
            ['fmlrc', 'module_fmlrc']
        ]
    )
    job.command="""
fmlrc \\
  -p {num_threads} \\
  {infile_npy} {long_fasta_reads} \\
  {outfile_fasta}""".format(
    infile_npy = infile_npy,
    long_fasta_reads = long_fasta_reads,
    outfile_fasta = outfile_fasta,
    num_threads = config.param("fmlrc", "num_threads", 1, "posint")
  )

    return job

def compile_assembly_stats(infile, outfile):
        
    job = Job(
        [infile],
        [outfile],
        [
            ['meco_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
    job.command="""
compileAssemblyResultsSingle.pl \\
  --infile {infile} > {outfile}""".format(
    infile = infile,
    outfile = outfile
     )

    return job

def get_insert_size(indir, outfile, flagstats): #flagstats just here for dependency purpose.
    
    job = Job(
        [flagstats],
        [outfile],
        [
            
            ['python3', 'module_python3'],
            ['meco_tools', 'module_tools'],
            ['samtools', 'module_samtools']
        ]
    )
    
    job.command="""
getInsertSize.pl --indir {indir} > {outfile}""".format(
    indir = indir,
    outfile = outfile
    )

    return job

def sspace(infile, libraries, outdir):
    # module load mugqic_dev/SSPACE/3.0 && SSPACE_Standard_v3.0.pl -l libraries.txt 
    # -s ./assembly/Ray_k_31/Scaffolds.fasta -x 0 -m 32 -o 20 -k 5 -a 0.70 -n 15 -p 0 
    # -v 0 -z 0 -g 0 -T 1 -S 0 -b standard_out_k31
    librairies = outdir + "/libraries.txt"

    job = Job(
        [infile, libraries],
        [outdir],
        [
            
            ['meco_tools', 'module_tools'],
            ['sspace', 'module_sspace']
        ]
    )
    
    job.command="""
cat {infile} > {outdir}/libraries.txt \\
SSPACE_Standard_v3.0.pl \\
  -l {libraries} \\
  -s {infile} \\
  -x 0 -m 32 -o 20 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 1 -S 0 \\
  -b {outdir}""".format(
    librairies = librairies,
    infile = infile,
    outdir = outdir
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
prodigal -i {infile} -f gff -p meta \\
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

def make_index(fasta, bwt):

    job = Job(
        [fasta],
        [bwt],
        [
            
            ['samtools', 'module_samtools'],
            ['meco_tools', 'module_tools'],
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
            ['meco_tools', 'module_tools'],
            ['bwa', 'module_bwa']
        ]
    )

    job.command="""
fastaToBed.pl --fasta {fasta} > {bed}""".format(
    fasta = fasta,
    bed = bed
    )

    return job

def gff_to_bed(gff, bed):
        
    job = Job(
        [gff],
        [bed],
        [
            ['samtools', 'module_samtools'],
            ['meco_tools', 'module_tools'],
            ['bwa', 'module_bwa']
        ]
    )

    job.command="""
gffToBed.pl --infile {gff} > {bed}""".format(
    gff = gff,
    bed = bed
    )

    return job

def flagstats(infile, outfile):
        
    job = Job(
        [infile],
        [outfile],
        [
            ['samtools', 'module_samtools'],
            ['meco_tools', 'module_tools']
        ]
    )

    job.command="""
samtools flagstat {infile} > {outfile}""".format(
    infile = infile,
    outfile = outfile
    )

    return job

def merge_flagstats(infiles_fs, infiles_qc, infiles_bbduk, outfile, infiles_bbduk_sub=False):
    
    infiles = ""
    if(infiles_bbduk_sub is not False):
        infiles = infiles_fs + infiles_qc + infiles_bbduk + infiles_bbduk_sub
    else:
        infiles = infiles_fs + infiles_qc + infiles_bbduk


    job = Job(
        infiles,
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )

    job.command="""
mergeFlagStatsAndQC.pl \\
   --infilesFlagstat {infile_flagstat} \\
   --infilesQC {infile_qc} \\
   --infilesBBduk {infile_bbduk} \\""".format(
    infile_flagstat = ",".join(infiles_fs),
    infile_qc = ",".join(infiles_qc),
    infile_bbduk = ",".join(infiles_bbduk)
    )

    if(infiles_bbduk_sub is False):
        job.command += """
   > {outfile}""".format(
    outfile = outfile)

    else:
        job.command += """
    --infilesBBdukSub {infile_bbduk_sub} \\
    > {outfile}""".format(
    infile_flagstat = ",".join(infiles_fs),
    infile_qc = ",".join(infiles_qc),
    infile_bbduk = ",".join(infiles_bbduk),
    infile_bbduk_sub = ",".join(infiles_bbduk_sub),
    outfile = outfile)

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
  | samtools view -Sbh -F 0x100 -f 0x2 - > {outfile}.tmp && \\
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

def minimap2_samtools_se(reference, infile1, outfile):
    job = Job(
        [infile1, reference],
        [outfile],
        [
            ['samtools', 'module_samtools'],
            ['minimap2', 'module_minimap2']
        ]
    )
    job.command="""
minimap2 -ax ava-ont \\
  -t {num_threads} \\
  {reference} \\
  {infile1} \\
  | samtools view -Sbh -F 0x100 - > {outfile}.tmp && \\
  samtools sort -@ {num_threads} -m {mem_per_thread} {outfile}.tmp -o {outfile}.tmp.sorted && \\
  mv {outfile}.tmp.sorted.bam {outfile} && \\
  rm {outfile}.tmp && \\
  samtools index {outfile}""".format(
    num_threads = config.param('minimap2', 'num_threads', type='int', required=True), 
    mem_per_thread = config.param('samtools', 'mem_per_thread', required=True),
    infile1 = infile1,
    reference = reference,
    outfile = outfile
    )

    return job

def bwa_mem_samtools_se(reference, bwt, infile1, outfile):
    job = Job(
        [infile1, reference, bwt],
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
  | samtools view -Sbh -F 0x100 - > {outfile}.tmp && \\
  samtools sort -@ {num_threads} -m {mem_per_thread} {outfile}.tmp -o {outfile}.tmp.sorted && \\
  mv {outfile}.tmp.sorted.bam {outfile} && \\
  rm {outfile}.tmp && \\
  samtools index {outfile}""".format(
    num_threads = config.param('bwa', 'num_threads', type='int', required=True), 
    mem_per_thread = config.param('samtools', 'mem_per_thread', required=True),
    infile1 = infile1,
    reference = reference,
    outfile = outfile
    )

    return job

def coverage_bed_v2_24(infile, bed, outfile, flag):
    job = Job(
        [infile, bed],
        [outfile],
        [
            ['bedtools', 'module_bedtools'],
            ['samtools', 'module_samtools']
        ]
    )

#    job.command="""
#samtools view -b -f {flag} {infile} | \\
#  bedtools coverage -b stdin \\
#  -a {bed} \\
#  -counts -sorted \\
#  > {outfile}""".format(
#    infile = infile,
#    bed = bed,
#    outfile = outfile,
#    flag = flag
#    )
    
    job.command="""
bedtools multicov \
 -bams {infile} \
 -bed {bed} > {outfile}""".format(
    infile = infile,
    bed = bed,
    outfile = outfile
    )

    return job

# V.2.23 and lower. Takes far less RAM to complete.
def coverage_bed(infile, bed, outfile, flag, prefix=""):
    
    infile_base = os.path.splitext(os.path.basename(infile))[0]
    outdir = os.path.dirname(outfile)
    
    job = Job(
        [infile, bed],
        [outfile],
        [
            ['bedtools', 'module_bedtools'],
            ['samtools', 'module_samtools']
        ]
    )

    job.command="""
samtools view -b -f {flag} {infile} > {outdir}/{infile_base}{prefix}.tmp && \\
  coverageBed -abam {outdir}/{infile_base}{prefix}.tmp \\
  -b {bed} \\
  -counts \\
  > {outfile} && \\
rm -f {outdir}/{infile_base}{prefix}.tmp""".format(
    infile = infile,
    bed = bed,
    outdir = outdir,
    infile_base = infile_base,
    outfile = outfile,
    flag = flag,
    prefix = prefix
    )

    return job

def merge_counts(infiles, outfile_raw, type_feature):

    job = Job(
        infiles,
        [outfile_raw],
        [
            ['tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
  
    job.command="""
mergeAbundance.pl \\
  --type {type_feature} \\
  --infiles {infile} \\
  > {outfile_raw}""".format(
            infile = ",".join(infiles),
            outfile_raw = outfile_raw,
            type_feature = type_feature
        )

    return job

def convert_tsv_to_hdf5(infile, outfile):

    job = Job(
        [infile],
        [outfile],
        [
            ['tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
  
    job.command="""
rm {outfile} && \
convertTsvToHdf5.R \\
  -i {infile} \\
  -o {outfile}""".format(
            infile = infile,
            outfile = outfile
        )

    return job

def normalize_counts(infile_h5, outfile_cpm):

    job = Job(
        [infile_h5],
        [outfile_cpm],
        [
            ['tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
  
    job.command="""
generateCPMsH5.R \\
  -i {infile_h5} \\
  -o {outfile_cpm}""".format(
            infile_h5 = infile_h5,
            outfile_cpm = outfile_cpm
        )

    return job

def blastn(infile, outfile, outdir, prefix):
    job = Job(
        [infile],
        [outfile],
        [
            ['blast', 'module_blast'],
            ['meco_tools', 'module_tools']
        ]
    )

    job.command="""
mkdir -p {blast_dir} && \\
if [[ -s {infile} ]] ; then
    blastn \\
     -db {db} \\
     -query {infile} \\
     -out {outfile} \\
     -word_size {word_size} \\
     -outfmt \\"{outfmt}\\" \\
     -max_target_seqs 1 \\
     -num_threads {num_threads} \\
     {other_params}
else
    touch {outfile}
fi ;""".format(
    db = config.param(prefix, 'db', 1, 'string'),
    infile = infile,
    outfile = outfile,
    word_size = str(config.param(prefix, 'word_size', 1, 'int')),
    outfmt = config.param(prefix, 'outfmt', required=True),
    num_threads = config.param(prefix, 'num_threads', 1, 'int'),
    other_params = config.param(prefix, 'other_params', 1, 'string'),
    blast_dir = outdir
    )

    return job


# Because array job : touch done instead of relying on the standard done mechanism.
def blastn_array_job(indir, inprefix, outdir, outprefix, prefix, infiles, outfiles, dones, job_scheduler):

    job_array_suffix = ""
    if job_scheduler == "sge":
        job_array_suffix = "SGE_TASK_ID"
    else:
        job_array_suffix = "SLURM_ARRAY_TASK_ID"

    infile = os.path.join(indir, inprefix + "\$" + job_array_suffix + "2")
    outfile = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.tsv")
    done = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.done")

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
            ['meco_tools', 'module_tools']
        ]
    )
    job.command = ""
    
    if job_scheduler == "sge":
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 1  ))
        """.format(job_array_suffix = job_array_suffix)
    else:
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 0  ))
        """.format(job_array_suffix = job_array_suffix)
    
    job.command+="""
{job_array_suffix}2=\$(printf "%07d" \${job_array_suffix})
mkdir -p {blast_dir} && \\
if [[ -s {infile} ]] ; then
    blastn \\
     -db {db} \\
     -query {infile} \\
     -out {outfile} \\
     -word_size {word_size} \\
     -outfmt \\"{outfmt}\\" \\
     -max_target_seqs 1 \\
     -num_threads {num_threads} \\
     {other_params} \\
    && touch {done} 
else
    touch {outfile}
fi ;""".format(
    job_array_suffix = job_array_suffix,
    infile = infile,
    outfile = outfile,
    db = config.param(prefix, 'db', 1, 'string'),
    infiles = infiles,
    outfiles = outfiles,
    word_size = str(config.param(prefix, 'word_size', 1, 'int')),
    outfmt = config.param(prefix, 'outfmt', required=True),
    num_threads = config.param(prefix, 'num_threads', 1, 'int'),
    other_params = config.param(prefix, 'other_params', 1, 'string'),
    blast_dir = outdir,
    SLURM_ARRAY_TASK_ID = job_array_suffix,
    done = done
    )

    return job

# Because array job : touch done instead of relying on the standard done mechanism.
def blastp_array_job(indir, inprefix, outdir, outprefix, prefix, infiles, outfiles, dones, job_scheduler):

    # For SLURM:
    #job_array_suffix = "SLURM_ARRAY_TASK_ID"
    # For SGE
    job_array_suffix = ""
    if job_scheduler == "sge":
        job_array_suffix = "SGE_TASK_ID"
    else:
        job_array_suffix = "SLURM_ARRAY_TASK_ID"

    infile = os.path.join(indir, inprefix + "\$" + job_array_suffix + "2")
    outfile = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.tsv")
    done = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.done")

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
            ['meco_tools', 'module_tools']
        ]
    )
    job.command = ""
    
    if job_scheduler == "sge":
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 1  ))
        """.format(job_array_suffix = job_array_suffix)
    else:
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 0  ))
        """.format(job_array_suffix = job_array_suffix)
    
    job.command+="""
{job_array_suffix}2=\$(printf "%07d" \${job_array_suffix})
mkdir -p {blast_dir} && \\
if [[ -s {infile} ]] ; then
    blastp \\
     -db {db} \\
     -query {infile} \\
     -out {outfile} \\
     -word_size {word_size} \\
     -outfmt \\"{outfmt}\\" \\
     -max_target_seqs 1 \\
     -num_threads {num_threads} \\
    && touch {done} 
else
    touch {outfile}
fi ;""".format(
    job_array_suffix = job_array_suffix,
    infile = infile,
    outfile = outfile,
    db = config.param(prefix, 'db', 1, 'string'),
    infiles = infiles,
    outfiles = outfiles,
    word_size = str(config.param(prefix, 'word_size', 1, 'int')),
    outfmt = config.param(prefix, 'outfmt', required=True),
    num_threads = config.param(prefix, 'num_threads', 1, 'int'),
    blast_dir = outdir,
    done = done
    )

    return job

def diamond_blastp_nr_array_job(indir, inprefix, outdir, outprefix, prefix, infiles, outfiles, db, dones, job_scheduler):
    
    job_array_suffix = ""
    if job_scheduler == "sge":
        job_array_suffix = "SGE_TASK_ID"
    else:
        job_array_suffix = "SLURM_ARRAY_TASK_ID"

    infile = os.path.join(indir, inprefix + "\$" + job_array_suffix + "2")
    outfile = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.tsv")
    done = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.done")

    all_outfiles = []
    for i in range(0, len(dones)):
        all_outfiles.append(dones[i])
    for i in range(0, len(outfiles)):
        all_outfiles.append(outfiles[i])
    

    job = Job(
        infiles,
        all_outfiles,
        [
            ['diamond', 'module_diamond'],
            ['meco_tools', 'module_tools']
        ]
    )
    job.command = ""
    
    if job_scheduler == "sge":
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 1  ))
        """.format(job_array_suffix = job_array_suffix)
    else:
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 0  ))
        """.format(job_array_suffix = job_array_suffix)
    
    job.command+="""
{job_array_suffix}2=\$(printf "%07d" \${job_array_suffix}2)
echo 'SLURM/SGE_TASK_ID2' \${job_array_suffix}2
mkdir -p {blast_dir} && \\
if [[ -s {infile} ]] ; then
    diamond blastp \\
     -d {db} \\
     -q {infile} \\
     -o {outfile} \\
     -k 10 \\
     -e {evalue} \\
     -p {num_threads} \\
    && touch {done} 
else
    touch {outfile}
fi ;""".format(
    infile = infile,
    outfile = outfile,
    db = db,
    evalue = config.param('blastp_nr', 'evalue', 1, 'string'),
    num_threads = config.param('blastp_nr', 'num_threads', 1, 'int'),
    blast_dir = outdir,
    done = done,
    job_array_suffix = job_array_suffix
    )

    return job

def diamond_blastp_kegg_array_job(indir, inprefix, outdir, outprefix, prefix, infiles, outfiles, db, dones, job_scheduler=""):

    job_array_suffix = ""
    if job_scheduler == "sge":
        job_array_suffix = "SGE_TASK_ID"
    else:
        job_array_suffix = "SLURM_ARRAY_TASK_ID"

    infile = os.path.join(indir, inprefix + "\$" + job_array_suffix + "2")
    outfile = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.tsv")
    done = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.done")

    all_outfiles = []
    for i in range(0, len(dones)):
        all_outfiles.append(dones[i])
    for i in range(0, len(outfiles)):
        all_outfiles.append(outfiles[i])

    job = Job(
        infiles,
        all_outfiles,
        [
            ['diamond', 'module_diamond'],
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command = ""
    #sys.stderr.write("Job scheduler: " + job_scheduler + "\n")
    if job_scheduler == "sge":
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 1  ))
        """.format(job_array_suffix = job_array_suffix)
    else:
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 0  ))
        """.format(job_array_suffix = job_array_suffix)
    
    job.command+="""
{job_array_suffix}2=\$(printf "%07d" \${job_array_suffix}2)
echo 'SLURM/SGE_TASK_ID2' \${job_array_suffix}2
mkdir -p {blast_dir} && \\
if [[ -s {infile} ]] ; then
    diamond blastp \\
     -d {db} \\
     -q {infile} \\
     -o {outfile} \\
     -k 1 \\
     -e {evalue} \\
     -p {num_threads} \\
    && touch {done} 
else
    touch {outfile}
fi ;""".format(
    job_array_suffix = job_array_suffix,
    infile = infile,
    outfile = outfile,
    db = db,
    evalue = config.param('blastp_kegg', 'evalue', 1, 'string'),
    num_threads = config.param('blastp_kegg', 'num_threads', 1, 'int'),
    blast_dir = outdir,
    done = done
    )

    return job

def ublastp(infile, outfile, outdir, db):
    job = Job(
        [infile],
        [outfile],
        [
            
            ['usearch', 'module_usearch'],
            ['meco_tools', 'module_tools']
        ]
    )

    job.command="""
mkdir -p {outdir} && \\
usearch -ublast \\
  {infile} \\
  -db {db} \\
  -evalue {evalue} \\
  -accel {accel} \\
  -blast6out {outfile} \\
  -maxhits 1 \\
  -threads {num_threads}""".format(
    db = db,
    infile = infile,
    outfile = outfile,
    accel = config.param('ublastp', 'accel', 1, 'float'),
    num_threads = config.param('ublastp', 'num_threads', 1, 'int'),
    outdir = outdir,
    evalue = config.param('ublastp', 'evalue', 1, 'string')
    )

    return job

def keep_ublast_best_hit(indir, outfile, number_of_chunks, prefix):

    infiles = []
    for i in range(1,number_of_chunks):
        infiles.append(os.path.join(indir, prefix + "_chunk_{:03d}.tsv".format(i)))

    job = Job(
        infiles, 
        [outfile],
        [
            
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
keepUblastpBestHit.pl \\
  --infiles {infiles} > {outfile}""".format(
    outfile = outfile,
    infiles = ",".join(infiles)
    )

    return job

def annotate_ublast(infile, outfile, annotations):

    job = Job(
        [infile],
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
fetchNRAnnotation.pl \\
  --blast {infile} \\
  --annotations {annotations} \\
  > {outfile}""".format(
    outfile = outfile,
    infile = infile,
    annotations = annotations
    )

    return job

# For KEGG
def blastp(infile, outfile, outdir, db):
    job = Job(
        [infile],
        [outfile],
        [
            ['blast', 'module_blast'],
            ['meco_tools', 'module_tools']
        ]
    )

    job.command="""
mkdir -p {blast_dir} && \\
blastp \\
  -db {db} \\
  -query {infile} \\
  -out {outfile} \\
  -outfmt {outfmt} \\
  -max_target_seqs 1 \\
  -evalue {evalue} \\
  -num_threads {num_threads}""".format(
    db = db,
    infile = infile,
    outfile = outfile,
    outfmt = config.param('blastp', 'outfmt', required=True),
    num_threads = config.param('blastp', 'num_threads', 1, 'int'),
    blast_dir = outdir,
    evalue = config.param('blastp', 'evalue', required=True)
    )

    return job

def parse_kegg(infile, outfile):

    job = Job(
        [infile], 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
parseKegg.pl \\
  --infile {infile} \\
  --ko {ko} \\
  --genes_desc {genes_desc} \\
  --genetoko {genetoko} \\
  > {outfile}""".format(
    outfile = outfile,
    infile = infile,
    ko = config.param('kegg', 'ko', required=True),
    genetoko = config.param('kegg', 'genetoko', required=True),
    genes_desc = config.param('kegg', 'genes_desc', required=True)
    )

    return job

def parse_cazy(infile, infile_reference, outfile):

    job = Job(
        [infile, infile_reference], 
        [outfile],
        [
            ['perl', 'module_perl'],
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
parseCazy.pl \\
  --infile {infile} \\
  --infile_reference {infile_reference} \\
  > {outfile}""".format(
    outfile = outfile,
    infile = infile,
    infile_reference = infile_reference
    )

    return job

def parse_kofam(infile, outfile, hmmsearch=False):

    job = Job(
        [infile], 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    if hmmsearch:
        job.command="""
parseKofam.pl --hmmsearch \\"""
    else:
        job.command="""
parseKofam.pl \\"""

    job.command+="""
  --infile {infile} \\
  --ref_database {ref_file} \\
  > {outfile}""".format(
    outfile = outfile,
    infile = infile,
    ref_file = config.param('parse_kofam', 'ref_database', required=True, type='filepath'),
    )

    return job

def keep_blast_best_hit(infile, outfile):

    job = Job(
        [infile], 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
keepBestBlastHit.pl \\
  --e {evalue} \\
  --length {length} \\
  --perc {perc} \\
  --infile {infile} > {outfile}""".format(
    outfile = outfile,
    infile = infile,
    evalue =  config.param('keep_best_hit', 'evalue', required=True),
    length =  config.param('keep_best_hit', 'length', required=True),
    perc =  config.param('keep_best_hit', 'align_perc', required=True)
    )

    return job

def extract_taxonomy(infile, outfile, accession_to_tax):

    job = Job(
        [infile],
        [outfile], 
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
getNcbiTax.pl \\
  --blast_table {infile} \\
  --accession_to_tax {accession_to_tax} \\
  --names {names} \\
  --cutoff {cutoff} \\
  --length {length} > {outfile}""".format(
    outfile = outfile,
    infile = infile,
    accession_to_tax = accession_to_tax,
    names = config.param('ncbi_tax', 'names', required=True),
    cutoff = config.param('ncbi_tax', 'cutoff', required=True),
    length = config.param('ncbi_tax', 'alignment_length', required=True)
    )

    return job

def generate_feature_table(infile_tax, infile_abundance, outfile):

    job = Job(
        [infile_tax, infile_abundance],
        [outfile], 
        [
            ['meco_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
generateOTUTable.R \\
  -t {infile_tax} \\
  -a {infile_abundance} \\
  -o {outfile}""".format(
    outfile = outfile,
    infile_tax = infile_tax,
    infile_abundance = infile_abundance
    )

    return job

def generate_feature_table_rrna(infile_tax, infile_abundance, outfile):

    job = Job(
        [infile_tax, infile_abundance],
        [outfile], 
        [
            ['meco_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
generateFeatureTableRrna.R \\
  -t {infile_tax} \\
  -a {infile_abundance} \\
  -o {outfile}""".format(
    outfile = outfile,
    infile_tax = infile_tax,
    infile_abundance = infile_abundance
    )

    return job

def generate_feature_table_reads_centric(infile_rdps, names, outfile):

    job = Job(
        infile_rdps,
        [outfile], 
        [
            ['meco_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
generateOTUTableReadsCentric.pl \\
  --infiles {infile_rdps} \\
  --names {names} \\
  > {outfile}""".format(
    outfile = outfile,
    infile_rdps = ",".join(infile_rdps),
    names = ",".join(names)
    )

    return job

# Because array job : touch done instead of relying on the standard done mechanism.
def rpsblast_array_job(indir, inprefix, outdir, outprefix, infiles, outfiles, dones, db, job_scheduler=""):

    job_array_suffix = ""
    if job_scheduler == "sge":
        job_array_suffix = "SGE_TASK_ID"
    else:
        job_array_suffix = "SLURM_ARRAY_TASK_ID"

    infile = os.path.join(indir, inprefix + "\$" + job_array_suffix + "2")
    outfile = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.tsv")
    done = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.done")

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
            ['meco_tools', 'module_tools']
        ]
    )

    job.command = ""
    
    if job_scheduler == "sge":
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 1  ))
        """.format(job_array_suffix = job_array_suffix)
    else:
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 0  ))
        """.format(job_array_suffix = job_array_suffix)
    
    job.command+="""
{job_array_suffix}2=\$(printf "%07d" \${job_array_suffix}2)
echo 'SLURM/SGE_TASK_ID2' \${job_array_suffix}2
mkdir -p {blast_dir} && \\
rpsblast \\
  -db {db} \\
  -query {infile} \\
  -out {outfile} \\
  -outfmt \\"{outfmt}\\" \\
  -max_target_seqs 1 \\
  -evalue {evalue} \\
  -num_threads {num_threads} \\
&& touch {done}""".format(
    db = db,
    infile = infile,
    outfile = outfile,
    outfmt = config.param('rpsblast', 'outfmt'),
    num_threads = config.param('rpsblast', 'num_threads', 1, 'int'),
    blast_dir = outdir,
    evalue = config.param('rpsblast', 'evalue'),
    job_array_suffix = job_array_suffix,
    done = done
    )

    return job

def rpsblast(infile, outfile, outdir, db):
    job = Job(
        [infile],
        [outfile],
        [
            ['blast', 'module_blast'],
            ['meco_tools', 'module_tools']
        ]
    )

    job.command="""
mkdir -p {blast_dir} && \\
rpsblast \\
  -db {db} \\
  -query {infile} \\
  -out {outfile} \\
  -outfmt \\"{outfmt}\\" \\
  -max_target_seqs 1 \\
  -evalue {evalue} \\
  -num_threads {num_threads}""".format(
    db = db,
    infile = infile,
    outfile = outfile,
    outfmt = config.param('rpsblast', 'outfmt'),
    num_threads = config.param('rpsblast', 'num_threads', 1, 'int'),
    blast_dir = outdir,
    evalue = config.param('rpsblast', 'evalue')
    )

    return job

def hmmsearch_array_job(indir, inprefix, outdir, outprefix, infiles, 
                      tblouts, domtblouts, pfamtblouts, 
                      dones, db, job_scheduler):
    
    job_array_suffix = ""
    if job_scheduler == "sge":
        job_array_suffix = "SGE_TASK_ID"
    else:
        job_array_suffix = "SLURM_ARRAY_TASK_ID"

    infile = os.path.join(indir, inprefix + "\$" + job_array_suffix + "2")
    tblout = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.tblout")
    domtblout = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.domtblout")
    #pfamtblout = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.pfamtblout")
    done = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.done")

    all_outfiles = []
    for i in range(0, len(dones)):
        all_outfiles.append(dones[i])
    for i in range(0, len(tblouts)):
        all_outfiles.append(tblouts[i])
   
    job = Job(
        infiles,
        all_outfiles,
        [
            ['meco_tools', 'module_tools'],
            ['hmmer', 'module_hmmer']
        ]
    )

    job.command = ""
    
    if job_scheduler == "sge":
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 1  ))
        """.format(job_array_suffix = job_array_suffix)
    else:
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 0  ))
        """.format(job_array_suffix = job_array_suffix)
    
    job.command+="""
{job_array_suffix}2=\$(printf "%07d" \${job_array_suffix}2)
echo 'SLURM/SGE_TASK_ID2' \${job_array_suffix}2
mkdir -p {hmm_dir} && \\
hmmsearch \\
  --tblout {tblout} \\
  --domtblout {domtblout} \\
  -E {evalue} \\
  --cpu {num_threads} \\
  {db} \\
  {infile} > /dev/null && \\
touch {done}""".format(
    db = db,
    infile = infile,
    tblout = tblout,
    domtblout =  domtblout,
    num_threads = config.param('hmmsearch', 'num_threads', 1, 'int'),
    evalue = config.param('hmmsearch', 'evalue'),
    hmm_dir = outdir,
    done = done,
    job_array_suffix = job_array_suffix
    )

    return job

def parse_hmms(infile, outfile):
    job = Job(
        [infile],
        [outfile],
        [
            ['perl', 'module_perl'],
            ['meco_tools', 'module_tools']
        ]
    )

    job.command="""
parseHmmsearch.pl \\
  --n 1 --e {evalue} --infile {infile} \\
  > {outfile} \\
  """.format(
    infile = infile,
    outfile = outfile,
    evalue = config.param('hmmsearch', 'evalue')
    )

    return job

def hmmscan(infile, prefix, outdir, db):
    
    tblout = os.path.join(prefix + ".tblout")
    domtblout = os.path.join(prefix + ".domtblout")
    pfamtblout = os.path.join(prefix + ".pfamtblout")

    job = Job(
        [infile],
        [tblout, domtblout, pfamtblout],
        [
            ['blast', 'module_blast'],
            ['meco_tools', 'module_tools'],
            ['hmmer', 'module_hmmer']
        ]
    )

    job.command="""
mkdir -p {hmm_dir} && \\
hmmscan \\
  --tblout {tblout} \\
  --domtblout {domtblout} \\
  --pfamtblout {pfamtblout} \\
  -E {evalue} \\
  --cpu {num_threads} \\
  {db} \\
  {infile} > /dev/null""".format(
    db = db,
    infile = infile,
    tblout = tblout,
    domtblout =  domtblout,
    pfamtblout = pfamtblout,
    num_threads = config.param('hmmscan', 'num_threads', 1, 'int'),
    evalue = config.param('hmmscan', 'evalue'),
    hmm_dir = outdir
    )

    return job

def split_rrna(infile_tab, infile_fasta,
               arc_5S, bac_5S, euk_5S, arc_16S, bac_16S, euk_18S, arc_23S, bac_23S, euk_28S, arc_bac_16S):

    job = Job(
        [infile_tab, infile_fasta], 
        [arc_5S, bac_5S, euk_5S, arc_16S, bac_16S, euk_18S, arc_23S, bac_23S, euk_28S, arc_bac_16S],
        [
            
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
splitRrna.pl \\
  --infile_tab {infile_tab} --infile_fasta {infile_fasta} \\
  --arc_5S {arc_5S} --bac_5S {bac_5S} --euk_5S {euk_5S} --arc_16S {arc_16S} --bac_16S {bac_16S} --euk_18S {euk_18S} --arc_23S {arc_23S} --bac_23S {bac_23S} --euk_28S {euk_28S} && \\
  cat {arc_16S} {bac_16S} > {arc_bac_16S}""".format(
    infile_tab = infile_tab,
    infile_fasta = infile_fasta,
    arc_5S = arc_5S,
    bac_5S = bac_5S,
    euk_5S = euk_5S,
    arc_16S = arc_16S,
    bac_16S = bac_16S,
    euk_18S = euk_18S,
    arc_23S = arc_23S,
    bac_23S = bac_23S,
    euk_28S = euk_28S,
    arc_bac_16S = arc_bac_16S
    )

    return job

def rrna_to_bed(infile, outfile):
        
    job = Job(
        [infile],
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    job.command="""
if [[ -s {infile} ]] ; then
    rnammerToBed.pl --infile {infile} > {outfile}
else
    touch {outfile}
fi""".format(
    infile = infile,
    outfile = outfile
    )

    return job

def barrnap_fna_to_bed(infile, outfile):
        
    job = Job(
        [infile],
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    job.command="""
if [[ -s {infile} ]] ; then
    barrnapToBed.pl --infile {infile} > {outfile}
else
    touch {outfile}
fi""".format(
    infile = infile,
    outfile = outfile
    )

    return job

def rdp_to_taxonomy(infile, outfile):

    job = Job(
        [infile],
        [outfile], 
        [
            ['meco_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""

if [[ -s {infile} ]] ; then
    rdpToTaxonomy.pl \\
     --infile {infile} \\
     > {outfile}
else
    touch {outfile}
fi""".format(
    outfile = outfile,
    infile = infile
    )

    return job

def estimate_number_of_chunks(infile, outfile, chunk_file_size):

    job = Job(
        [infile],
        [outfile], 
        [
            ['meco_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
estimateChunkFileSize.pl \\
  --infile {infile} \\
  --targeted_chunk_size {chunk_file_size} > \\
  {outfile}""".format(
    outfile = outfile,
    infile = infile,
    chunk_file_size = chunk_file_size
    )

    return job

def exonerate(infile, outdir, number_of_chunks_file, prefix):

    outfile = os.path.join(outdir, "exonerate.done")

    filename, file_extension = os.path.splitext(infile)
    #sys.stderr.write("file_extension: " + file_extension + "\n")
    to_delete = ""
    if file_extension == ".faa":
        to_delete = outdir + "/*.faa"
        prefix2 = prefix + "_chunk_"
    elif file_extension == ".fna":
        to_delete = outdir + "/*.fna"
        prefix2 = prefix + "_chunk_"
    elif file_extension == ".fasta":
        to_delete = outdir + "/*.fasta"
        prefix2 = prefix + "_chunk_"

    job = Job(
        [infile, number_of_chunks_file], 
        [outfile],
        [
            ['exonerate', 'module_exonerate'],
            ['meco_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
   
    job.command="""
rm -rf {to_delete} && mkdir -p {to_delete} && \\
fastasplit -f {infile} \\
  -o {outdir} \\
  -c \\`cat {number_of_chunks_file}\\` &&
doublecheckExonerate.pl --indir {outdir} --prefix {prefix2} --no_chunks \\`cat {number_of_chunks_file}\\` && \\
 touch {outfile} && \\
 rm -rf {to_delete}""".format(
    outdir = outdir,
    number_of_chunks_file = number_of_chunks_file,
    infile = infile,
    to_delete = to_delete,
    prefix2 = prefix2,
    outfile = outfile
    )

    return job

def exonerate_direct(infile, outdir, num_chunks, prefix):

    outfile = os.path.join(outdir, "exonerate.done")
    outfiles = []
    for i in range(num_chunks):
        outfiles.append(os.path.join(outdir, prefix + "_chunk_{:07d}".format(i)))

    filename, file_extension = os.path.splitext(infile)
    #sys.stderr.write("file_extension: " + file_extension + "\n")
    #to_delete = ""
    if file_extension == ".faa":
        to_delete = outdir + "/*.faa"
        prefix2 = prefix + "_chunk_"
    elif file_extension == ".fna":
        to_delete = outdir + "/*.fna"
        prefix2 = prefix + "_chunk_"
    elif file_extension == ".fasta":
        to_delete = outdir + "/*.fasta"
        prefix2 = prefix + "_chunk_"

    job = Job(
        [infile], 
        [outfile] + outfiles,
        [
            ['exonerate', 'module_exonerate'],
            ['meco_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
   
    job.command="""
rm -rf {to_delete} && mkdir -p {to_delete} && \\
fastasplit -f {infile} \\
  -o {outdir} \\
  -c {num_chunks} &&
doublecheckExonerate.pl --indir {outdir} --prefix {prefix2} --no_chunks {num_chunks} && \\
 touch {outfile} && \\
 rm -rf {to_delete}""".format(
    outdir = outdir,
    num_chunks = num_chunks,
    infile = infile,
    to_delete = to_delete,
    prefix2 = prefix2,
    outfile = outfile
    )

    return job

def merge_chunks(indir, outfile, number_of_chunks, prefix):

    infiles = []
    for i in range(number_of_chunks):
        infiles.append(os.path.join(indir, prefix + "_chunk_{:07d}.tsv".format(i)))

    job = Job(
        infiles, 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
cat {infiles} \\
  > {outfile}""".format(
    outfile = outfile,
    infiles = " ".join(infiles)
    )

    return job

def edger(abundance, design_file, outdir):

    dummy_outfile = os.path.join(outdir, "DDA.done")

    job = Job(
        [abundance, design_file], 
        [dummy_outfile],
        [
            
            ['meco_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
edger.R \\
  -i {abundance} \\
  -o {outdir} \\
  -d {design} \\
  -p {pvalue} \\
  -f {fdr} \\
  -l {logfc} \\
  && touch {dummy_outfile}""".format(
    abundance = abundance,
    outdir = outdir,
    design = design_file,
    dummy_outfile = dummy_outfile,
    logfc = config.param('DDA', 'logfc', 1, 'float'),
    fdr = config.param('DDA', 'fdr', 1, 'float'),
    pvalue = config.param('DDA', 'pvalue', 1, 'float')
    )

    return job

def edger_glm(abundance, mapping_file, outdir):

    basename = os.path.splitext(os.path.basename(mapping_file))[0]
    dummy_outfile = os.path.join(outdir, basename + "_DDA_GLM.done")

    job = Job(
        [abundance, mapping_file], 
        [dummy_outfile],
        [
            ['meco_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
edgerFeaturesGLM_SMG.R \\
  -i {abundance} \\
  -o {outdir} \\
  -m {mapping_file} \\
  -p {pvalue} \\
  -f {fdr} \\
  -l {logfc} \\
  -b {blocks} \\
  -t {treatments} \\
  && touch {dummy_outfile}""".format(
    abundance = abundance,
    outdir = outdir,
    dummy_outfile = dummy_outfile,
    logfc = config.param('DDA', 'logfc', 1, 'float'),
    fdr = config.param('DDA', 'fdr', 1, 'float'),
    pvalue = config.param('DDA', 'pvalue', 1, 'float'),
    blocks = config.param('DDA', 'blocks', 1, 'string'),
    treatments = config.param('DDA', 'treatments', 1, 'string'),
    mapping_file = mapping_file
    )

    return job

def generate_gff(infile_gff, infile_fasta, kegg, pfam, cog, kog, 
                 taxonomy, ublast_nr, cazy_parsed, 
                 outfile_gff, outfile_fasta, outfile_annotations):

    job = Job(
        [infile_gff, infile_fasta, kegg, pfam, cog, kog, taxonomy, ublast_nr, cazy_parsed], 
        [outfile_gff, outfile_fasta, outfile_annotations],
        [
            ['meco_tools', 'module_tools'],
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
  --cazy {cazy_parsed} \\
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
    cazy_parsed = cazy_parsed,
    outfile_gff = outfile_gff,
    outfile_fasta = outfile_fasta,
    outfile_annotations = outfile_annotations
    )

    return job

def generate_gc_table(infile, contigs, abundance, gc_table):

    job = Job(
        [infile, contigs, abundance], 
        [gc_table],
        [
            
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
generateGCtable.pl \\
  --infile {infile} \\
  --contigs {contigs} \\
  --abundance {abundance} \\
  > {gc_table}""".format(
    infile = infile,
    contigs = contigs,
    abundance = abundance,
    gc_table = gc_table
  )

    return job

def merge_chunks_hmms(indir, outdir, number_of_chunks, prefix, out_prefix):

    tblout = []
    domtblout = []
    pfamtblout = []
    for i in range(number_of_chunks):
        tblout.append(os.path.join(indir, prefix + "_chunk_{:07d}.tblout".format(i)))
        domtblout.append(os.path.join(indir, prefix + "_chunk_{:07d}.domtblout".format(i)))
        #pfamtblout.append(os.path.join(indir, prefix + "_chunk_{:07d}.pfamtblout".format(i)))

    tblout_out = os.path.join(outdir, out_prefix + "_tblout.tsv") 
    domtblout_out = os.path.join(outdir, out_prefix + "_domtblout.tsv") 
    #pfamtblout_out = os.path.join(outdir, out_prefix + "_pfamtblout.tsv") 

    job = Job(
        tblout, 
        [tblout_out, domtblout_out],
        [
            ['exonerate', 'module_exonerate'],
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
cat {tblout} \\
  > {tblout_out}.tmp && \\
cat {domtblout} \\
  > {domtblout_out}.tmp && \\
awk 'BEGIN{{OFS=FS=\\" \\"}} NR<=3{{print}}; NR>3{{tmp=\$1; \$1=\$3; \$3=tmp; tmp=\$2; \$2=\$4; \$4=tmp; print}}' {tblout_out}.tmp > {tblout_out} && rm {tblout_out}.tmp && \\
awk 'BEGIN{{OFS=FS=\\" \\"}} NR<=3{{print}}; NR>3{{tmp=\$1; \$1=\$4; \$4=tmp; tmp=\$2; \$2=\$5; \$5=tmp; print}}' {domtblout_out}.tmp > {domtblout_out} && rm {domtblout_out}.tmp""".format(
    tblout = " ".join(tblout),
    domtblout = " ".join(domtblout),
    tblout_out = tblout_out,
    domtblout_out = domtblout_out
    )

    return job

def kegg_overrep(infile, blastp_table, type, outfile):

    script = ""
    dummy_outfile = ""
    if type == "pathways":
        script = "getKeggPathways.py"
    if type == "modules":
        script = "getKeggModules.py"
    if type == "KO":
        script = "getKeggKO.py"

    job = Job(
        [blastp_table],
        [outfile],
        [
            ['tools', 'module_tools'],
            #['lapack', 'module_lapack'],
            ['python3', 'module_python3']
        ]
    )
    job.command="""
{script} \\
  --infile-blastp {blastp_table} \\
  --infile-gene-abundance {infile} \\
  > {outfile}""".format(
    blastp_table = blastp_table,
    infile = infile,
    script = script,
    outfile = outfile
    )

    return job

def cog_overrep(infile, blastp_table, outfile):

    job = Job(
        [blastp_table],
        [outfile],
        [
            ['tools', 'module_tools'],
            ['python3', 'module_python3']
        ]
    )
    job.command="""
getCOG.py \\
  --infile-blastp {blastp_table} \\
  --infile-gene-abundance {infile} \\
  > {outfile}""".format(
    blastp_table = blastp_table,
    infile = infile,
    outfile = outfile
    )

    return job

def pfam_overrep(infile, pfam_table, outfile):

    job = Job(
        [pfam_table],
        [outfile],
        [
            ['tools', 'module_tools'],
            ['python3', 'module_python3']
        ]
    )
    job.command="""
getPfam.py \\
  --infile-pfam {pfam_table} \\
  --infile-gene-abundance {infile} \\
  > {outfile}""".format(
    pfam_table = pfam_table,
    infile = infile,
    outfile = outfile
    )

    return job

def kegg_matrix(blastp_table, type, outfile):

    script = ""
    dummy_outfile = ""
    if type == "pathways":
        script = "generateKeggPathwaysMatrix.py"
    if type == "modules":
        script = "generateKeggModulesMatrix.py"
    if type == "K":
        script = "generateKeggKsMatrix.py"

    job = Job(
        [blastp_table],
        [outfile],
        [
            ['tools', 'module_tools'],
            ['lapack', 'module_lapack'],
            ['qiime-dependencies', 'module_qiime-dependencies']
        ]
    )
    job.command="""
{script} \\
  --infile-blastp {blastp_table} \\
  > {outfile}""".format(
    blastp_table = blastp_table,
    script = script,
    outfile = outfile
    )

    return job

def get_abundance_for_metabat(infiles, outfile, paired_contigs):

    job = Job(
        infiles, 
        [outfile],
        [
            ['metabat', 'module_metabat2'],
            ['tools','module_tools']
        ]
    )
    job.command="""
export OMP_NUM_THREADS={num_threads} && \\
jgi_summarize_bam_contig_depths \\
  --outputDepth {outfile} \\
  --pairedContigs {paired_contigs} \\
  --minContigLength {min_contig_length} \\
  --minContigDepth {min_contig_depth} \\
  --percentIdentity {perc_id} \\
  {infiles}""".format(
    infiles = " ".join(infiles),
    outfile = outfile,
    min_contig_length = config.param('metabat_abundance', 'min_contig_length', 1, 'posint'),
    min_contig_depth = config.param('metabat_abundance', 'min_contig_depth', 1, 'posint'),
    paired_contigs = paired_contigs,
    perc_id =  config.param('metabat_abundance', 'perc_id', 1, 'posint'),
    num_threads =  config.param('metabat_abundance', 'num_threads', 1, 'posint')
    )

    return job

def metabat2(contigs, abundance, outdir, dummy_outfile):

    saved_tnf = os.path.join(outdir, "saved.tnf")
    saved_dist = os.path.join(outdir, "saved.dist")

    job = Job(
        [contigs, abundance], 
        [dummy_outfile],
        [
            ['metabat','module_metabat2']
        ]
    )
    
    job.command="""
rm -rf {outdir} && \\
mkdir -p {outdir} && \\
metabat2 \\
  -i {contigs} \\
  -a {abundance} \\
  -o {outdir} \\
  --maxP {max_p} \\
  -m {min_contig} \\
  -t {num_threads} \\
  --saveTNF {saved_tnf} --saveDistance {saved_dist} \\
&& touch {dummy_outfile}""".format(
    contigs = contigs,
    abundance = abundance,
    outdir = outdir,
    dummy_outfile = dummy_outfile,
    saved_tnf = saved_tnf,
    saved_dist = saved_dist,
    num_threads = config.param('metabat2', 'num_threads', 1, 'posint'),
    max_p = config.param('metabat2', 'max_p', 1, 'posint'),
    min_contig = config.param('metabat2', 'min_contig', 1, 'posint')
    )

    return job

def maxbin2(contigs, abundance, outdir, out_prefix, dummy_outfile):

    job = Job(
        [contigs, abundance], 
        [dummy_outfile],
        [
            ['maxbin2','module_maxbin2'],
            ['fraggenescan','module_fraggenescan'],
            ['hmmer3','module_hmmer'],
            ['meco_tools','module_tools']
        ]
    )
   #TODO rename .fasta to .fa after maxin2 has run 
    job.command="""
rm -rf {outdir} && \\
mkdir -p {outdir} && \\
splitAbundanceForMaxbin.R \\
    -i {abundance} \\
    -o {outdir} && \\
find ./{outdir}/ -name \\\"*.txt\\\" | grep .txt > {outdir}/abundance_files.tmp && mv {outdir}/abundance_files.tmp {outdir}/abundance_files.txt && \\
run_MaxBin.pl \\
    -thread {num_threads} \\
    -contig {contigs} \\
    -out {out_prefix} \\
    -abund_list {outdir}/abundance_files.txt \\
    -min_contig_length {min_contig} \\
    -max_iteration {max_iteration} \\
    -prob_threshold {prob_threshold} \\
    -plotmarker && \\
for f in {outdir}/*.fasta; do mv -- \\\"\$f\\\" \\\"\${{f%.fasta}}.fa\\\"; done && \\
touch {dummy_outfile}""".format(
    contigs = contigs,
    abundance = abundance,
    outdir = outdir,
    out_prefix = out_prefix,
    dummy_outfile = dummy_outfile,
    num_threads = config.param('maxbin2', 'num_threads', 1, 'posint'),
    max_iteration = config.param('maxbin2', 'max_iteration', 1, 'posint'),
    min_contig = config.param('maxbin2', 'min_contig', 1, 'posint'),
    prob_threshold = config.param('maxbin2', 'prob_threshold', 1, 'float')
    )

    return job

def checkm(dummy_infile, indir, outdir, outfile):

    job = Job(
        [dummy_infile], 
        [outfile],
        [
            
            ['hmmer','module_hmmer'],
            ['prodigal','module_prodigal'],
            ['pplacer','module_pplacer'],
            ['checkm','module_checkm']
        ]
    )
    
    job.command="""
bash -c 'set +u && source \$CHECKM_HOME/checkm_venv/bin/activate && set -u &&
rm -rf {outdir} && \\
mkdir -p {outdir} && \\
\$CHECKM_HOME/checkm_venv/bin/checkm lineage_wf \\
  -f {outfile} \\
  -x fa \\
  -t {num_threads} \\
  --tab_table""".format(
        outfile = outfile,
        num_threads = config.param('checkm', 'num_threads', 1, 'posint'),
        outdir = outdir,
        indir = indir
    )
  
    if config.param('DEFAULT', 'tmpdir', 0, 'string') != "":
        sys.stderr.write("tmpdir:" + config.param('DEFAULT', 'tmpdir', 0, 'string') + "\n")
        job.command+=""" --tmpdir {tmpdir} {indir} {outdir}'""".format(
            tmpdir = config.param('DEFAULT', 'tmpdir', 0, 'string'),
            outdir = outdir,
            indir = indir
        )
    else:
        job.command+=""" {indir} {outdir}'""".format(
            outdir = outdir,
            indir = indir
        )
    
    return job

def parse_bins(dummy_infile, indir, outdir, outfile, outfile_raw, taxonomy):

    job = Job(
        [dummy_infile, taxonomy], 
        [outfile, outfile_raw],
        [
            ['tools', 'module_tools']
        ]
    )
    
    job.command="""
rm -rf {outdir} && \\
mkdir -p {outdir} && \\
validateMetabatBins.pl \\
  --indir {indir} \\
  --taxonomy {taxonomy} \\
  --split {split} \\
  --outdir {outdir} --outfile_link_raw {outfile_raw} > {outfile}""".format(
    outdir = outdir,
    indir = indir,
    outfile = outfile,
    taxonomy = taxonomy,
    outfile_raw = outfile_raw,
    split = config.param('parse_bins', 'split', 1, 'string')
    )

    return job

def summarize_bins(checkm, parsed_bins, contigs_abundance, contigs_indir, outfile, bins_list):

    job = Job(
        [checkm, parsed_bins, contigs_abundance], 
        [outfile, bins_list],
        [
            ['meco_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
    
    job.command="""
addReadsCountPerBinsCheckM.pl \\
  --checkm {checkm} \\
  --parsed_bins {parsed_bins} \\
  --contigs_abundance {contigs_abundance} \\
  --contigs_indir {contigs_indir} \\
  > {outfile} && \\
cat {outfile} | awk '{{print $1}}' | sort | uniq > {bins_list}""".format(
    checkm = checkm,
    parsed_bins = parsed_bins,
    contigs_abundance = contigs_abundance,
    contigs_indir = contigs_indir,
    outfile = outfile,
    bins_list = bins_list
  )

    return job

def generate_link_file(parsed_bins, annotations, bins_list, link):

    job = Job(
        [parsed_bins, parsed_bins, bins_list, annotations], 
        [link],
        [
            
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
linkBinsToGenes.pl \\
  --link {parsed_bins} \\
  --link2 {annotations} \\
  --selected_bins {bins_list} \\
  > {link}""".format(
    parsed_bins = parsed_bins,
    annotations = annotations,
    bins_list = bins_list,
    link = link
  )

    return job

def bins_feature_table(summarized_bins, parsed_bins, contigs_abundance, outfile):

    job = Job(
        [summarized_bins, parsed_bins, contigs_abundance], 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
binsFeatureTable.R \\
  -s {summarized_bins} \\
  -p {parsed_bins} \\
  -a {contigs_abundance} \\
  -o {outfile}""".format(
    summarized_bins = summarized_bins,
    parsed_bins = parsed_bins,
    contigs_abundance = contigs_abundance,
    outfile = outfile
  )

    return job

def cleanup(tmpdir):
        
    job = Job(
        [""], 
        [""],
        [
            ['memtime', 'module_memtime']
        ]
    )
    
    job.command="""
\\
rm  {tmpdir} -rf""".format(
    tmpdir = tmpdir
    )
    return job

def fastq_to_fasta(fastq, fasta):
    
    job = Job(
        [fastq],
        [fasta],
        [
            
            ['tools', 'module_tools']
        ]
    )
    job.command="""
fastqToFasta.pl \\
  --fastq {fastq} | gzip > {fasta}""".format(
    fastq = fastq,
    fasta = fasta
    )
        
    return job

def fastqsgz_to_fastas(fastqs_gz, fastas, done):
    
    job = Job(
        fastqs_gz,
        fastas + [done],
        [
            ['tools', 'module_tools']
        ]
    )
    job.command="""
fastqsGzToFastas.pl \\
  --num_threads {num_threads} \\
  --infiles {fastqs_gz} \\
  --outfiles {fastas} && \\
  touch {done}""".format(
    fastqs_gz = ",".join(fastqs_gz),
    fastas = ",".join(fastas),
    done = done,
    num_threads = config.param("fastqs_to_fastas", "num_threads", 1, "posint")
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
 -r {r} -f {f} \\
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
    taxonomy_folder = config.param("CAT", "taxonomy_folder", 1, "string"),
    r = config.param("CAT", "r", 1, "int"),
    f = config.param("CAT", "f", 1, "float")
    )

    return job

def generate_feature_table_consensus(infile_taxonomy, infile_abundance, outfile_taxonomy, outfile):

    job = Job(
        [infile_taxonomy, infile_abundance],
        [outfile, outfile_taxonomy], 
        [
            ['meco_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
    
    job.command="""
generateFeatureTableFromCAT.pl \\
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

def spades_mag(infile_paired_R1, infile_paired_R2, trusted_contigs, outdir):
    
    #temporary fix; change fna to fasta
    #trusted_contigs = os.path.splitext(trusted_contigs)[0]
    #trusted_contigs = trusted_contigs + ".fa"
    min_contig_length = str(config.param("spades_mag", "min_contig_length", 1, "posint"))

    job = Job(
        [infile_paired_R1, infile_paired_R2, trusted_contigs],
        [os.path.join(outdir, "scaffolds.fasta"), os.path.join(outdir, "scaffolds_gt" + min_contig_length+ ".fasta")],
        [
            ['spades', 'module_spades'],
            ['tools','module_tools']
        ]
    )
  
    job.command="""
rm {outdir} -rf && \\
spades.py --careful -t {num_threads} -k {kmers} \\
  -1 {infile_paired_R1} -2 {infile_paired_R2} \\
  --untrusted-contigs {trusted_contigs} \\
  -o {outdir} && \\
filterFastaByLength.pl \\
  --infile {outdir}/scaffolds.fasta \\
  --length {length} > \\
  {outdir}/scaffolds_gt{length}.fasta""".format(
    num_threads = config.param("spades_mag", "num_threads", 1, "posint"),
    kmers = config.param("spades_mag", "kmers", 1, "string"),
    length = config.param("spades_mag", "min_contig_length", 1, "posint"),
    infile_paired_R1 = infile_paired_R1,
    infile_paired_R2 = infile_paired_R2,
    trusted_contigs = trusted_contigs,
    outdir = outdir
  )

    return job

def opera_lg(contigs_file, bam_file, outdir):
     
    job = Job(
        [contigs_file, bam_file],
        [os.path.join(outdir, "scaffoldSeq.fasta")],
        [
            ['opera-lg', 'module_opera-lg'],
            ['samtools', 'module_samtools']
        ]
    )
    job.command="""
rm {outdir}/* -rf && \\
OPERA-LG \\
  {contigs_file} \\
  {bam_file} \\
  {outdir}""".format(
    contigs_file = contigs_file,
    bam_file = bam_file,
    outdir = outdir
  )

    return job

def sspace(contigs_file, fastq_R1, fastq_R2, outdir):
     
    job = Job(
        [contigs_file, fastq_R1, fastq_R2],
        [os.path.join(outdir, "standard_output.final.scaffolds.fasta")],
        [
            ['sspace', 'module_sspace'],
        ]
    )
    job.command="""
rm {outdir}/* -rf && \\
echo 'lib1 bowtie {fastq_R1} {fastq_R2} 200 0.25 FR' > {outdir}/libraries.txt && \\
SSPACE_Standard_v3.0.pl \\
  -s {contigs_file} \\
  -x 1 -T {num_threads} \\
  -l {outdir}/libraries.txt && \\
mv standard_output/* {outdir}/ && \\
rm -rf standard_output""".format(
    contigs_file = contigs_file,
    outdir = outdir,
    num_threads = config.param("sspace", "num_threads", 1, "posint"),
    fastq_R1 = fastq_R1,
    fastq_R2 = fastq_R2
  )

    return job

def gapfiller(fastq1, fastq2, contigs_file, outdir):
     
    job = Job(
        [fastq1, fastq2, contigs_file],
        [os.path.join(outdir, "GapFiller_output.fasta")],
        [
            ['gapfiller', 'module_gapfiller'],
            ['libboost','module_libboost']
        ]
    )
    job.command="""
rm {outdir}/* -rf && \\
GapFiller \\
  --seed1 {fastq1} \\
  --seed2 {fastq2} \\
  --seed-ins 250 \\
  --seed-var 100 \\
  --query {contigs_file} \\
  --output-prefix {outdir}/GapFiller_output""".format(
    contigs_file = contigs_file,
    fastq1 = fastq1,
    fastq2 = fastq2,
    outdir = outdir
  )

    return job

def baseclear_gapfiller(fastq1, fastq2, contigs_file, outdir, sample_name):
     
    job = Job(
        [fastq1, fastq2, contigs_file],
        [os.path.join(outdir, "outgf.gapfilled.final.fa")],
        [
            ['gapfiller', 'module_gapfiller']
            #['libboost','module_libboost']
        ]
    )
    job.command="""
rm {outdir}/* -rf && \\
echo 'lib1 bowtie {fastq1} {fastq2} 200 0.25 FR' > {outdir}/libraries.txt && \\
GapFiller.pl \\
  -l {outdir}/libraries.txt \\
  -s {contigs_file} \\
  -m 30 -o 2 -r 0.7 -n 10 -d 50 -t 10 -g 0 -i 1 \\
  -b {sample_name} \\
  -T {num_threads} && \\
mv {sample_name}/* {outdir}/ && \\
rm -rf {sample_name} && \\
mv {outdir}/{sample_name}.summaryfile.final.txt {outdir}/outgf.summaryfile.final.txt && \\
mv {outdir}/{sample_name}.gapfilled.final.fa {outdir}/outgf.gapfilled.final.fa && \\
mv {outdir}/{sample_name}.filled.final.txt {outdir}/outgf.filled.final.txt && \\
mv {outdir}/{sample_name}.closed.evidence.final.txt {outdir}/outgf.closed.evidence.final.txt""".format(
    contigs_file = contigs_file,
    fastq1 = fastq1,
    fastq2 = fastq2,
    outdir = outdir,
    sample_name = sample_name,
    num_threads = config.param("gapfiller", "num_threads", 1, "posint")
  )

    return job

def extract_silva_taxonomy(infile_blastn, abundance, accession_to_tax, outfile):

    job = Job(
        [infile_blastn, abundance, accession_to_tax],
        [outfile], 
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
  #--gi_to_tax {gi_to_tax} \\
    job.command="""
getSilvaTax.pl \\
  --infile_blastn {infile_blastn} \\
  --infile_accession {accession_to_tax} \\
  --infile_abundance {abundance} \\
  --evalue {evalue} \\
  --length {length} > {outfile}""".format(
    outfile = outfile,
    infile_blastn = infile_blastn,
    accession_to_tax = accession_to_tax,
    abundance = abundance,
    evalue = config.param('silva_tax', 'cutoff', required=True),
    length = config.param('silva_tax', 'alignment_length', required=True)
    )

    return job

def star_paired(infile_R1, infile_R2, infile_ref, index_dir, outdir, sample_id):
    
    outfile_mapped = os.path.join(outdir, sample_id + "_Aligned.sortedByCoord.out.bam")

    job = Job(
        [infile_R1, infile_R2, infile_ref, os.path.join(infile_ref)],
        [outfile_mapped], 
        [
            ['star', 'module_star']
        ]
    )
    
    job.command="""
STAR \\
  --runThreadN {num_threads} \\
  --genomeDir {index_dir} \\
  --readFilesIn {infile_R1} {infile_R2} \\
  --readFilesCommand zcat \\
  --outSAMtype BAM SortedByCoordinate \\
  --outFileNamePrefix {outdir}/{sample_id}_ \\
  --quantMode GeneCounts""".format(
      infile_R1 = infile_R1,
      infile_R2 = infile_R2,
      infile_ref = infile_ref,
      #outfile_mapped = outfile_mapped,
      outdir = outdir,
      index_dir = index_dir,
      sample_id = sample_id,
      num_threads = config.param('star', 'num_threads', 1, 'posint')
    )

    return job

def feature_count(infile_bam, outfile_counts):
    
    job = Job(
        [infile_bam],
        [outfile_counts], 
        [
            ['subread', 'module_subread']
        ]
    )
    
    job.command="""
featureCounts \\
  -p -B -C \\
  -T {num_threads} \\
  -a {gtf} \\
  -o {outfile_counts} \\
  {infile_bam}""".format(
      infile_bam = infile_bam,
      outfile_counts = outfile_counts,
      gtf = config.param('star', 'gtf', 1, 'filepath'),
      num_threads = config.param('feature_counts', 'num_threads', 1, 'posint')
    )

    return job

def bbmap_index(infile_ref, index_dir):
    
    job = Job(
        [infile_ref],
        [os.path.join(index_dir, "bbmap_index.done")], 
        [
            ['bbmap', 'module_bbmap']
        ]
    )
    
    job.command="""
bbmap.sh \\
  -threads={num_threads} \\
  -Xmx{ram} \\
  ref={infile_ref} \\
  path={index_dir} \\
&& touch {done_file}""".format(
      infile_ref = infile_ref,
      index_dir = index_dir,
      done_file = os.path.join(index_dir, "bbmap_index.done"),
      num_threads = config.param('bbmap', 'num_threads', 1, 'posint'),
      ram = config.param('bbmap', 'ram', 1, 'string')
    )

    return job

def bbmap_paired(infile_R1, infile_R2, infile_ref, outfile_mapped, index_dir):
    
    job = Job(
        [infile_R1, infile_R2, infile_ref, os.path.join(index_dir, "bbmap_index.done")],
        [outfile_mapped], 
        [
            ['bbmap', 'module_bbmap'],
            ['samtools', 'module_samtools'],
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )
    
    job.command="""
bbmap.sh \\
  -threads={num_threads} \\
  -Xmx{ram} \\
  in={infile_R1} in2={infile_R2} \\
  outm={outfile_mapped} \\
  ref={infile_ref} \\
  minid={minid} \\
  sam=1.3 nhtag=t mdtag=t xmtag=t amtag=t nmtag=t xstag=us \\
  path={index_dir} && \\
samtools view -Sbh -F 0x100 {outfile_mapped} > {outfile_mapped}.tmp && \\
  samtools sort -@ {num_threads} -m {mem_per_thread} {outfile_mapped}.tmp -o {outfile_mapped}.tmp.sorted.bam && \\
  mv {outfile_mapped}.tmp.sorted.bam {outfile_mapped} && \\
  rm {outfile_mapped}.tmp && \\
  samtools index {outfile_mapped}""".format(
      infile_R1 = infile_R1,
      infile_R2 = infile_R2,
      infile_ref = infile_ref,
      outfile_mapped = outfile_mapped,
      index_dir = index_dir,
      num_threads = config.param('bbmap', 'num_threads', 1, 'posint'),
      mem_per_thread = config.param('samtools', 'mem_per_thread', 1, 'string'),
      ram = config.param('bbmap', 'ram', 1, 'string'),
      minid = config.param('bbmap', 'min_id', 1, 'float')
    )

    return job

def bbmap_subtract(infile_R1, infile_R2, outfile_unmapped_R1, outfile_unmapped_R2, outfile_mapped_R1, outfile_mapped_R2, index_dir, log):
   
    #outdir = os.path.dirname(outfile_unmapped_R1) 

    job = Job(
        [infile_R1, infile_R2],
        [outfile_unmapped_R1, outfile_unmapped_R2, log], 
        [
            ['bbmap', 'module_bbmap']
        ]
    )
    
    job.command="""
bbmap.sh \\
  -threads={num_threads} \\
  -Xmx{ram} \\
  in={infile_R1} in2={infile_R2} \\
  outu={outfile_unmapped_R1} outu2={outfile_unmapped_R2} \\
  outm={outfile_mapped_R1} outm2={outfile_mapped_R2} \\
  minid={minid} \\
  sam=1.3 nhtag=t mdtag=t xmtag=t amtag=t nmtag=t xstag=us \\
  maxindel=3 bwr=0.16 bw=12 fast=t overwrite=t minhits=2 \\
  path={index_dir} \\
  statsfile={log}""".format(
      infile_R1 = infile_R1,
      infile_R2 = infile_R2,
      outfile_unmapped_R1 = outfile_unmapped_R1,
      outfile_unmapped_R2 = outfile_unmapped_R2,
      outfile_mapped_R1 = outfile_mapped_R1,
      outfile_mapped_R2 = outfile_mapped_R2,
      index_dir = index_dir,
      log = log,
      num_threads = config.param('bbmap_sub', 'num_threads', 1, 'posint'),
      ram = config.param('bbmap_sub', 'ram', 1, 'string'),
      minid = config.param('bbmap_sub', 'min_id', 1, 'float')
    )

    return job

def bbmap_mag(infiles1, infiles2, outfile_mapped, infile_ref, outfile_R1, outfile_R2):

    job = Job(
        infiles1 + infiles2 + [infile_ref],
        [outfile_mapped, outfile_R1, outfile_R2],
        [
            ['bbmap', 'module_bbmap'],
            ['samtools', 'module_samtools'],
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )
    job.command="""
bbmapMags.pl \\
  --infiles1 {infiles1} \\
  --infiles2 {infiles2} \\
  --min_id {min_id} \\
  --outfile_mapped {outfile_mapped} \\
  --reference {infile_ref} \\
  --ram {ram} \\
  --num_threads {num_threads} && \\
reformat.sh \\
  overwrite=true \\
  in={outfile_mapped} \\
  out1={outfile_R1} \\
  out2={outfile_R2} \\
  interleaved addcolon""".format(
    num_threads = config.param('bbmap_mag', 'num_threads', type='posint', required=True),
    min_id = config.param('bbmap_mag', 'min_id', type='float', required=True),
    ram = config.param('bbmap_mag', 'ram', type='string', required=True),
    infiles1 = ",".join(infiles1),
    infiles2 = ",".join(infiles2),
    infile_ref = infile_ref,
    outfile_mapped = outfile_mapped,
    outfile_R1 = outfile_R1,
    outfile_R2 = outfile_R2
    )

    return job

def quast(contigs_file, outdir):
    
    job = Job(
        [contigs_file],
        [os.path.join(outdir, "report.tsv")],
        [
            ['python3', 'module_python3'],
            ['quast', 'module_quast']
        ]
    )
    
    job.command="""
quast.py {contigs_file} \\
  -o {outdir}""".format(
    contigs_file = contigs_file,
    outdir = outdir
    )

    return job

def barrnap(contigs_file, kingdom, outfile_fna, outfile_tsv):
    
    job = Job(
        [contigs_file],
        [outfile_fna, outfile_tsv],
        [
            ['python', 'module_perl'],
            ['bedtoolsforbarrnap', 'module_bedtoolsforbarrnap'],
            ['hmmer', 'module_hmmer'],
            ['barrnap', 'module_barrnap']
        ]
    )
    
    job.command="""
barrnap \\
    --kingdom {kingdom} \\
    --outseq {outfile_fna} \\
    --thread {num_threads} \\
    {contigs_file} \\
    > {outfile_tsv}""".format(
    contigs_file = contigs_file,
    kingdom = kingdom,
    outfile_fna = outfile_fna,
    outfile_tsv = outfile_tsv,
    num_threads = config.param('barrnap', 'num_threads', 1, 'posint')
    )

    return job

def trnascanse_array_job(indir, inprefix, outdir, outprefix, prefix, infiles, outfiles, dones, job_scheduler):

    job_array_suffix = ""
    if job_scheduler == "sge":
        job_array_suffix = "SGE_TASK_ID"
    else:
        job_array_suffix = "SLURM_ARRAY_TASK_ID"

    infile = os.path.join(indir, inprefix + "\$" + job_array_suffix + "2")
    outfile = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.tsv")
    done = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.done")

    all_outfiles = []
    for i in range(0, len(dones)):
        all_outfiles.append(dones[i])
    for i in range(0, len(outfiles)):
        all_outfiles.append(outfiles[i])
    
    job = Job(
        infiles,
        all_outfiles,
        [
            ['trnascanse', 'module_trnascanse'],
            ['infernal', 'module_infernal'],
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command = ""
    
    if job_scheduler == "sge":
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 1  ))
        """.format(job_array_suffix = job_array_suffix)
    else:
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 0  ))
        """.format(job_array_suffix = job_array_suffix)
    
    job.command+="""
{job_array_suffix}2=\$(printf "%07d" \${job_array_suffix}2)
mkdir -p {outdir} && \\
if [[ -s {infile} ]] ; then
    tRNAscan-SE \\
     -o {outfile} \\
     -G \\
     {infile} \\
     --thread {num_threads} --forceow \\
    && touch {done} 
else
    touch {outfile}
fi ;""".format(
    job_array_suffix = job_array_suffix,
    infile = infile,
    outfile = outfile,
    infiles = infiles,
    outfiles = outfiles,
    outdir = outdir,
    num_threads = config.param(prefix, 'num_threads', 1, 'int'),
    done = done
    )

    return job

def split_barrnap_fna(infile_fasta, prefix, fiveS, SSU, LSU, twelveS, fiveeightS):

    outfiles = []
    if prefix == "bac" or prefix == "arc":
        outfiles = [fiveS, SSU, LSU]
    elif prefix == "euk":
        outfiles = [fiveS, SSU, LSU]
    elif prefix == "mito":
        outfiles = [SSU, twelveS]

    job = Job(
        [infile_fasta], 
        outfiles,
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
splitRrnaBarrnap.pl \\
  --infile_fasta {infile_fasta} \\
  --fiveS {fiveS} --SSU {SSU} --LSU {LSU} --twelveS {twelveS} --fiveeightS {fiveeightS}""".format(
    infile_fasta = infile_fasta,
    fiveS = fiveS,
    SSU = SSU,
    LSU = LSU,
    twelveS = twelveS,
    fiveeightS = fiveeightS
    )

    return job

def concat_fasta(infiles, outfile):

    job = Job(
        infiles, 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
cat {infiles_string} > {outfile}""".format(
    infiles_string = " ".join(infiles),
    outfile = outfile
    )

    return job

def rename_rrna_bed(infile, outfile):

    job = Job(
        [infile], 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
renameBed.pl --infile {infile} > {outfile}""".format(
    infile = infile,
    outfile = outfile
    )

    return job

def anvio_simplify_fasta_headers(infile, outfile):

    job = Job(
        [infile], 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
simplifyFastaHeaders.pl --infile {infile} > {outfile}""".format(
    infile = infile,
    outfile = outfile
    )

    return job

def anvio_gene_calling(infile_faa, infile_gff, infile_rrna, outfile):

    job = Job(
        [infile_faa, infile_gff, infile_rrna], 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
generateGeneCallingFileForAnvio.pl \\
  --infile_faa {infile_faa} \\
  --infile_gff {infile_gff} \\
  --infile_rrna_bac {infile_rrna} \\
  > {outfile}""".format(
    infile_faa = infile_faa,
    infile_gff = infile_gff,
    infile_rrna = infile_rrna,
    outfile = outfile
    )

    return job

def anvio_make_contigs_db(infile_fna, name, called_genes, outfile, done_file):

    job = Job(
        [infile_fna, called_genes], 
        [done_file],
        [
            ['anvio', 'module_anvio']
        ]
    )
    
    job.command="""
rm -rf {outfile} && \\
bash -c 'set +u && source \$ANVIO_VENV/activate && set -u && \\
  anvi-gen-contigs-database \\
    -f {infile_fna} \\
    -o {outfile} \\
    -n '{name}' \\
    --external-gene-calls {called_genes} \\
    --skip-predict-frame' && \\
touch {done_file}""".format(
    infile_fna = infile_fna,
    name = name,
    called_genes = called_genes,
    outfile = outfile,
    done_file = done_file
    )

    return job

def anvio_generate_kegg_file(infile_blast, infile_blast_parsed, outfile):

    job = Job(
        [infile_blast, infile_blast_parsed], 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
generateKeggFileForAnvio.pl \\
  --infile_kegg_diamond_blastp {infile_blast} \\
  --infile_kegg_annotations {infile_blast_parsed} \\
  > {outfile}""".format(
    infile_blast = infile_blast,
    infile_blast_parsed = infile_blast_parsed,
    outfile = outfile
    )

    return job

def anvio_generate_taxonomy_file(infile_gene_calls, infile_taxonomy, outfile):

    job = Job(
        [infile_gene_calls, infile_taxonomy], 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
generateTaxonomyFileForAnvio.pl \\
  --infile_gene_calls {infile_gene_calls} \\
  --infile_taxonomy {infile_taxonomy} \\
  > {outfile}""".format(
    infile_gene_calls = infile_gene_calls,
    infile_taxonomy = infile_taxonomy,
    outfile = outfile
    )

    return job

def anvio_generate_cog_file(infile_cog, outfile):

    job = Job(
        [infile_cog], 
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
generateCOGsFileForAnvio.pl \\
  --infile {infile_cog} \\
  > {outfile}""".format(
    infile_cog = infile_cog,
    outfile = outfile
    )

    return job

def anvio_import_annotations(infile_done, infile_cog, infile_taxonomy, infile_kegg, contigs_db, done_file):

    job = Job(
        [infile_cog, infile_taxonomy, infile_kegg, infile_done], 
        [done_file],
        [
            ['meco_tools', 'module_tools'],
            ['anvio', 'module_anvio']
        ]
    )
    
    job.command="""
bash -c 'set +u && source \$ANVIO_VENV/activate && set -u && \\
  anvi-import-taxonomy-for-genes -c {contigs_db} -i {infile_taxonomy} -p default_matrix && \\
  anvi-import-functions -c {contigs_db} -i {infile_kegg} && \\
  anvi-import-functions -c {contigs_db} -i {infile_cog}' && \\
touch {done_file}""".format(
    contigs_db = contigs_db,
    infile_kegg = infile_kegg,
    infile_cog = infile_cog,
    infile_taxonomy = infile_taxonomy,
    done_file = done_file
    )

    return job

def anvio_import_profiles(done_file, infile_bam, contigs_db, outdir, outfile):

    job = Job(
        [done_file, infile_bam],
        [outfile],
        [
            ['meco_tools', 'module_tools'],
            ['anvio', 'module_anvio']
        ]
    )
    
    job.command="""
bash -c 'set +u && source \$ANVIO_VENV/activate && set -u && \\
  anvi-profile -W -T {num_threads} -i {infile_bam} -c {contigs_db} -o {outdir}'""".format(
    contigs_db = contigs_db,
    done_file = done_file,
    infile_bam = infile_bam,
    outdir = outdir,
    num_threads = config.param('anvio_profile', 'num_threads', 1, 'posint')
    )

    return job

def anvio_merge(profiles, contigs_db, outdir, outfile_done):

    job = Job(
        profiles,
        [outfile_done, os.path.join(outdir, "PROFILE.db")],
        [
            ['meco_tools', 'module_tools'],
            ['anvio', 'module_anvio']
        ]
    )
    
    job.command="""
bash -c 'set +u && source \$ANVIO_VENV/activate && set -u && \\
  anvi-merge {profiles} \\
  -o {outdir} \\
  -c {contigs_db} \\
  --enforce-hierarchical-clustering \\
  --overwrite-output-destinations' && \\
  touch {outfile_done}""".format(
    contigs_db = contigs_db,
    outfile_done = outfile_done,
    outdir = outdir,
    profiles = " ".join(profiles)
    )

    return job

def anvio_run_hmm(infile_done, merged_profiles, contigs_db, done_file):

    job = Job(
        [infile_done],
        [done_file],
        [
            ['meco_tools', 'module_tools'],
            ['prodigal', 'module_prodigal'],
            ['hmmer', 'module_hmmer'],
            ['anvio', 'module_anvio']
        ]
    )
    
    job.command="""
bash -c 'set +u && source \$ANVIO_VENV/activate && set -u && \\
  anvi-run-hmms -T {num_threads} \\
  -c {contigs_db} \\
  --just-do-it' && \\
  touch {done_file}""".format(
    contigs_db = contigs_db,
    done_file = done_file,
    merged_profiles = merged_profiles,
    infile_done = infile_done,
    num_threads = config.param('anvio_hmms', 'num_threads', 1, 'posint')
    )

    return job

def anvio_import_metadata(infile_done, profiles_db, outfile_done):

    job = Job(
        [infile_done],
        [outfile_done],
        [
            ['meco_tools', 'module_tools'],
            ['anvio', 'module_anvio']
        ]
    )
    
    job.command="""
bash -c 'set +u && source \$ANVIO_VENV/activate && set -u && \\
  anvi-import-misc-data {mapping_file_for_anvio} \\
  -p {profiles_db} \\
  --target-data-table layers \\
  --just-do-it' && \\
  touch {outfile_done}""".format(
    infile_done = infile_done,
    outfile_done = outfile_done,
    profiles_db = profiles_db,
    mapping_file_for_anvio = config.param('DEFAULT', 'mapping_file_for_anvio', 1, 'filepath')
    )

    return job

def anvio_generate_mag_file(infile_done, indir, outfile):

    job = Job(
        [infile_done],
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
generateMAGFileForAnvio.pl --indir {indir} > {outfile}""".format(
    infile_done = infile_done,
    outfile = outfile,
    indir = indir
    )

    return job

def anvio_import_mags(infile_done, contigs_db, profiles_db, name, mag_file, outfile_done):

    job = Job(
        [infile_done, profiles_db, mag_file],
        [outfile_done],
        [
            ['meco_tools', 'module_tools'],
            ['anvio', 'module_anvio']
        ]
    )
    
    job.command="""
bash -c 'set +u && source \$ANVIO_VENV/activate && set -u && \\
  anvi-import-collection \\
  -C {name} \\
  -p {profiles_db} \\
  -c {contigs_db} \\
  --contigs-mode \\
  {mag_file}' && \\
  touch {outfile_done}""".format(
    infile_done = infile_done,
    outfile_done = outfile_done,
    profiles_db = profiles_db,
    contigs_db = contigs_db,
    mag_file = mag_file,
    name = name
    )

    return job

def cat_bins_into_single_file(infile_done, indir, outfile):

    job = Job(
        [infile_done],
        [outfile],
        [
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command="""
cat {indir}/*.fa > {outfile}""".format(
    infile_done = infile_done,
    outfile = outfile,
    indir = indir
    )

    return job

def calc_kmer_freq(infile, outfile):

    job = Job( 
        [infile],
        [outfile],
        [
            ['meco_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )    
    
    job.command="""
calcKmerFreq.pl \\
  -i {infile} \\
  -k 5 \\
  -m 1000 \\
  -o {outfile}""".format(
    infile = infile,
    outfile = outfile
    )    

    return job

def clr_transform(infile, outfile):

    job = Job( 
        [infile],
        [outfile],
        [
            ['meco_tools', 'module_tools'],
            ['R', 'module_R'],
            ['perl', 'module_perl']
        ]
    )    
    
    job.command="""
clrTransform.R \\
  -i {infile} \\
  -o {outfile}""".format(
    infile = infile,
    outfile = outfile
    )    

    return job


def bhtsne(infile, outfile):

    job = Job(
        [infile],
        [outfile],
        [
            ['bhtsne', 'module_bhtsne'],
            ['python3', 'module_python3']
        ]
    )
    
    job.command="""
bhtsne.py \\
  --use_pca \\
  -d {no_dim} \\
  -p {perplexity} \\
  -t {theta} \\
  -n {init_dims} \\
  -i {infile} \\
  -o {outfile}.tmp && \\
tail -n +2 {infile} | awk '{{print \$1}}' | paste - {outfile}.tmp > {outfile} && \\
rm {outfile}.tmp""".format(
    infile = infile,
    outfile = outfile,
    perplexity = config.param('bhtsne', 'perplexity', 1, 'posint'),
    theta = config.param('bhtsne', 'theta', 1, 'float'),
    init_dims = config.param('bhtsne', 'initial_dims', 1, 'posint'),
    no_dim = config.param('bhtsne', 'no_dim', 1, 'posint')
    )

    return job

def get_contigs_length_and_gc(infile, outfile):

    job = Job(
        [infile],
        [outfile],
        [
            ['meco_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
getFastaLength.pl \\
  --infile {infile} \\
  --gc \\
  > {outfile}""".format(
    infile = infile,
    outfile = outfile
    )

    return job

def parse_viral_contigs(infile, infile_fasta, infile_contigs_length_gc, infile_annotations, 
                        outfile, outfile_viral_contigs_list, outfile_viral_contigs_fasta):

    job = Job(
        [infile, infile_fasta, infile_contigs_length_gc, infile_annotations],
        [outfile, outfile_viral_contigs_list, outfile_viral_contigs_fasta],
        [
            ['meco_tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
parseViromeHmms.pl \\
  --length_cutoff {length_cutoff} \\
  --infile_contigs_length {infile_contigs_length_gc} \\
  --infile {infile} \\
  --infile_annotations {infile_annotations} \\
  > {outfile}.raw && \\
head -n 1 {outfile}.raw > {outfile} && \\
tail +2 {outfile}.raw | sort -n -r -k2,2 >> {outfile} && \\
cat {outfile} | awk '\$2 >= 5' | awk '\$8 <= 20' | awk '\$6 <= 40' | awk '\$4 >= 10' | cut -f 1 > {outfile}.filter1 && \\
cat {outfile} | awk '\$2 >= 5' | awk '\$2 >= \$5' | cut -f 1 > {outfile}.filter2 && \\
cat {outfile} | awk '\$2 >= 5' | awk '\$4 >= 60' | cut -f 1 > {outfile}.filter3 && \\
cat {outfile}.filter* | sort | uniq > {outfile_viral_contigs_list} && \\
getFastaSeqsByListOfGenes.pl --gene_list {outfile_viral_contigs_list} --infile {infile_fasta} > {outfile_viral_contigs_fasta}""".format(
    infile = infile,
    infile_annotations = infile_annotations,
    infile_fasta = infile_fasta,
    infile_contigs_length_gc = infile_contigs_length_gc,
    outfile = outfile,
    outfile_viral_contigs_list = outfile_viral_contigs_list,
    outfile_viral_contigs_fasta = outfile_viral_contigs_fasta,
    length_cutoff = config.param("virome", "length_cutoff", 1, "posint")
    )

    return job

def kofamscan_array_job(indir, inprefix, outdir, outprefix, prefix, infiles, outfiles, dones, job_scheduler=""):

    job_array_suffix = ""
    if job_scheduler == "sge":
        job_array_suffix = "SGE_TASK_ID"
    else:
        job_array_suffix = "SLURM_ARRAY_TASK_ID"

    infile = os.path.join(indir, inprefix + "\$" + job_array_suffix + "2")
    outfile = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.tsv")
    done = os.path.join(outdir, outprefix + "\$" + job_array_suffix + "2.done")

    all_outfiles = []
    for i in range(0, len(dones)):
        all_outfiles.append(dones[i])
    for i in range(0, len(outfiles)):
        all_outfiles.append(outfiles[i])

    job = Job(
        infiles,
        all_outfiles,
        [
            ['kofamscan', 'module_kofamscan'],
            ['hmmer', 'module_hmmer'],
            ['parallel', 'module_parallel'],
            ['meco_tools', 'module_tools']
        ]
    )
    
    job.command = ""
    #sys.stderr.write("Job scheduler: " + job_scheduler + "\n")
    if job_scheduler == "sge":
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 1  ))
        """.format(job_array_suffix = job_array_suffix)
    else:
        job.command="""
{job_array_suffix}2=\$(( {job_array_suffix} - 0  ))
        """.format(job_array_suffix = job_array_suffix)
    
    job.command+="""
{job_array_suffix}2=\$(printf "%07d" \${job_array_suffix}2)
echo 'SLURM/SGE_TASK_ID2' \${job_array_suffix}2
if [[ -s {infile} ]] ; then
exec_annotation \\
    -o {outfile} \\
    -p {profiles_dir} \\
    -k {ko_list} \\
    --cpu {num_threads} \\
    --tmp-dir {outdir}/\${job_array_suffix}2 \\
    -E {evalue} \\
    --keep-tabular \\
    {infile} \\
    && touch {done} \\
    && rm -rf {outdir}/\${job_array_suffix}2
else
    touch {outfile}
fi ;""".format(
    job_array_suffix = job_array_suffix,
    infile = infile,
    outfile = outfile,
    evalue = config.param('kofamscan', 'evalue', 1, 'string'),
    num_threads = config.param('kofamscan', 'num_threads', 1, 'int'),
    ko_list = config.param('kofamscan', 'ko_list', 1, 'filepath'),
    #tmpdir = config.param('DEFAULT', 'tmpdir', 1, 'string'),
    profiles_dir = config.param('kofamscan', 'profiles', 1, 'dirpath'),
    outdir = outdir,
    done = done
    )

    return job
