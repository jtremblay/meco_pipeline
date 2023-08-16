#!/usr/bin/env python

# Python Standard Modules

# MECO Modules
from core.config import *
from core.job import *

def validate_map_and_barcodes(barcodes, outfile_done):
    
    job = Job(
        [barcodes], 
        [outfile_done],
        [['tools', 'module_tools']]
    )
    
    job.command = """
validateMapAndBarcodes.pl \\
  --infile_barcodes {barcodes} \\
  --infile_mapping_file {map} && \\
  touch {outfile_done}""".format(
        barcodes = barcodes,
        map = config.param('default', 'mapping_file', 1, 'filepath'),
        outfile_done = outfile_done
    )

    return job

def bbduk(infile,
          contam, 
          ncontam, 
          ncontam_single, 
          log, 
          db, 
          infile_done=False
         ):

    if(infile_done == False):
        infiles = [infile]
    else:
        infiles = [infile, infile_done]

    job = Job(
        infiles, 
        [contam, ncontam, log, ncontam_single],
        [
            ['bbmap', 'module_bbmap']
        ]
    )
        
    job.command="""
bbduk.sh \\
  interleaved=1 \\
  in={infile} \\
  stats={log} \\
  out={ncontam} \\
  outm={contam} \\
  outs={ncontam_single} \\
  k={k} \\
  minkmerhits={c} \\
  ref={db} \\
  overwrite=true \\
  threads=1""".format(
    infile = infile,
    log = log,
    ncontam = ncontam,
    ncontam_single = ncontam_single,
    contam = contam,
    k = config.param('bbduk', 'k', 'int'),
    c = config.param('bbduk', 'c', 'int'),
    db = db
    ) 
    return job

def split_barcodes(infile, barcodes, outfile, log):
    
    job = Job(
        [infile],
        [outfile],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
barcodes.pl \\
--infile {infile} \\
--barcodes {barcodes} \\
--outfile {outfile} \\
--num_threads {num_threads}  \\
--log {log}""".format(
    infile = infile,
    outfile = outfile,
    barcodes = barcodes,
    num_threads = config.param( 'barcodes', 'num_threads', 1, 'int'),
    log = log
    )
    
    return job

def split_barcodes_dir(infile, barcodes, outfile, log, outdir, suffix=None, rm_old_files=True):
    
    job = Job(
        [infile],
        [outfile],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
#rm -rf {outdir}/* && \\
    job.command=""""""
    if(rm_old_files is True): 
        job.command+="""rm -rf {outdir}/*{suffix}.fastq && \\""".format(
            outdir = outdir,
            suffix = suffix
        )
    if(suffix is not None): 
        job.command+="""
barcodes.pl \\
--infile {infile} \\
--barcodes {barcodes} \\
--outdir {outdir} \\
--num_threads {num_threads} \\
--suffix {suffix} \\
--log {log} && \\
touch {outfile}""".format(
    infile = infile,
    outfile = outfile,
    barcodes = barcodes,
    outdir = outdir,
    suffix = suffix,
    num_threads = config.param( 'barcodes', 'num_threads', 1, 'int'),
    #trim_length = config.param( 'barcodes', 'trim_length', 1, 'int'),
    log = log
    )
    else:
        job.command+="""
barcodes.pl \\
--infile {infile} \\
--barcodes {barcodes} \\
--outdir {outdir} \\
--num_threads {num_threads} \\
--log {log} && \\
touch {outfile}""".format(
    infile = infile,
    outfile = outfile,
    barcodes = barcodes,
    outdir = outdir,
    suffix = suffix,
    num_threads = config.param( 'barcodes', 'num_threads', 1, 'int'),
    #trim_length = config.param( 'barcodes', 'trim_length', 1, 'int'),
    log = log
    )

    
    return job

def split_barcodes_dir_sorted(infile, barcodes, outdir, log):
#$barcodes_splitter_tool." --infile_fastq ".$illumina_infile." --infile_barcodes ".$barcodes_infile." --outdir ".$tmpdir_fastq 
    job = Job(
        [infile],
        [log],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
barcodesSorted.pl \\
 --infile_fastq {infile} \\
 --infile_barcodes {barcodes} \\
 --outdir {outdir} \\
 --log {log}""".format(
    infile = infile,
    log = log,
    barcodes = barcodes,
    outdir = outdir
    )
    
    return job

def remove_unpaired_reads(infile, outfile_paired, unpairedR1, unpairedR2):
    job = Job(
        [infile], 
        [outfile_paired],
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
--num_threads {num_threads}""".format(
    infile = infile,
    outfile_paired = outfile_paired,
    unpairedR1 = unpairedR1,
    unpairedR2 = unpairedR2,
    num_threads =  config.param('remove_unpaired', 'num_threads', 'int')
    )
                
    return job

def split_pairs(infile, outfileR1, outfileR2):
        
    job = Job(
        [infile], 
        [outfileR1, outfileR2],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
splitPairs.pl \\
--infile {infile} \\
--outfile_1 {outfileR1} \\
--outfile_2 {outfileR2} \\
--num_threads {num_threads}""".format(
    infile = infile,
    outfileR1 = outfileR1,
    outfileR2 = outfileR2,
    num_threads = config.param( 'split_pairs', 'num_threads', 'int')
    )
                
    return job

def generate_qscore_sheet(infile, prefix, log, outfile, barcodes):
    suffix = "suffix"

    job = Job(
        [infile],
        [outfile],
        [
            ['fastx', 'module_fastx'],
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
qscoreSheets.pl \\
 --fastq {infile} \\
 --tmp {tmp} \\
 --prefix {prefix} \\
 --suffix {suffix} \\
 --log {log} \\
 --outfile {outfile} \\
 --phred {phred} \\
 --barcodes {barcodes} \\
 --num_threads {num_threads}""".format(
        infile = infile,
        suffix = suffix,
        prefix = prefix,
        barcodes = barcodes,
        log = log,
        outfile = outfile,
        tmp = config.param('default', 'tmpdir', 1, 'dirpath'),
        phred = config.param('default', 'qual', 1, 'int'),
        num_threads = config.param('qscore_sheet', 'num_threads', 1, 'int')
    )
                    
    return job

def generate_qscore_graph_single(infile, prefix, outfile):
        
    job = Job(
        [infile] , 
        [outfile],
        [
            ['R', 'module_R'],
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
    job.command="""
qscorePlots.pl \\
 --infile_1 {infile} \\
 --name {prefix} \\
 --pdf {outfile} \\
 --display 1 \\
 --single""".format(
    infile = infile,
    prefix = prefix,
    outfile = outfile
    )

    return job

def generate_qscore_graph_paired(infileR1, infileR2, outfile):
                
    job = Job(
        [infileR1, infileR2], 
        [outfile],
        [
            ['R', 'module_R'],
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
qscorePlots.pl \\
  --infile_1 {infileR1} \\
  --infile_2 {infileR2} \\
  --name qual_stats \\
  --pdf {outfile} \\
  --display 1 \\
  --paired""".format(
    infileR1 = infileR1,
    infileR2 = infileR2,
    outfile = outfile
    )

    return job

def cut_reads(infile, begin, end, outfile):
        
    job = Job(
        [infile],
        [outfile],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
cutFastqSeq.pl \\
--infile {infile} \\
--begin {begin} \\
--end {end} \\
--outfile {outfile}""".format(
    infile = infile,
    begin = begin,
    end = end,
    outfile = outfile
    )

    return job

def flash(infileR1, infileR2, prefix, outdir):
                
    job = Job(
        [infileR1, infileR2], 
        [outdir + "/assembly_complete/ncontam_nphix_trimmed.extendedFrags.fastq"],
        [
            ['flash', 'module_flash'],
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
    job.command="""
flash.pl \\
--infile_1 {infileR1} \\
--infile_2 {infileR2} \\
--prefix {prefix} \\
--outdir {outdir} \\
--n {n} \\
--mismatches {m} \\
--p {p} \\
--num_threads {num_threads}""".format(
    infileR1 = infileR1,
    infileR2 = infileR2,
    prefix = prefix,
    outdir = outdir,
    n = config.param('flash', 'sampling', 'int'),
    m = config.param('flash', 'minOverlap', 'int'),
    p = config.param('flash', 'phred', 'int'),
    num_threads = config.param( 'flash', 'num_threads', 'int')
    )
    return job

def old_tags_qc(infile, fwd_primer, rev_primer, outfile, outfile_failed, log, reject_unmatched=False, correct_orientation=False):
        
    job = Job(
        [infile],
        [outfile, outfile_failed, log],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
tagsQC.pl \\
 --infile {infile} \\""".format(
     infile = infile
    )

    if(fwd_primer != "null"):
        job.command+="""
 --primer_5_prime {fwd_primer} \\
 --length_5_prime {length_5_prime} \\""".format(
            fwd_primer = fwd_primer,
            length_5_prime = config.param( 'tags_QC', 'length5Prime', 'int')
        )
    
    if(rev_primer != "null"):
        job.command+="""
 --primer_3_prime {rev_primer} \\
 --length_3_prime {length_3_prime} \\""".format(
            rev_primer = rev_primer,
            length_3_prime = config.param('tags_QC', 'length3Prime', 'int') 
        )
    if reject_unmatched is True or config.param('tags_QC', 'reject_unmatched', 'string') == "true":
        job.command+="""
 --reject_unmatched \\"""
    
    if correct_orientation is True:
        job.command+="""
 --correct_orientation \\"""

    job.command+="""
 --qscore_1 {qscore1} \\
 --qscore_2 {qscore2} \\
 --outfile {outfile} \\
 --outfile_failed {outfile_failed} \\
 --num_threads {num_threads} \\
 --qual {qual} \\
 --lq_threshold {lq_threshold} \\
 --primer_mismatch {primer_mismatch} \\
 --min_length {min_length} \\
 --max_length {max_length} \\
 --N {N} \\
 --log {log}""".format(
        qscore1 = config.param('tags_QC', 'qscore1', 'int'),
        qscore2 = config.param('tags_QC', 'qscore2', 'int'),
        outfile = outfile,
        outfile_failed = outfile_failed,
        num_threads =  config.param('tags_QC', 'num_threads', 'int'),
        qual = config.param('default', 'qual', 'int'),
        lq_threshold = config.param('tags_QC', 'lq_threshold', 'int'),
        primer_mismatch = config.param( 'tags_QC', 'primerMismatch', 'float'),
        min_length = config.param('tags_QC', 'minlength', 'int'),
        max_length = config.param('tags_QC', 'maxlength', 'int'),
        N = config.param('tags_QC', 'N', 'int'),
        log = log
    )

                
    return job

def legacy_tags_qc_R1R2(infile_R1, infile_R2, fwd_primer, rev_primer, 
                 outfile_R1, outfile_R2, outfile_failed_R1, outfile_failed_R2, 
                 log, reject_unmatched=False):
        
    job = Job(
        [infile_R1, infile_R2],
        [outfile_R1, outfile_R2, outfile_failed_R1, outfile_failed_R2, log],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
tagsQC.pl \\
 --infile {infile_R1} \\
 --infile_R2 {infile_R2} \\""".format(
     infile_R1 = infile_R1,
     infile_R2 = infile_R2
    )

    if(fwd_primer != "null"):
        job.command+="""
 --primer_5_prime {fwd_primer} \\
 --length_5_prime {length_5_prime} \\""".format(
            fwd_primer = fwd_primer,
            length_5_prime = config.param( 'tags_QC', 'length5Prime', 'int')
        )

    if(rev_primer != "null"):
        job.command+="""
 --primer_3_prime {rev_primer} \\
 --length_3_prime {length_3_prime} \\""".format(
            rev_primer = rev_primer,
            length_3_prime = config.param('tags_QC', 'length3Prime', 'int') 
        )
    
    if reject_unmatched is True or config.param('tags_QC', 'reject_unmatched', 'string') == "true":
        job.command+="""
 --reject_unmatched \\"""

    job.command+="""
 --qscore_1 {qscore1} \\
 --qscore_2 {qscore2} \\
 --outfile {outfile_R1} \\
 --outfile_failed {outfile_failed_R1} \\
 --outfile_R2 {outfile_R2} \\
 --outfile_failed_R2 {outfile_failed_R2} \\
 --num_threads {num_threads} \\
 --qual {qual} \\
 --log {log} \\
 --lq_threshold {lq_threshold} \\
 --primer_mismatch {primer_mismatch} \\
 --min_length {min_length} \\
 --max_length {max_length} \\
 --N {N}""".format(
        qscore1 = config.param('tags_QC', 'qscore1', 'int'),
        qscore2 = config.param('tags_QC', 'qscore2', 'int'),
        outfile_R1 = outfile_R1,
        outfile_failed_R1 = outfile_failed_R1,
        outfile_R2 = outfile_R2,
        outfile_failed_R2 = outfile_failed_R2,
        num_threads =  config.param('tags_QC', 'num_threads', 'int'),
        qual = config.param('default', 'qual', 'int'),
        lq_threshold = config.param('tags_QC', 'lq_threshold', 'int'),
        primer_mismatch = config.param( 'tags_QC', 'primerMismatch', 'float'),
        min_length = config.param('tags_QC', 'minlength', 'int'),
        max_length = config.param('tags_QC', 'maxlength', 'int'),
        N = config.param('tags_QC', 'N', 'int'),
        log = log
    )
                
    return job

def bbduk_qc_R1R2(infile_R1, infile_R2, 
                  outfile_R1_noN, outfile_R2_noN, 
                  fwd_and_rev_primer_file, 
                  outfile_R1, outfile_R2, 
                  outfile_failed_R1, outfile_failed_R2, 
                  log):
        
    job = Job(
        [infile_R1, infile_R2],
        [outfile_R1_noN, outfile_R2_noN, outfile_R1, outfile_R2, outfile_failed_R1, outfile_failed_R2, log],
        [
            ['bbmap', 'module_bbmap'],
            ['java', 'module_java']
        ]
    )

    job.command="""
reformat.sh \\
   in={infile_R1} \\
   out={outfile_R1_noN} \\
   qtrim=rl \\
   overwrite=true \\
   trimq=1 && \\
reformat.sh \\
   in={infile_R2} \\
   out={outfile_R2_noN} \\
   qtrim=rl \\
   overwrite=true \\
   trimq=1 && \\
bbduk.sh \\
    in={outfile_R1_noN} \\
    in2={outfile_R2_noN} \\
    ref={fwd_and_rev_primer_file} \\
    out={outfile_R1} \\
    out2={outfile_R2} \\
    outm={outfile_failed_R1} \\
    outm2={outfile_failed_R2} \\
    ktrim=l \\
    hdist={hdist} \\
    k={k} \\
    tbo=true \\
    overwrite=true \\
    restrictleft={restrictleft} \\
    maq={maq} \\
    stats={log} \\
    minlength={minlength} \\
    maxlength={maxlength} \\
    threads={num_threads}""".format(
        maq = config.param('bbduk_QC', 'maq', 'int'),
        infile_R1 = infile_R1,
        outfile_R1 = outfile_R1,
        outfile_R1_noN = outfile_R1_noN,
        outfile_failed_R1 = outfile_failed_R1,
        infile_R2 = infile_R2,
        outfile_R2 = outfile_R2,
        outfile_R2_noN = outfile_R2_noN,
        outfile_failed_R2 = outfile_failed_R2,
        num_threads =  config.param('bbduk_QC', 'num_threads', 'int'),
        hdist = config.param( 'bbduk_QC', 'hdist', 'int'),
        restrictleft = config.param( 'bbduk_QC', 'restrictleft', 'int'),
        k = config.param( 'bbduk_QC', 'k', 'int'),
        log = log,
        minlength = config.param('bbduk_QC', 'minlength', 'posint'),
        maxlength = config.param('bbduk_QC', 'maxlength', 'posint'),
        fwd_and_rev_primer_file = fwd_and_rev_primer_file
    )
                
    return job

def tags_qc_R1R2(infile_R1, infile_R2, 
                  outfile_R1_noN, outfile_R2_noN, 
                  primer_file, 
                  outfile_R1, outfile_R2, 
                  log):
        
    job = Job(
        [infile_R1, infile_R2],
        [outfile_R1_noN, outfile_R2_noN, outfile_R1, outfile_R2, log],
        [
            ['bbmap', 'module_bbmap'],
            ['java', 'module_java'],
            ['ptrimmer', 'module_ptrimmer']
        ]
    )

#reformat.sh \\
#   in={infile_R2} \\
#   out={outfile_R2_noN} \\
#   maxns={maxNs} \\
#   qtrim=rl \\
#   overwrite=true \\
#   maxlength={maxlength} \\
#   minlength={minlength} \\
#   trimq=1 && \\
    job.command="""
reformat.sh \\
   in={infile_R1} \\
   in2={infile_R2} \\
   out={outfile_R1_noN} \\
   out2={outfile_R2_noN} \\
   maxns={maxNs} \\
   qtrim=rl \\
   overwrite=true \\
   maxlength={maxlength} \\
   minlength={minlength} \\
   trimq=1 && \\
pTrimmer \\
   -t pair \\
   -a {primer_file} \\
   -f {outfile_R1_noN} \\
   -d {outfile_R1} \\
   -r {outfile_R2_noN} \\
   -e {outfile_R2} \\
   -k {k} \\
   -m {m} \\
   -q {maq} \\
   -s {log}""".format(
        infile_R1 = infile_R1,
        outfile_R1 = outfile_R1,
        outfile_R1_noN = outfile_R1_noN,
        infile_R2 = infile_R2,
        outfile_R2 = outfile_R2,
        outfile_R2_noN = outfile_R2_noN,
        m = config.param( 'separate_reads_qc', 'm', 1, 'int'),
        k = config.param( 'separate_reads_qc', 'k', 1, 'int'),
        maq = config.param('separate_reads_qc', 'maq', 1, 'int'),
        minlength = config.param('separate_reads_qc', 'minlength', 1, 'posint'),
        maxlength = config.param('separate_reads_qc', 'maxlength', 1, 'posint'),
        maxNs = config.param( 'separate_reads_qc', 'maxNs', 1, 'int'),
        log = log,
        primer_file = primer_file
    )
                
    return job

def tags_qc_se(infile_R1, 
               outfile_R1_noN, 
               primer_file, 
               outfile_R1, 
               log):
        
    job = Job(
        [infile_R1, primer_file],
        [outfile_R1_noN, outfile_R1, log],
        [
            ['bbmap', 'module_bbmap'],
            ['java', 'module_java'],
            ['ptrimmer', 'module_ptrimmer']
        ]
    )

    job.command="""
reformat.sh \\
   in={infile_R1} \\
   out={outfile_R1_noN} \\
   maxns={maxNs} \\
   qtrim=rl \\
   overwrite=true \\
   maxlength={maxlength} \\
   minlength={minlength} \\
   trimq=1 && \\
pTrimmer \\
   -t single \\
   -a {primer_file} \\
   -f {outfile_R1_noN} \\
   -d {outfile_R1} \\
   -k {k} \\
   -m {m} \\
   -q {maq} \\
   -s {log}""".format(
        infile_R1 = infile_R1,
        outfile_R1 = outfile_R1,
        outfile_R1_noN = outfile_R1_noN,
        m = config.param( 'separate_reads_qc', 'm', 1, 'int'),
        k = config.param( 'separate_reads_qc', 'k', 1, 'int'),
        maq = config.param('separate_reads_qc', 'maq', 1, 'int'),
        minlength = config.param('separate_reads_qc', 'minlength', 1, 'posint'),
        maxlength = config.param('separate_reads_qc', 'maxlength', 1, 'posint'),
        maxNs = config.param( 'separate_reads_qc', 'maxNs', 1, 'int'),
        log = log,
        primer_file = primer_file
    )
                
    return job

def tags_qc_ex(infile_R1, infile_R2, 
               outfile_R1_noN, outfile_R2_noN, 
               primer_file, 
               outfile_R1, outfile_R2, 
               log):
        
    job = Job(
        [infile_R1, infile_R2],
        [outfile_R1_noN, outfile_R2_noN, outfile_R1, outfile_R2, log],
        [
            ['bbmap', 'module_bbmap'],
            ['java', 'module_java'],
            ['ptrimmer', 'module_ptrimmer']
        ]
    )

    job.command="""
reformat.sh \\
   in={infile_R1} \\
   in2={infile_R2} \\
   out={outfile_R1_noN} \\
   out2={outfile_R2_noN} \\
   overwrite=true \\
   maxns=20 \\
   qtrim=rl \\
   maxns=20 \\
   trimq=1 && \\
pTrimmer \\
   -t pair \\
   -a {primer_file} \\
   -f {outfile_R1_noN} \\
   -d {outfile_R1} \\
   -r {outfile_R2_noN} \\
   -e {outfile_R2} \\
   -k {k} \\
   -m {m} \\
   -q 1 \\
   -s {log}""".format(
        infile_R1 = infile_R1,
        outfile_R1 = outfile_R1,
        outfile_R1_noN = outfile_R1_noN,
        infile_R2 = infile_R2,
        outfile_R2 = outfile_R2,
        outfile_R2_noN = outfile_R2_noN,
        m = config.param( 'separate_reads_qc', 'm', 1, 'int'),
        k = config.param( 'separate_reads_qc', 'k', 1, 'int'),
        maxNs = config.param( 'separate_reads_qc', 'maxNs', 1, 'int'),
        log = log,
        primer_file = primer_file
    )
                
    return job

def tags_qc_ex_merged(infile, outfile, log):
        
    job = Job(
        [infile],
        [outfile, log],
        [
            ['bbmap', 'module_bbmap'],
            ['java', 'module_java']
        ]
    )

    job.command="""
reformat.sh \\
   in={infile} \\
   out={outfile} \\
   maxns={maxns} \\
   overwrite=true \\
   minavgquality={maq} 2> {log}""".format(
        infile = infile,
        outfile = outfile,
        maq = config.param( 'merged_reads_qc', 'maq', 1, 'int'),
        maxns = config.param( 'merged_reads_qc', 'maxNs', 1, 'int'),
        log = log
    )
                
    return job

def count_report(rA_files, rA_names, analysis_type, barcodes_dist, feature_table, obs_table, outfile):
        
    job = Job(
        [feature_table],
        [outfile],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
    
    cmd = "countReport.pl \\\n"

    for file in rA_files:
        cmd += " --file " + file + " \\\n"

    for name in rA_names:
        cmd += " --name " + name + " \\\n"
    
    cmd += " --analysisType " + analysis_type + " \\\n"
    cmd += " --barcodesDist " + barcodes_dist + " \\\n"
    cmd += " --OTUtable " + feature_table + " \\\n"
    cmd += " --obsTable " + obs_table + " \\\n"
    cmd += " > " + outfile
     
    job.command = cmd
                    
    return job

def clustering_dnaclust(infile, barcodes, outdir, reads_type):
       
    script = ""
    if reads_type == "short_reads":
        script = "clusteringShortReadsDnaclust.pl"
    elif reads_type == "long_reads":
        script = "clusteringLongReadsDnaclust.pl"
    else:
        sys.stderr.write("[clustering] reads_type = short_reads OR long_reads\n")

    job = Job(
        [infile],
        [os.path.join(outdir, "obs.fasta"), os.path.join(outdir, "obs.tsv")],
        [
            ['vsearch', 'module_vsearch'],
            ['dnaclust', 'module_dnaclust'],
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
{script} \\
 --infile_fastq {infile} \\
 --ref_db {ref_db} \\
 --barcodes {barcodes} \\
 --outdir {outdir} \\
 --first_round {first_round} \\
 --second_round {second_round} \\
 --lowAbunCutOff {lowAbunCutOff} \\
 --num_threads {num_threads}""".format(
    infile = infile,
    ref_db =  config.param('DB', 'chimeras', 1, 'filepath'),
    barcodes = barcodes,
    outdir = outdir,
    lowAbunCutOff = config.param('clustering', 'min_reads', 'int'),
    first_round = config.param('clustering', 'first_round', 'float'),
    second_round = config.param('clustering', 'second_round', 'float'),
    num_threads =  config.param('clustering', 'num_threads', 'int'),
    script = script
    )
    
    return job

def clustering_vsearch(infile, barcodes, outdir, reads_type):
       
    script = ""
    if reads_type == "short_reads":
        script = "clusteringShortReadsVsearch.pl"
    elif reads_type == "long_reads":
        script = "clusteringLongReadsVsearch.pl"
    else:
        sys.stderr.write("[clustering] reads_type = short_reads OR long_reads\n")

    job = Job(
        [infile],
        [outdir + "/obs.fasta", outdir + "/obs.tsv"],
        [
            ['vsearch', 'module_vsearch'],
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""
{script} \\
 --infile_fastq {infile} \\
 --ref_db {ref_db} \\
 --barcodes {barcodes} \\
 --outdir {outdir} \\
 --first_round {first_round} \\
 --second_round {second_round} \\
 --lowAbunCutOff {lowAbunCutOff} \\
 --num_threads {num_threads}""".format(
    infile = infile,
    ref_db =  config.param('DB', 'chimeras', 1, 'filepath'),
    barcodes = barcodes,
    outdir = outdir,
    lowAbunCutOff = config.param('clustering', 'min_reads', 'int'),
    first_round = config.param('clustering', 'first_round', 'float'),
    second_round = config.param('clustering', 'second_round', 'float'),
    num_threads =  config.param('clustering', 'num_threads', 'int'),
    script = script
    )
    
    return job

def clustering_deblur(infile, indir, outdir, outfile_biom, outfile_fasta, done_file):
       

    job = Job(
        [infile],
        [outfile_biom, outfile_fasta, done_file],
        [
            ['deblur', 'module_deblur'],
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

#deblur workflow --seqs-fp ./infiles/ --output-dir ./deblur_outdir/  --overwrite -t 100
    job.command="""
deblur workflow --seqs-fp {indir} \\
 --overwrite \\
 -t {trim_length} \\
 --output-dir {outdir} \\
 -O {num_threads}  --indel-max {indel_max}  --min-reads {min_reads} --min-size {min_size} && \\
 touch {done_file}""".format(
    indir = indir,
    outdir = outdir,
    num_threads = config.param('clustering', 'num_threads', 1, 'int'),
    trim_length = config.param('clustering', 'trim_length', 1, 'int'),
    indel_max = config.param('clustering', 'indel_max', 1, 'int'),
    min_reads = config.param('clustering', 'min_reads', 1, 'int'),
    min_size = config.param('clustering', 'min_size', 1, 'int'),
    done_file = done_file
 )
    
    return job

def clustering_dada2(infiles_done, indir, outdir, outfile, done_file, data_type):

    job = Job(
        infiles_done,
        [outfile, done_file],
        [
            ['R', 'module_R'],
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )
    
    min_overlap = config.param('clustering', 'min_overlap', 0, 'int')
    if not isinstance(min_overlap, int):
        min_overlap = 10

    job.command="""
mkdir -p {outdir} && \\
dada2.R -i {indir} \\
 -o {outdir} \\
 -t {data_type} \\
 -p {num_threads} \\
 -l {trimLeft} \\
 -m {maxEE} \\
 -c {truncQ} \\
 -v {min_overlap} \\
 -q {minQ} &> {outdir}/log.txt && \\
 cat {outdir}/log.txt | parseDADA2Log.pl > {outdir}/log_merged_pairs.txt && \\
 touch {done_file}""".format(
    indir = indir,
    outdir = outdir,
    data_type = data_type,
    num_threads = config.param('clustering', 'num_threads', 1, 'int'),
    trimLeft = config.param('clustering', 'trimLeft', 1, 'int'),
    maxEE = config.param('clustering', 'maxEE', 1, 'int'),
    truncQ = config.param('clustering', 'truncQ', 1, 'int'),
    minQ = config.param('clustering', 'minQ', 1, 'int'),
    min_overlap = min_overlap,
    done_file = done_file
 )
    
    return job


def postprocess_asvs(done_file, dada2_tsv, outfile_tsv, outfile_fasta):
       
    job = Job(
        [dada2_tsv, done_file],
        [outfile_tsv, outfile_fasta],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl']
        ]
    )

#deblur workflow --seqs-fp ./infiles/ --output-dir ./deblur_outdir/  --overwrite -t 100
    job.command="""
convertASVIds.pl \\
  --min_reads {min_reads} \\
  --infile_tsv {dada2_tsv} \\
  --outfile_tsv {outfile_tsv} \\
  --outfile_fasta {outfile_fasta}""".format(
    min_reads = config.param('clustering', 'min_reads', 1, 'int'),
    dada2_tsv = dada2_tsv,
    outfile_tsv = outfile_tsv,
    outfile_fasta = outfile_fasta
 )
    
    return job

def convert_biom_to_tsv(infile_biom, outfile_tsv):
       
    job = Job(
        [infile_biom],
        [outfile_tsv],
        [
            ['qiime', 'module_qiime'],
            ['qiime-dependencies', 'module_qiime-dependencies'],
            ['tools', 'module_tools'],
            ['perl', 'module_perl'],
            ['lapack', 'module_lapack']
        ]
    )

#deblur workflow --seqs-fp ./infiles/ --output-dir ./deblur_outdir/  --overwrite -t 100
    job.command="""
biom convert \\
  -i {infile_biom} \\
  -o {outfile_tsv} \\
  --to-tsv""".format(
    infile_biom = infile_biom,
    outfile_tsv = outfile_tsv
 )
    
    return job

def remove_chimeras(infile_tsv, infile_fasta, outdir, outfile_tsv, outfile_fasta, method="denovo_and_ref"):
       
    job = Job(
        [infile_fasta, infile_tsv],
        [outfile_fasta, outfile_tsv],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl'],
            ['vsearch', 'module_vsearch']
        ]
    )

    job.command="""
scanningForChimeras.pl \\
  --infile_tsv {infile_tsv} \\
  --infile_fasta {infile_fasta} \\
  --outdir {outdir} \\
  --ref_db {ref_db} \\
  --num_threads {num_threads}""".format(
    infile_fasta = infile_fasta,
    infile_tsv = infile_tsv,
    ref_db = config.param('DB', 'chimeras', 1, 'filepath'),
    num_threads =  config.param('clustering', 'num_threads', 'int'),
    outfile_tsv = outfile_tsv,
    outfile_fasta = outfile_fasta,
    outdir = outdir
 )
    if(method == "ref_only"):
        job.command += " --ref_only"
    
    return job

def cleanup(tmpdir):
        
    job = Job(
        [""], 
        [""],
        [
            ['tools', 'module_tools']
        ]
    )
    
    job.command="""
rm  {tmpdir} -rf""".format(
    tmpdir = tmpdir
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


