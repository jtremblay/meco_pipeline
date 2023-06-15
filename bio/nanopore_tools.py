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

# MUGQIC Modules
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
    headcrop = config.param('trim', 'headcrop', required=False, type='posint')
    average_quality = config.param('trim', 'average_quality', required=False, type='posint')

    job.command = """
java -XX:ParallelGCThreads={threads} -Xmx2G -jar \$TRIMMOMATIC_JAR {mode} \\
  -threads {threads} \\
  -phred{quality_offset} \\
  {input1} {output1} \\
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
        min_length = min_length,
        headcrop = str(headcrop),
        average_quality = average_quality
    )
    #ILLUMINACLIP:{adapter_file}{illumina_clip_settings} \\

    job.command += " \\\n  2> " + trim_log

    # Compute statistics
    job.command += " && \\\ngrep ^Input " + trim_log + " | perl -pe 's/^Input Read Pairs: (\\d+).*Both Surviving: (\\d+).*Forward Only Surviving: (\\d+).*$/Raw Fragments,\\1\\nFragment Surviving,\\2\\nSingle Surviving,\\3/' > " + trim_stats

    return job

#def filtlong(infile, outfile):
#        
#    job = Job(
#        [infile],
#        [outfile],
#        [
#            ['filtlong', 'module_filtlong']
#        ]
#    )
#    job.command="""
#    filtlong \\
#       --min_length 1000 
#       --keep_percent 90 
#       --target_bases 500000000 
#       input.fastq.gz | gzip > output.fastq.gz""".format(
#        min_length = config.param('filtlong', 'min_length', 'posint'),
#        keep_percent = config.param('filtlong', 'keep_percent', 'posint'),
#        target_bases = config.param('filtlong', 'target_bases', 'posint')
#    )
#
#    return job


#getNumberOfBasesInFasta.pl --infile {output_dir}/data/filtered_subreads.fasta > {output_dir}/data/number_of_bp.txt""".format(


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


def cat_fastqs_and_convert_to_fasta(infiles, outfile, outfile_bp):
     
    job = Job(
        infiles,
        [outfile, outfile_bp],
        [
            ['perl', 'module_perl'],
            ['caf_tools', 'module_tools']
        ]
    )

    job.command="""
zcat {infiles} | sed -n '1~4s/^@/>/p;2~4p' > {outfile} && \\
getNumberOfBasesInFasta.pl --infile {outfile} > {outfile_bp}""".format(
    infiles = " ".join(infiles),
    outfile = outfile,
    outfile_bp = outfile_bp
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

def miniasm(infile_fastq, infile_alignment, outfile):
     
    job = Job(
        [infile_fastq, infile_alignment],
        [outfile],
        [
            ['miniasm', 'module_miniasm']
        ]
    )

#miniasm -f qced_reads/all_qced_reads.fastq qced_reads/all_qced_reads.paf.gz > qced_reads/reads.gfa
#awk '/^S/{print ">"$2"\n"$3}' qced_reads/reads.gfa > assembly/Contigs.fasta

    job.command="""
miniasm -f {type} \\
  {infile_fasta} \\
  {infile_alignment} \\
  > {outfile}
""".format(
    infile_alignment,
    infile_fastq,
    outfile = outfile
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
compileRayResultsSingle.pl \\
  --infile {infile} > {outfile}""".format(
    infile = infile,
    outfile = outfile
     )

    return job

