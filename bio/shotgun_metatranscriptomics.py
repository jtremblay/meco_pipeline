#!/usr/bin/env python

#LICENSE AND COPYRIGHT

#Copyright (C) 2023 Julien Tremblay

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

def htseq_count(input, gff, output):
    stranded="no"
    input_type="sam"
    options = config.param("htseq", "options", 1, "string")

    job = Job(
        [input],
        [output],
        [
            ['qiime-dependencies', 'module_qiime-dependencies'],
            ['python', 'module_python'],
            ['htseq', 'module_htseq'],
            ['samtools', 'module_samtools'],
            ['lapack', 'module_lapack'],
            ['memtime', 'module_memtime']
        ]
    )
    
    job.command="""\
samtools view {input} | htseq-count {options} - {gff} \\
  --stranded={stranded} \\
  --format={input_type} > {output}""".format(
    options = options,
    stranded = stranded,
    input = input,
    gff = gff,
    output = output,
    input_type = input_type
    )
    
    return job

def merge_bedtools_counts(infiles, outfile_raw, outfile_normalized):
    job = Job(
        infiles,
        [outfile_raw, outfile_normalized],
        [
            ['memtime', 'module_memtime'],
            ['tools', 'module_meco_tools'],
            ['R', 'module_R']
        ]
    )
    job.command="""
memtime mergeGeneAbundance.R \\
  -i {infile} \\
  -o {outfile_raw} \\
  -n {outfile_normalized} && \\
  sed -i -E 's/.genes//g' {outfile_raw} && \\
  sed -i -E 's/.genes//g' {outfile_normalized}""".format(
    infile = ",".join(infiles),
    outfile_raw = outfile_raw,
    outfile_normalized = outfile_normalized
    )

    return job

def merge_htseq_counts(infiles_htseq, outfile_raw, outfile_normalized):
    job = Job(
        infiles_htseq,
        [outfile_raw, outfile_normalized],
        [
            ['memtime', 'module_memtime'],
            ['tools', 'module_meco_tools'],
            ['R', 'module_R']
        ]
    )

    job.command="""
memtime mergeHtseqAbundance.R \\
  -i {infiles_htseq} \\
  -o {outfile_raw} \\
  -n {outfile_normalized}""".format(
    infiles_htseq = ",".join(infiles_htseq),
    outfile_raw = outfile_raw,
    outfile_normalized = outfile_normalized
    )

    return job

def edger(abundance, design_file, outdir):
    
    basename =  os.path.splitext(os.path.basename(design_file))[0]

    dummy_outfile = os.path.join(outdir, "DEG_" + basename + ".done")

    job = Job(
        [abundance, design_file], 
        [dummy_outfile],
        [
            ['memtime', 'module_memtime'],
            ['meco_tools', 'module_meco_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
memtime edgerRNA.R \\
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
    logfc = config.param('DEG', 'logfc', 1, 'float'),
    fdr = config.param('DEG', 'fdr', 1, 'float'),
    pvalue = config.param('DEG', 'pvalue', 1, 'float')
    )

    return job

def edger_glm(abundance, mapping_file, outdir):

    basename = os.path.splitext(os.path.basename(mapping_file))[0]
    dummy_outfile = os.path.join(outdir, basename + "_DEG_GLM.done")

    job = Job(
        [abundance, mapping_file], 
        [dummy_outfile],
        [
            ['meco_tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
edgerGLM.R \\
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
    logfc = config.param('DEG', 'logfc', 1, 'float'),
    fdr = config.param('DEG', 'fdr', 1, 'float'),
    pvalue = config.param('DEG', 'pvalue', 1, 'float'),
    blocks = config.param('DEG', 'blocks', 1, 'string'),
    treatments = config.param('DEG', 'treatments', 1, 'string'),
    mapping_file = mapping_file
    )

    return job

def edger_single(infile, outfile):

    job = Job(
        [infile], 
        [outfile],
        [
            ['memtime', 'module_memtime'],
            ['meco_tools', 'module_meco_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
memtime edgerSingleRNA.R \\
  -i {infile} \\
  -o {outfile}""".format(
    infile = infile,
    outfile = outfile
    )

    return job

def rmarkdown_duk(log1, log2):

    job = Job(
        log1 + log2, 
        [outfile],
        [
            ['memtime', 'module_memtime'],
            ['meco_tools', 'module_meco_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
memtime edgerSingleRNA.R \\
  -i {infile} \\
  -o {outfile}""".format(
    infile = infile,
    outfile = outfile
    )

    return job

def explore_heatmaps(log_fc, normalized_significant, indir, outdir, date, gene_list):
   
    dummy_outfile = os.path.join(outdir, "explore_heatmaps.done")

    job = Job(
        [normalized_significant], 
        [dummy_outfile],
        [
            ['memtime', 'module_memtime'],
            ['meco_tools', 'module_meco_tools'],
            ['R', 'module_R'],
            ['pandoc', 'module_pandoc']
        ]
    )
    
    job.command="""
memtime R -e "thedate='{date}'; rmarkdown::render('$RMD_PATH/metatranscriptomics_exploratory_heatmaps.Rmd', output_dir=file.path('{outdir}'))" --args {log_fc} {normalized_significant} {annotations} {mapping_file} {indir} {outdir} {gene_list}""".format(
    date = date,
    log_fc = log_fc,
    normalized_significant = os.path.abspath(normalized_significant),
    annotations = os.path.abspath(config.param('DEFAULT', 'annotations', 1, 'filepath')),
    mapping_file = os.path.abspath(config.param('DEFAULT', 'mappingFile', 1, 'filepath')),
    indir = os.path.abspath(indir),
    outdir = os.path.abspath(outdir),
    gene_list = os.path.abspath(gene_list)
    )

    return job

# Methods for rmarkdowns of each steps.


