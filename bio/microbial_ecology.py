#!/usr/bin/env python

# Python Standard Modules

# MECO Modules
from core.config import *
from core.job import *

def rdp_wrapper(infile, outfile, forced_exit_status=1):
    job = Job(
        [infile],
        [outfile],
        [
            ['tools', 'module_tools'],
            ['perl', 'module_perl'],
            ['rdp_classifier', 'module_rdp_classifier'],
            ['java', 'module_java']
        ]
    )

    job.command="""

if [[ -s {infile} ]] ; then
    parallelRDP.pl \\
     --infile {infile} \\
     --RDP_training_set {rdp_training_set} \\
     --outfile {outfile} \\
     --minWords {min_words} \\
     --num_threads {num_threads} \\
     --ram {ram} \\
     --fasta &&
    if [[ -s {outfile} ]] ; then
        echo "{outfile} was sucessfully generated: exists and not empty."
    else
        echo "{outfile} was not sucessfully generated: exists, but empty."
        echo "Most probable cause: not enough RAM was requested for the java heap."
        echo "Increase the RAM=\<int\>g in the \[RDP\] parameters in the .ini file and try again".
        exit {forced_exit_status}
    fi

else
    touch {outfile}
    echo "{infile} is empty! ... Something went wrong in an upstream step."
    echo "Verify that proper parameters were specified for the OTU/ASV tables generation."
    exit {forced_exit_status}
fi""".format(
    infile = infile,
    rdp_training_set = config.param('DB', 'rdp_training_set', 1, 'filepath'),
    outfile = outfile,
    min_words = config.param('RDP', 'minWords', 1, 'int'),
    num_threads = config.param('RDP', 'num_threads', 1, 'int'),
    ram = config.param('RDP', 'RAM', 1, 'string'),
    forced_exit_status=forced_exit_status
    )

    return job

def generate_feature_table(obs, rdp, outfile, outfile_failed):
    job = Job(
        [obs, rdp],
        [outfile],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
addTaxToObs.pl \\
  --seqobs {obs} \\
  --rdp {rdp} \\
  --cutoff {cutoff} \\
  --outfile {outfile} \\
  --outfile_failed {outfile_failed} \\
  --tax_level {tax_level}""".format(
        obs = obs,
        rdp = rdp,
        cutoff = config.param('add_taxonomy', 'cutoff', 1, 'float'),
        outfile = outfile,
        outfile_failed = outfile_failed,
        tax_level = config.param('add_taxonomy', 'tax_level', 1, 'string')
    )
    
    return job

def generate_feature_table_closed(obs, blast_output, outfile, outfile_failed, link):
    job = Job(
        [obs, blast_output],
        [outfile, link],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
addTaxToObsClosed.pl \\
  --seqobs {obs} \\
  --blast_output {blast_output} \\
  --evalue {evalue} \\
  --al_length {al_length} \\
  --perc_id {perc_id} \\
  --outfile {outfile} \\
  --outfile_failed {outfile_failed} \\
  --link {link}""".format(
        obs = obs,
        blast_output = blast_output,
        perc_id = config.param('add_taxonomy', 'perc_id_blast', 1, 'float'),
        evalue = config.param('add_taxonomy', 'evalue_blast', 1, 'float'),
        al_length = config.param('add_taxonomy', 'al_length_blast', 1, 'int'),
        outfile = outfile,
        outfile_failed = outfile_failed,
   	link = link 
   )
    
    return job

def collapse_duplicates(infile, outfile):
    job = Job(
        [infile],
        [outfile],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
collapseFeatureTableDuplicates.R \\
  -i {infile} \\
  -o {outfile}""".format(
        outfile = outfile,
   	infile = infile 
   )
    
    return job

def predict_traits(traits, tree, counts, feature_table, feature_table_predicted):
    job = Job(
        [traits, tree, counts, feature_table],
        [feature_table_predicted],
        [
            ['perl', 'module_perl'],
            ['qiime', 'module_qiime'],
            ['qiime-dependencies', 'module_qiime-dependencies'],
	    ['picrust', 'module_picrust']
        ]
    )

    job.command="""
predict_traits.py \\
  -i {traits} \\
  -t {tree} \\
  -r {counts} \\
  -l {feature_table} \\
  -o {feature_table_predicted}""".format(
        traits = traits,
        tree = tree,
        counts = counts,
        feature_table = feature_table,
        feature_table_predicted = feature_table_predicted
    )
    
    return job

def predict_metagenomes(infile, outfile_biom, outfile_tsv):
    job = Job(
        [infile],
        [outfile_biom, outfile_tsv],
        [
            ['perl', 'module_perl'],
            ['qiime', 'module_qiime'],
            ['qiime-dependencies', 'module_qiime-dependencies'],
	    ['picrust', 'module_picrust']
        ]
    )

    job.command="""
predict_metagenomes.py \\
  -i {infile} \\
  -o {outfile_biom} && rm -f {outfile_tsv} && \\
biom convert \\
  -i {outfile_biom} \\
  -o {outfile_tsv} \\
  --to-tsv """.format(
        infile = infile,
        outfile_biom = outfile_biom,
        outfile_tsv = outfile_tsv
    )
    
    return job

def normalize_by_copy_number(infile, outfile_biom, outfile_tsv):
    job = Job(
        [infile],
        [outfile_biom, outfile_tsv],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools'],
            ['qiime', 'module_qiime'],
            ['qiime-dependencies', 'module_qiime-dependencies'],
	    ['picrust', 'module_picrust']
        ]
    )

    job.command="""
normalize_by_copy_number.py \\
  -i {infile} \\
  -o {outfile_biom} && rm -f {outfile_tsv} && \\
biom convert \\
  -i {outfile_biom} \\
  -o {outfile_tsv} \\
  --table-type 'OTU table' --to-tsv --output-metadata-id='taxonomy' --header-key='taxonomy'""".format(
        outfile_biom = outfile_biom,
        outfile_tsv = outfile_tsv,
   	infile = infile 
   )
    
    return job

def filter_samples(infile, outfile_biom, outfile_tsv):
    job = Job(
        [infile],
        [outfile_biom, outfile_tsv],
        [
            ['perl', 'module_perl'],
            ['qiime', 'module_qiime'],
            ['qiime-dependencies', 'module_qiime-dependencies'],
            ['python', 'module_python']
        ]
    )

    job.command="""
filter_samples_from_otu_table.py \\
  -i {infile} \\
  -o {outfile_biom} \\
  -n {min_count} && \\
rm -f {outfile_tsv} && \\
   biom convert \\
    -i {outfile_biom} \\
    -o {outfile_tsv} \\
  --table-type 'OTU table' --to-tsv --output-metadata-id='taxonomy' --header-key='taxonomy'""".format(
        infile = infile,
        outfile_biom = outfile_biom,
        outfile_tsv = outfile_tsv,
        min_count = config.param('multiple_rarefaction', 'rarefaction_threshold', 1, 'string')
    )

    return job

# Deprecetated?...
def filter_feature_table(infile, outfile):
    job = Job(
        [infile],
        [outfile],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
rmEmptyCol.pl \\
  --infile {infile} \\
  --outfile {outfile}""".format(
        infile = infile,
        outfile = outfile
    )

    return job

def feature_table_multiplier(infile, outfile, multiplier):
    job = Job(
        [infile],
        [outfile],
        [
            ['R', 'module_R'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
featureTableMultiplier.R \\
  -i {infile} \\
  -o {outfile} \\
  -n {multiplier}""".format(
        infile = infile,
        outfile = outfile,
        multiplier = multiplier
    )

    return job

def filter_and_sort_feature_table(infile, outfile):
    job = Job(
        [infile],
        [outfile],
        [
            ['perl', 'module_perl'],
            ['R', 'module_R'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
rmEmptyCol.pl \\
  --infile {infile} \\
  --outfile {outfile_tmp} && \\
sortFeatureTable.R -i {outfile_tmp} -o {outfile}""".format(
        infile = infile,
        outfile_tmp = str(outfile + ".tmp"),
        outfile = outfile
    )

    return job

def split_feature_table(infile, matched, unmatched, selection):
    job = Job(
        [infile],
        [matched, unmatched],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
splitFeatureTable.pl \\
  --infile {infile} \\
  --matched {matched} \\
  --unmatched {unmatched} \\
  --keepChloro {chloro} \\
  --keepMito {mito} \\
  --select {selection}""".format(
        infile = infile,
        matched = matched,
        unmatched = unmatched,
        selection = selection,
        chloro = config.param('split_feature_table', 'keep_chloroplast', 1, 'string'),
        mito = config.param('split_feature_table', 'keep_mitochondria', 1, 'string')
    )

    return job

def split_feature_table_multi(infile, outdir, unmatched, selection):
    
    outfiles = []
    basename = os.path.splitext(infile)[0]
    for sel in selection:
        outfiles.append(os.path.join(outdir, basename + "_" + sel + ".tsv"))
    
    job = Job(
        [infile],
        [unmatched] + outfiles,
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
splitFeatureTableMulti.pl \\
  --infile {infile} \\
  --outdir {outdir} \\
  --unmatched {unmatched} \\
  --select {selection}""".format(
        infile = infile,
        unmatched = unmatched,
        selection = selection,
        outdir = outdir
    )

    return job
# This is good for converting OTU table in tsv to biom hdf5 format (>biom v2.0)
def convert_tsv_to_biom_hdf5(tsv, biom):
    job = Job(
        [tsv],
        [biom],
        [
            ['qiime', 'module_qiime'],
            ['qiime-dependencies', 'module_qiime-dependencies'],
            ['python', 'module_python'],
            ['lapack', 'module_lapack']
        ]
    )

    job.command="""
rm -f {biom} && \\
biom convert \\
  -i {tsv} \\
  -o {biom} \\
  --table-type 'OTU table' --to-hdf5 --process-obs-metadata taxonomy""".format(
  biom = biom,
  tsv = tsv
    )
    
    return job


def convert_biom_to_tsv_multi(indir, dummy_infile, dummy_outfile):
    job = Job(
        [dummy_infile],
        [dummy_outfile],
        [
            ['perl', 'module_perl'],
            ['qiime-dependencies', 'module_qiime-dependencies'],
            ['qiime', 'module_qiime'],
            ['python', 'module_python'],
            ['tools', 'module_tools'],
            ['lapack', 'module_lapack']
        ]
    )

    job.command="""
    convertBiomToTsvMulti.pl --indir {indir} --num_threads {num_threads} && \\
       touch {dummy_outfile}""".format(
        indir = indir,
        num_threads = config.param('convert_feature_tables', 'num_threads', 1, 'int'),
        dummy_outfile = dummy_outfile
    )

    return job

def merge_rarefied_tables(indir, infile_tax, dummy_infile, outfile_mean):
    job = Job(
        [dummy_infile, infile_tax],
        [outfile_mean],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
mergeRarefiedTables.pl \\
  --indir {indir} \\
  --infile {infile_tax} \\
  > {outfile_mean}""".format(
        indir = indir,
        infile_tax = infile_tax,
        outfile_mean = outfile_mean
    )

    return job

def normalize_feature_table(feature_table, tsv_normalized):
    job = Job(
        [feature_table],
        [tsv_normalized],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
#-m {m} && \\
    job.command="""
normalizeFeatureTable.R \\
  -i {feature_table} \\
  -o {tsv_normalized} \\
  -c {n}""".format(
        feature_table = feature_table,
        tsv_normalized = tsv_normalized,
        n = config.param('normalization', 'normalization_threshold', 1, 'int')
    )

    return job

def percentage_feature_table(feature_table, tsv_percentaged):
    job = Job(
        [feature_table],
        [tsv_percentaged],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
#-m {m} && \\
    job.command="""
percentageFeatureTable.R \\
  -i {feature_table} \\
  -o {tsv_percentaged}""".format(
        feature_table = feature_table,
        tsv_percentaged = tsv_percentaged
    )

    return job

#def summarize_taxonomy_absolute(biom, i, outdir):
#    basename =  os.path.splitext(os.path.basename(biom))[0]
#    
#    outbiom = os.path.join(outdir, basename + "_L" + str(i) + ".txt")
#    outtsv = os.path.join(outdir, basename + "_L" + str(i) + ".biom")
#
#    job = Job(
#        [biom],
#        [outbiom, outtsv],
#        [
#            ['qiime-dependencies', 'module_qiime-dependencies'],
#            ['qiime', 'module_qiime'],
#            ['python', 'module_python'],
#            ['lapack', 'module_lapack']
#        ]
#    )
#
#    job.command="""
#summarize_taxa.py \\
#  -i {biom} \\
#  -L {i} \\
#  -o {outdir} \\
#  -a && \\
#  sed -i 's/#OTU ID/Taxon/' {outdir}/{basename}_L{i}.txt && \\
#  perl -ni -e 'print unless(\$_ =~ m/Constructed from biom file/)' {outdir}/{basename}_L{i}.txt""".format(
#        biom = biom,
#        i = i,
#        outdir = outdir,
#        basename = basename
#    )
#
#    return job
#
#def summarize_taxonomy_relative(biom, i, outdir):
#    basename =  os.path.splitext(os.path.basename(biom))[0]
#
#    outbiom = os.path.join(outdir, basename + "_L" + str(i) + ".txt")
#    outtsv = os.path.join(outdir, basename + "_L" + str(i) + ".biom")
#
#    job = Job(
#        [biom],
#        [outbiom, outtsv],
#        [
#            ['qiime-dependencies', 'module_qiime-dependencies'],
#            ['qiime', 'module_qiime'],
#            ['python3', 'module_python3'],
#            ['lapack', 'module_lapack']
#        ]
#    )
#
#    job.command="""
#summarize_taxa.py \\
#  -i {biom} \\
#  -L {i} \\
#  -o {outdir} && \\
#  sed -i 's/#OTU ID/Taxon/' {outdir}/{basename}_L{i}.txt && \\
#  perl -ni -e 'print unless(\$_ =~ m/Constructed from biom file/)' {outdir}/{basename}_L{i}.txt""".format(
#        biom = biom,
#        i = i,
#        outdir = outdir,
#        basename = basename
#    )
#
#    return job

def summarize_taxonomy_absolute(tsv, i, outdir):
    basename =  os.path.splitext(os.path.basename(tsv))[0]
    
    outtsv = os.path.join(outdir, basename + "_L" + str(i) + ".tsv")

    job = Job(
        [tsv],
        [outtsv],
        [
            ['python3', 'module_python3'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
microbiomeutils.py taxsum \\
  -i {tsv} \\
  -l {i} \\
  -t absolute \\
  > {outtsv}""".format(
        tsv = tsv,
        i = i,
        outtsv = outtsv
    )

    return job

def summarize_taxonomy_relative(tsv, i, outdir):
    basename =  os.path.splitext(os.path.basename(tsv))[0]

    outtsv = os.path.join(outdir, basename + "_L" + str(i) + ".tsv")

    job = Job(
        [tsv],
        [outtsv],
        [
            ['python3', 'module_python3'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
microbiomeutils.py taxsum \\
  -i {tsv} \\
  -l {i} \\
  -t relative \\
  > {outtsv}""".format(
        tsv = tsv,
        i = i,
        outtsv = outtsv
  )

    return job

def plot_taxa_single(infile, outdir, prefix):
    job = Job(
        [infile],
        [outdir + "/" + prefix + ".pdf"],
        [
            ['qiime', 'module_R'],
            ['python', 'module_tools']
        ]
    )

    job.command="""
taxBarplot.R \\
  -i {infile} \\
  -o {outdir} \\
  -p {prefix}""".format(
    infile = infile,
    outdir = outdir + "/",
    prefix = prefix
    )

    return job

def plot_taxa_single_with_mapping_file(infile, outdir, prefix):

    dummy_outfile = os.path.join(outdir, prefix + ".done")

    job = Job(
        [infile],
        [dummy_outfile],
        [
            ['qiime', 'module_R'],
            ['python', 'module_tools']
        ]
    )

    job.command="""
taxBarplotWithMappingFile.R \\
  -i {infile} \\
  -o {outdir} \\
  -p {prefix} \\
  -m {mapping_file} && touch {dummy_outfile}""".format(
    infile = infile,
    outdir = outdir + "/",
    prefix = prefix,
    mapping_file = config.param('default', 'mapping_file', 1, 'filepath'),
    dummy_outfile = dummy_outfile
    )

    return job

def plot_pca(infile, outdir, prefix):
    dummy_outfile = outdir + "/pcaPlots.done"
    job = Job(
        [infile],
        [dummy_outfile],
        [
            ['qiime', 'module_R'],
            ['python', 'module_tools']
        ]
    )

    job.command="""
plotPCA.R \\
  -i {infile} \\
  -o {outdir} \\
  -p {prefix} \\
  -m {mapping_file} &&
  touch {dummy_outfile}""".format(
    infile = infile,
    outdir = outdir + "/",
    prefix = prefix,
    mapping_file = config.param('default', 'mapping_file', 1, 'filepath'),
    dummy_outfile = dummy_outfile
    )

    return job


def filter_for_mafft(tsv, fasta, outfile):
    job = Job(
        [tsv, fasta],
        [outfile],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
filterForAlignment.pl \\
  --infile_feature_table {tsv} \\
  --infile_fasta {fasta} \\
  > {outfile}""".format(
    tsv = tsv,
    fasta = fasta,
    outfile = outfile
    )

    return job

def pynast(infile, log, outfile):
    job = Job(
        [infile],
        [outfile, log],
        [
            ['qiime-dependencies', 'module_qiime-dependencies'],
            ['qiime', 'module_qiime'],
            ['python', 'module_python'],
            ['pynast', 'module_pynast'],
            ['openmpi', 'module_openmpi'],
            ['lapack', 'module_lapack']
        ]
    )

#mpirun -np {num_threads} \\
    job.command="""
  pynast \\
  -i {infile} \\
  -p 10 \\
  -l 50 \\
  -g {log} \\
  -a {outfile} \\
  -t {db}""".format(
    num_threads = config.param('pynast', 'num_threads', 1, 'int'),
    infile = infile,
    log = log,
    outfile = outfile,
    db = config.param('DB', 'core', 1, 'filepath')
    )

    return job

def mafft(infile, outfile):
    job = Job(
        [infile],
        [outfile],
        [
            ['mafft', 'module_mafft']
        ]
    )

    job.command="""
mafft \\
  --progress /dev/null \\
  {infile} \\
  > {outfile}""".format(
    infile = infile,
    outfile = outfile
    )

    return job

def filter_pynast_alignment(infile, outdir, suppress_lane=False):
    job = Job(
        [infile],
        [outdir + "/pynast_pfiltered.fasta"],
        [
            ['qiime-dependencies', 'module_qiime-dependencies'],
            ['qiime', 'module_qiime'],
            ['python', 'module_python'],
            ['lapack', 'module_lapack']
        ]
    )
    if(suppress_lane == False):
        job.command="""
    filter_alignment.py \\
    -i {infile} \\
    -o {outdir}""".format(
        infile = infile,
        outdir = outdir
        )
    else:
        job.command="""
    filter_alignment.py \\
    -i {infile} \\
    -o {outdir} -s""".format(
        infile = infile,
        outdir = outdir
        )


    return job

def fasttree(infile, outfile):
    job = Job(
        [infile],
        [outfile],
        [
            ['fasttree', 'module_fasttree']
        ]
    )

    job.command="""
FastTree \\
  -nt {infile} \\
  > {outfile}""".format(
    infile = infile,
    outfile = outfile
    )

    return job

#def beta_diversity(biom, distance, outdir, tree, prefix):
#    
#    outfile = os.path.join(outdir, prefix + ".txt")
#    sys.stderr.write("betadiv outfile: " + outfile + "\n")
#
#    job = Job(
#        [biom, tree],
#        [outfile],
#        [
#            ['qiime-dependencies', 'module_qiime-dependencies'],
#            ['qiime', 'module_qiime'],
#            ['python', 'module_python'],
#            ['lapack', 'module_lapack']
#        ]
#    )
##unset LD_LIBRARY_PATH;
#
#    job.command="""
#export OMPI_MCA_mtl=^psm && \\
#beta_diversity.py \\
#  -i {biom} \\
#  -m {distance} \\
#  -o {outdir} \\
#  -t {tree} \\
#  && unset OMPI_MCA_mtl \\
#  && sed -i 's/nan/0.0/g' {outfile}""".format( # Usually this is necessary if one sample have very few reads.
#    biom = biom,
#    distance = distance,
#    outdir = outdir,
#    tree = tree,
#    outfile = outfile
#    )
#
#    return job

def beta_diversity(tsv, distance, outdir, tree, prefix):
    
    outfile = os.path.join(outdir, prefix + ".tsv")
    #sys.stderr.write("betadiv outfile: " + outfile + "\n")

    job = Job(
        [tsv, tree],
        [outfile],
        [
            ['python3', 'module_python3'],
            ['tools', 'module_tools']
        ]
    )
#unset LD_LIBRARY_PATH;

    job.command="""
microbiomeutils.py betadiv \\
  -i {tsv} \\
  -m {distance} \\
  -t {tree} \\
  > {outfile} \\
  && sed -i 's/nan/0.0/g' {outfile}""".format( # Usually this is necessary if one sample have very few reads.
    tsv = tsv,
    distance = distance,
    outdir = outdir,
    tree = tree,
    outfile = outfile
    )

    return job

def beta_diversity_R(biom, metric, outdir, tree, prefix):
    job = Job(
        [biom, tree],
        [os.path.join(outdir, prefix + ".txt"), os.path.join(outdir, prefix + "_coords.txt")],
        [
            ['R', 'module_R'],
            ['tools', 'module_tools']
        ]
    )
    job.command="""
computeDistances.R \\
  -i {biom} \\
  -m {metric} \\
  -o {outdir} \\
  -t {tree} \\
  -p {prefix}""".format(
    biom = biom,
    metric = metric,
    outdir = outdir,
    tree = tree,
    prefix = prefix
    )

    return job

#def principal_coordinates(infile, outfile):
#    job = Job(
#        [infile],
#        [outfile],
#        [
#            ['qiime-dependencies', 'module_qiime-dependencies'],
#            ['qiime', 'module_qiime'],
#            ['python', 'module_python'],
#            ['lapack', 'module_lapack']
#        ]
#    )
#
#    job.command="""
#principal_coordinates.py \\
#  -i {infile} \\
#  -o {outfile}""".format(
#    infile = infile,
#    outfile = outfile
#    )
#
#    return job

def principal_coordinates(infile, outfile):
    job = Job(
        [infile],
        [outfile],
        [
            ['python3', 'module_python3'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
microbiomeutils.py pcoa \\
  -i {infile} \\
  > {outfile}""".format(
    infile = infile,
    outfile = outfile
    )

    return job

#def pca_2d_plot(infile, outdir):
#    base = os.path.basename(infile)
#    basename = os.path.splitext(base)[0]
#    outfile = os.path.join(outdir, basename + "_2D_PCoA_plots.html")
#
#    #sys.stderr.write("basename : "  + basename + "\n")
#    #sys.stderr.write("outfile : "  + outfile + "\n")
#    #sys.stderr.write("2d plot outfile : "  + outdir + "" + basename + "_2D_PCoA_plots.html" + "\n")
#    job = Job(
#        [infile],
#        [outfile],
#        [
#            ['qiime-dependencies', 'module_qiime-dependencies'],
#            ['qiime', 'module_qiime'],
#            ['python', 'module_python'],
#            ['lapack', 'module_lapack']
#        ]
#    )
#
#    job.command="""
#make_2d_plots.py \\
#  -i {infile} \\
#  -m {mapping_file} \\
#  -o {outdir}""".format(
#    infile = infile,
#    mapping_file = config.param('default', 'mapping_file', 1, 'filepath'),
#    outdir = outdir
#    )
#
#    return job

def pca_3d_plot(infile, outdir):
    base = os.path.basename(infile)
    basename = os.path.splitext(base)[0]
    #outfile = outdir + "" + basename + "_3D_PCoA_plots.html"
    outfile = os.path.join(outdir, "index.html")
    #regex = "'s/@group \({.*(.*)}\) \(.*\)$/@group \\1 \\2\\n@group \\1 animate/'"

    job = Job(
        [infile],
        [outfile],
        [
            ['python3', 'module_python3'],
            ['tools', 'module_tools']
        ]
    )
    job.command="""
microbiomeutils.py emperor \\
  -i {infile} \\
  -m {mapping_file} \\
  -o {outdir}""".format(
    infile = infile,
    mapping_file = config.param('default', 'mapping_file', 1, 'filepath'),
    outdir = outdir,
    outfile = outfile
    )

    return job

# Rarefaction for OTU table normalization.
def rarefy_hdf5_multi_value(biom, tsv, outdir, root_outdir, n_value, dummy_outfile):

    job = Job(
        [biom],
        [dummy_outfile],
        [
            ['perl', 'module_perl'],
            ['qiime-dependencies', 'module_qiime-dependencies'],
            ['qiime', 'module_qiime'],
            ['python', 'module_python'],
            ['tools', 'module_tools'],
            ['lapack', 'module_lapack']
        ]
    )

    job.command="""
rm -rf {root_outdir}/rarefaction/* && \\
multipleRarefactionsEvenDepth.pl \\
  --infile {biom} \\
  --infile_tsv {tsv} \\
  --outdir {outdir} \\
  --n {rarefaction_threshold} \\
  --step {step} \\
  --permutations {perm} \\
  --num_threads {num_threads} \\
  && touch {dummy_outfile}""".format(
    root_outdir = root_outdir,
    biom = biom,
    tsv = tsv,
    outdir = outdir,
    #min_fraction_threshold = config.param('rarefaction', 'minFractionThreshold', 1, 'float'),
    rarefaction_threshold = n_value,
    step = 1,
    perm = config.param('rarefaction', 'perm_normalization', 1, 'int'),
    num_threads = config.param('rarefaction', 'num_threads', 1, 'int'),
    dummy_outfile = dummy_outfile
    )

    return job

# Rarefy OTU table for alpha div and/or alpha div saturation
def multiple_rarefaction(biom, tsv, outdir, root_outdir, num_reads, step, perm, dummy_outfile):

    job = Job(
        [biom],
        [dummy_outfile],
        [
            ['perl', 'module_perl'],
            ['qiime-dependencies', 'module_qiime-dependencies'],
            ['qiime', 'module_qiime'],
            ['python', 'module_python'],
            ['tools', 'module_tools'],
            ['lapack', 'module_lapack']
        ]
    )
    
    if num_reads == "auto" and step == "auto": #for saturation
        job.command="""
rm -rf {outdir}/* && \\
rm -rf {root_outdir}/alpha_rarefaction/* && \\
rm -rf {root_outdir}/collated/* && \\
rm -rf {root_outdir}/plots/* && \\
multipleRarefactions.pl \\
  --infile {biom} \\
  --infileTsv {tsv} \\
  --outdir {outdir} \\
  --permutations {perm} \\
  --num_threads {num_threads} && \\
  touch {dummy_outfile} && \\
  rm -rf {outdir}/RARIF_*""".format(
    root_outdir = root_outdir,
    biom = biom,
    tsv = tsv,
    outdir = outdir,
    perm = perm,
    num_threads = config.param('multiple_rarefaction', 'num_threads', 1, 'int'),
    dummy_outfile = dummy_outfile
    )

    else: # for fixed rarefaction
        job.command="""
rm -rf {root_outdir}/rarefaction/* && \\
rm -rf {root_outdir}/alpha_rarefaction/* && \\
rm -rf {root_outdir}/collated/* && \\
rm -rf {root_outdir}/plots/* && \\
multipleRarefactions.pl \\
  --infile {biom} \\
  --infileTsv {tsv} \\
  --outdir {outdir} \\
  --n {rarefaction_threshold} \\
  --step {step} \\
  --permutations {perm} \\
  --num_threads {num_threads} && \\
  touch {dummy_outfile} && \\
  rm -rf {outdir}/RARIF_*""".format(
    root_outdir = root_outdir,
    biom = biom,
    tsv = tsv,
    outdir = outdir,
    rarefaction_threshold = num_reads,
    step = step,
    perm = perm,
    num_threads = config.param('multiple_rarefaction', 'num_threads', 1, 'int'),
    dummy_outfile = dummy_outfile
    )
#--n {rarefaction_threshold} \\
    return job


#def alpha_diversity(indir, outdir, dummy_in_base="rarefactions.done"):
#    dummy_outdir = os.path.relpath(os.path.join(outdir, os.pardir))
#    dummy_outdir = "./" + dummy_outdir
#    dummy_outfile = dummy_outdir + "/alphaDiv.done"
#    dummy_infile = dummy_outdir + "/" + dummy_in_base
#
#    #sys.stderr.write("alphadiv infile: " + dummy_infile + "\n")
#    #sys.stderr.write("alphadiv outfile: " + dummy_outfile + "\n")
#    
#    m = config.param('alpha_diversity', 'm', 1, 'string')
#    m = m.replace(":",",")
#
#    job = Job(
#        [dummy_infile],
#        [dummy_outfile],
#        [
#            ['qiime-dependencies', 'module_qiime-dependencies'],
#            ['qiime', 'module_qiime'],
#            ['python', 'module_python'],
#            #['gcc', 'module_gcc'],
#            #['openmpi', 'module_openmpi'],
#            ['lapack', 'module_lapack']
#        ]
#    )
#
##unset LD_LIBRARY_PATH && \\
#    job.command="""
#rm -rf {outdir}/* && \\
#export OMPI_MCA_mtl=^psm && \\
#parallel_alpha_diversity.py \\
#  -i {indir} \\
#  -o {outdir} \\
#  -m {m} \\
#  -O {num_threads} && \\
#  touch {dummy_outfile} && \\
#  unset OMPI_MCA_mtl && \\
#  rm -rf {outdir}/ALDIV_*""".format(
#    indir = indir,
#    outdir = outdir,
#    m = m,
#    num_threads = config.param('alpha_diversity', 'num_threads', 1, 'int'),
#    dummy_outfile = dummy_outfile
#    )
#
#    return job
#
#def alpha_diversity_serial(indir, outdir):
#    dummy_outdir = os.path.relpath(os.path.join(outdir, os.pardir))
#    dummy_outdir = "./" + dummy_outdir
#    dummy_outfile = dummy_outdir + "/alphaDiv.done"
#    dummy_infile = dummy_outdir + "/rarefactions.done"
#
#    #sys.stderr.write("alphadiv infile: " + dummy_infile + "\n")
#    #sys.stderr.write("alphadiv outfile: " + dummy_outfile + "\n")
#    
#    m = config.param('alpha_diversity', 'm', 1, 'string')
#    m = m.replace(":",",")
#
#    job = Job(
#        [dummy_infile],
#        [dummy_outfile],
#        [
#            ['qiime-dependencies', 'module_qiime-dependencies'],
#            ['qiime', 'module_qiime'],
#            ['python', 'module_python'],
#            #['gcc', 'module_gcc'],
#            #['openmpi', 'module_openmpi'],
#            ['lapack', 'module_lapack']
#        ]
#    )
#
#    job.command="""
#rm -rf {outdir}* && \\
#export OMPI_MCA_mtl=^psm && \\
#alpha_diversity.py \\
#  -i {indir} \\
#  -o {outdir} \\
#  -m {m} \\
#  && touch {dummy_outfile} && \\
#  unset OMPI_MCA_mtl && \\
#  rm -rf {outdir}/ALDIV_*""".format(
#    indir = indir,
#    outdir = outdir,
#    m = m,
#    num_threads = config.param('alpha_diversity', 'num_threads', 1, 'int'),
#    dummy_outfile = dummy_outfile
#    )
#
#    return job
#
#def collate_alpha(indir, outdir):
#    dummy_outdir = os.path.relpath(os.path.join(outdir, os.pardir))
#    dummy_outdir = "./" + dummy_outdir
#    dummy_infile = dummy_outdir + "/alphaDiv.done"
#    dummy_outfile = dummy_outdir + "/collateAlpha.done"
#    #sys.stderr.write("collate infile: " + dummy_infile + "\n")
#
#    job = Job(
#        [dummy_infile],
#        [dummy_outfile],
#        [
#            ['qiime-dependencies', 'module_qiime-dependencies'],
#            ['qiime', 'module_qiime'],
#            ['python', 'module_python'],
#            ['lapack', 'module_lapack']
#        ]
#    )
#
#    job.command="""
#collate_alpha.py \\
#  -i {indir} \\
#  -o {outdir} \\
#  && touch {dummy_outfile}""".format(
#    indir = indir,
#    outdir = outdir,
#    dummy_outfile = dummy_outfile
#    )
#
#    return job

def rarefaction_plots(indir, outdir):
    dummy_outdir = os.path.relpath(os.path.join(outdir, os.pardir))
    dummy_outdir = "./" + dummy_outdir
    dummy_infile = indir + "/rtk.done"
    dummy_outfile = outdir + "/rarefactionPlots.done"
    observed_otus = dummy_outdir + "/alpha_richness.tsv"
    shannon = dummy_outdir + "/alpha_shannon.tsv"
    #simpson = dummy_outdir + "/collated/simpson.txt"
    #chao1 = dummy_outdir + "/collated/chao1.txt"
    
    job = Job(
        [dummy_infile],
        [dummy_outfile],
        [
            ['R', 'module_R'],
            ['tools', 'module_tools']
        ]
    )
    job.command = """
  plotAlphaDivWithFactors.R \\
   -i {observed_otus} \\
   -o {outdir} \\
   -p observedFeatures \\
   -m {mapping_file} \\
  && \\
  plotAlphaDivWithFactors.R \\
   -i {shannon} \\
   -o {outdir} \\
   -p shannon \\
   -m {mapping_file} \\
  && touch {dummy_outfile}""".format(
    indir = indir,
    mapping_file = config.param('default', 'mapping_file', 1, 'filepath'),
    outdir = outdir,
    dummy_outfile = dummy_outfile,
    observed_otus = observed_otus,
    shannon = shannon
    )

    return job

#def upgma_clustering(infile, outfile):
#    outdir = os.path.dirname(os.path.relpath(infile))
#    outdir = os.path.join(outdir, "upgma")
#    outdir = "./" + outdir + "/"
#    
#    base = os.path.basename(outfile)
#    prefix = os.path.splitext(base)[0]
#    
#    job = Job(
#        [infile],
#        [outfile],
#        [
#            ['qiime-dependencies', 'module_qiime-dependencies'],
#            ['qiime', 'module_qiime'],
#            ['python', 'module_python'],
#            ['R', 'module_R'],
#            ['tools', 'module_tools'],
#            ['lapack', 'module_lapack']
#        ]
#    )
#
#    job.command="""
#export OMPI_MCA_mtl=^psm && \\
#upgma_cluster.py \\
#  -i {infile} \\
#  -o {outfile} \\
#  && plotTree.R \\
#  -t {outfile} \\
#  -o {outdir} \\
#  -p {prefix} \\
#  && unset OMPI_MCA_mtl""".format(
#    infile = infile,
#    outfile = outfile,
#    prefix = prefix,
#    outdir = outdir
#  )
#
#    return job

def otu_heatmap(infile, outdir, prefix, n):
    job = Job(
        [infile],
        [outdir + "/" + prefix + "_" + str(n) + ".pdf"],
        [
            ['R', 'module_R'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
featureHeatmap.R \\
  -t {infile} \\
  -o {outdir} \\
  -p {prefix} \\
  -n {n} \\
  -m {mapping_file}""".format(
    infile = infile,
    outdir = outdir,
    prefix = prefix + "_" + n,
    n = n,
    mapping_file = config.param('default', 'mapping_file', 1, 'filepath')
    )

    return job

def blast(infile, outfile, db):
    job = Job(
        [infile],
        [outfile, outfile + ".besthit"],
        [
            ['blast', 'module_blast'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
blastn \\
  -db {db} \\
  -query {infile} \\
  -out {outfile} \\
  -outfmt \\"{outfmt}\\" \\
  -max_target_seqs 1 \\
  -num_threads {num_threads} \\
  && keepBestBlastHit.pl --infile {outfile} > {outfile}.besthit""".format(
    db = db,
    infile = infile,
    outfile = outfile,
    outfmt = config.param('blast', 'outfmt'),
    num_threads = config.param('blast', 'num_threads', 1, 'int')
    )

    return job

def edger_otu_pairwise(abundance, design_file, outdir):

    dummy_outfile = os.path.join(outdir, "DOA.done")

    job = Job(
        [abundance, design_file], 
        [dummy_outfile],
        [
            ['tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
edgerOTUsPairwise.R \\
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
    logfc = config.param('DOA', 'logfc', 1, 'float'),
    fdr = config.param('DOA', 'fdr', 1, 'float'),
    pvalue = config.param('DOA', 'pvalue', 1, 'float')
    )

    return job

def edger_feature_glm(abundance, outdir, mapping_file):

    dummy_outfile = os.path.join(outdir, "DOA_GLM.done")

    job = Job(
        [abundance], 
        [dummy_outfile],
        [
            ['tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )
    
    job.command="""
edgerFeaturesGLM.R \\
  -i {abundance} \\
  -o {outdir} \\
  -m {mapping_file} \\
  -p {pvalue} \\
  -f {fdr} \\
  -l {logfc} \\
  -b {blocks} \\
  -t {treatment} \\
  && touch {dummy_outfile}""".format(
    abundance = abundance,
    outdir = outdir,
    dummy_outfile = dummy_outfile,
    logfc = config.param('DOA_edger', 'logfc', 1, 'float'),
    fdr = config.param('DOA_edger', 'fdr', 1, 'float'),
    pvalue = config.param('DOA_edger', 'pvalue', 1, 'float'),
    blocks = config.param('DOA_edger', 'blocks', 1, 'string'),
    treatment = config.param('DOA_edger', 'treatment', 1, 'string'),
    mapping_file = mapping_file
    )

    return job

def ancom_otus(feature_table, outdir, mapping_file):

    dummy_outfile = os.path.join(outdir, "DOA_ancom.done")

    job = Job(
        [feature_table], 
        [dummy_outfile],
        [
            ['tools', 'module_tools'],
            ['python3', 'module_python3']
        ]
    )
#~/build/meco_tools/RRNATagger/ancom.py -i ./feature_table_filtered_bacteriaArchaea_rarefied_1000.tsv 
# -m mapping_file3.tsv -o ./test_outdir/ -v Treatment -t 16
    job.command="""
export OMP_NUM_THREADS=1 && \\
export MKL_NUM_THREADS=1 && \\
ancom.py \\
  -i {feature_table} \\
  -o {outdir} \\
  -m {mapping_file} \\
  -v {treatment} \\
  -t {num_threads} \\
  && touch {dummy_outfile}""".format(
    feature_table = feature_table,
    outdir = outdir,
    dummy_outfile = dummy_outfile,
    treatment = config.param('DOA_ancom', 'treatment', 1, 'string'),
    num_threads = config.param('DOA_ancom', 'num_threads', 1, 'string'),
    mapping_file = mapping_file
    )

    return job

def round_feature_table(infile, outfile):
    job = Job(
        [infile], 
        [outfile],
        [
            ['R', 'module_R'],
            ['tools', 'module_tools']
        ]
    )
    
    job.command="""
roundFeatureTable.R -i {infile} -o {outfile}""".format(
    infile = infile,
    outfile = outfile
    )
    return job

def cleanup(indir, dummy_infile):
    dummy_outfile = os.path.join(indir, "cleanup.done")
    job = Job(
        [dummy_infile], 
        [dummy_outfile],
        [
            ['pigz', 'module_pigz'],
            ['tools', 'module_tools']
        ]
    )
    
    job.command="""
fastqToGz.pl --indir {indir} --num_threads {num_threads} \\
&& touch {dummy_outfile}""".format(
    indir = indir,
    num_threads = config.param('cleanup', 'num_threads', 1, 'posint'),
    dummy_outfile = dummy_outfile
    )
    return job

def generate_tarball(blast, ini_file_path, project_path, 
                     report_path, dependency_list):
        
    outfile = os.path.join(report_path + ".tar.gz")
    mapping_file = config.param('default', 'mapping_file', 1, 'filepath')

    dependency_list.append(blast)

    job = Job(
        dependency_list,
        [outfile],
        [
            ['R', 'module_R']
        ]
    )
#  && rm -rf *.png
    job.command="""
tar -zcvf {outfile} --dereference --exclude='{project_path}/fastqs' \\
                                  --exclude='{project_path}/../duk' \\
                                  --exclude='{project_path}/fastqs' \\
                                  --exclude='{project_path}/alpha_div/alpha_rarefaction' \\
                                  --exclude='{project_path}/alpha_div/rarefaction' \\
                                  --exclude='{project_path}/feature_tables/rarefactions' \\
                                  {mapping_file} {project_path}/""".format(
        project_path = project_path,
        outfile = outfile,
        mapping_file = mapping_file
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

def rtk_multi_even(infile, outdir, header_first_element, remove_last_col=False):
        
    job = Job(
        [infile],
        [os.path.join(outdir, "rtk.done")],
        [
            ['tools', 'module_tools'],
            ['rtk', 'module_rtk']
        ]
    )
    
    job.command="""
mkdir -p {outdir} && \\""".format(outdir = outdir)
    
#for i in {{0..{perm}}}; do echo -n {n_normalization},; done | sed 's/.$//' > {outdir}/depth_list.txt && \\
#sed -i '1 s/^\\\\t/{header_first_element}\\\\t/' {infile} && \\
    job.command+="""
sed '1 s/^\\\S\\\+\\\\t/\\\\t/' {infile} > {outdir}/tmp.tsv && \\
rtk {mode} \\
  -i {outdir}/tmp.tsv \\
  -o {outdir}/ \\
  -r {perm} -w {perm} \\
  -d {n_normalization} && \\
rm {outdir}/tmp.tsv && \\
touch {outdir}/rtk.done""".format(
        infile = infile,
        outdir = outdir,
        header_first_element = header_first_element,
        n_normalization = config.param('rarefaction', 'n_normalization', 1, 'posint'),
        mode = "memory -ns",#config.param('alphadiv', 'mode', 1, 'string'),
        perm = config.param('rarefaction', 'perm_normalization', 1, 'int')
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

def normalize_cpm_feature_table(abundance_cpm, feature_table, tsv_normalized):
    job = Job(
        [feature_table, abundance_cpm],
        [tsv_normalized],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools'],
            ['R', 'module_R']
        ]
    )

    job.command="""
normalizeFeatureTableCPM.R \\
  -i {feature_table} \\
  -o {tsv_normalized} \\
  -c {abundance_cpm}""".format(
        feature_table = feature_table,
        tsv_normalized = tsv_normalized,
        abundance_cpm = abundance_cpm
    )

    return job

def generate_cog_abundance(abundance_file, cog_term, rpsblast_file, outfile):
    job = Job(
        [abundance_file, rpsblast_file],
        [outfile],
        [
            ['perl', 'module_perl'],
            ['tools', 'module_tools']
        ]
    )

    job.command="""
generateCOGMatrix.pl \\
  --infile {abundance_file} \\
  --cog_term {cog_term} \\
  --cog_file {rpsblast_file} \\
  > {outfile}""".format(
        abundance_file = abundance_file,
        cog_term = cog_term,
        rpsblast_file = rpsblast_file,
        outfile = outfile
    )

    return job
