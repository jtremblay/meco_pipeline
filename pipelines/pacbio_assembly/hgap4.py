#!/usr/bin/env python

################################################################################
# Copyright (C) 2023 INRS Centre Armand-Frappier
#
# This file is part of MECO Pipelines.
#
# MECO Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MECO Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MECO Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import logging
import os
import sys

# Append meco_pipelines directory to Python library path
#sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MECO Modules
from core.config import *
from core.job import *
from bio.readset import *

from bio import pacbio_tools
from bio import smrtlink
from pipelines import common

log = logging.getLogger(__name__)

class HGAP4(common.MECOPipeline):
    """
    PacBio Assembly Pipeline - HGAP4
    Written by Julien Tremblay
    ========================

    Contigs assembly with PacBio reads is done using what is refer as the HGAP workflow.
    Briefly, raw subreads generated from raw .ba(s|x).h5 PacBio data files are filtered for quality.
    A subread length cutoff value is extracted from subreads, depending on subreads distribution,
    and used into the preassembly (aka correcting step) (BLASR) step which consists of aligning
    short subreads on long subreads.
    Since errors in PacBio reads is random, the alignment of multiple short reads on longer reads
    allows to correct sequencing error on long reads.
    These long corrected reads are then used as seeds into assembly (Celera assembler) which gives contigs.
    These contigs are then *polished* by aligning raw reads on contigs (BLASR) that are then processed
    through a variant calling algorithm (Quiver) that generates high quality consensus sequences
    using local realignments and PacBio quality scores.

    Prepare your readset file as follow (all separated by a tab). EstimatedGenomeSize field should be the same in each row. 
    Sample    Smartcell     EstimatedGenomesize   BAS                BAX
    LB501T    1             6500000               <leave empty>      ./raw_reads/m170517_223644_42158_c101204432550000001823285710241730_s1_p0.3.bax.h5
    LB501T    2             6500000               <leave empty>      ./raw_reads/m170517_223644_42158_c101204432550000001823285710241730_s1_p0.2.bax.h5
    LB501T    3             6500000               <leave empty>      ./raw_reads/m170517_223644_42158_c101204432550000001823285710241730_s1_p0.1.bax.h5
    TODO: replace Smartcell with FileNumber.
    """

    def __init__(self):
        self.argparser.add_argument("-r", "--readsets", help="readset file", type=file)
        super(HGAP4, self).__init__()

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            if self.args.readsets:
                self._readsets = parse_pacbio_readset_file(self.args.readsets.name)
            else:
                self.argparser.error("argument -r/--readsets is required!")

        return self._readsets

    def bax5_to_subreads(self):
        """
        Filter reads and subreads based on their length and QVs, using smrtpipe.py (from the SmrtAnalysis package).

        1. fofnToSmrtpipeInput.py
        2. modify RS_Filtering.xml files according to reads filtering values entered in .ini file
        3. smrtpipe.py with filtering protocol
        4. prinseq-lite.pl: write fasta file based on fastq file

        Informative run metrics such as loading efficiency, read lengths and base quality are generated in this step as well.
        """


        jobs = []

        for sample in self.samples:
            if not os.path.exists(os.path.join(sample.name, "fofns")):
                os.makedirs(os.path.join(sample.name, "fofns"))
            if not os.path.exists(os.path.join(sample.name, "subreads")):
                os.makedirs(os.path.join(sample.name, "subreads"))
            
            log.info("Sample: " + sample.name)
            #sample_nb_base_pairs = sum([readset.nb_base_pairs for readset in sample.readsets])
            #log.info("nb_base_pairs: " + str(sample_nb_base_pairs))
            estimated_genome_size = sample.readsets[0].estimated_genome_size
            log.info("estimated_genome_size: " + str(estimated_genome_size))
            #estimated_coverage = sample_nb_base_pairs / estimated_genome_size
            #log.info("estimated_coverage: " + str(estimated_coverage))
            job = pacbio_tools.generate_fofns(
                sample.name,
                os.path.join("readset.tsv"),
                os.path.join(sample.name, "fofns")
            )
            job.name = "generate_fofns"
            job.subname = "generate_fofns"
            jobs.append(job)
            
            job = smrtlink.bax5_to_subreads(
                os.path.join(sample.name, "fofns"),
                os.path.join(sample.name, "subreads", sample.name)
            )
            job.name = "bax5_to_subreads"
            job.subname = "bax5_to_subreads"
            jobs.append(job)



            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)

                log.info("COVERAGE_CUTOFF: " + coverage_cutoff + "_X_coverage")

                #jobs.append(concat_jobs([
                #    Job(command="mkdir -p " + os.path.join(coverage_directory, "preassembly")),
                #    pacbio_tools.get_cutoff(
                #        os.path.join(sample.name, "filtering", "data", "filtered_subreads.fasta"),
                #        os.path.join(sample.name, "filtering", "data", "number_of_bp.txt"),#estimated_coverage,
                #        estimated_genome_size,
                #        coverage_cutoff,
                #        os.path.join(coverage_directory, "preassemblyMinReadSize.txt")
                #    )
                #], name="pacbio_tools_get_cutoff." + sample.name + ".coverage_cutoff" + cutoff_x, subname="pacbio_tools_get_cutoff"))

        return jobs


    @property
    def steps(self):
        return [
            self.bax5_to_subreads, #1
        ]

if __name__ == '__main__':
    HGAP4()

HGAP4().submit_jobs()
