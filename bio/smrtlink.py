#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015, 2016, 2017, 2018, 2019 Julien Tremblay

# This file is part of MECO Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MECO Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# MECO Modules
from core.config import *
from core.job import *


def bax5_to_subreads(indir, outdir):

    job = Job(
        [os.path.join(indir, "generate_fofns.done")],
        [(outdir + "_generate_bam_subreads.done")],
        [
            ['smrtlink', 'module_smrtlink'],
            ['meco_tools', 'module_tools']
        ]
    )

#bax2bam -f fofn_2.txt -o subreads/  --subread --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,MergeQV,SubstitutionQV,PulseWidth,SubstitutionTag
    job.command="""\
pacBioBax5ToSubreads.pl \\
 --indir {indir} \\
 --outdir {outdir} && \\
 touch {done_file}""".format(
    indir = indir,
    outdir = outdir,
    done_file = outdir + "_generate_bam_subreads.done"
    )

    return job

