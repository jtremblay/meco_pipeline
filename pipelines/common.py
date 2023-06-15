#!/usr/bin/env python

################################################################################
# Copyright (C) 2023 INRS - Centre Armand-Frappier
#
# This file is part of CAF Pipelines.
#
# CAF Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CAF Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CAF Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import logging
import os
import re
import socket
import string
import sys
#import hashlib
from hashlib import md5

# Append caf_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))

# CAF Modules
from core.job import *
from core.pipeline import *
from bio.readset import *


log = logging.getLogger(__name__)

# Abstract pipeline gathering common features of all CAF pipelines (readsets, samples, remote log, etc.)
class CAFPipeline(Pipeline):

    def __init__(self):
        self.version = open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "VERSION"), 'r').read().split('\n')[0]

        # Add pipeline specific arguments
        self.argparser.description = "Version: " + self.version + "\n\nFor more documentation, visit our website: https://bitbucket.org/jtremblay514/caf_pipeline_public/"
        self.argparser.add_argument("-v", "--version", action="version", version="caf_pipeline " + self.version, help="show the version information and exit")

        super(CAFPipeline, self).__init__()

    @property
    def readsets(self):
        return self._readsets

    @property
    def samples(self):
        if not hasattr(self, "_samples"):
            self._samples = list(collections.OrderedDict.fromkeys([readset.sample for readset in self.readsets]))
        return self._samples

    #def mugqic_log(self):
    #    server = "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi"
    #    listName = {}
    #    for readset in self.readsets:
    #        if listName.has_key(readset.sample.name) :
    #            listName[readset.sample.name]+="."+readset.name
    #        else: 
    #            listName[readset.sample.name]=readset.sample.name+"."+readset.name
    #    request = \
    #        "hostname=" + socket.gethostname() + "&" + \
    #        "ip=" + socket.gethostbyname(socket.gethostname()) + "&" + \
    #        "pipeline=" + self.__class__.__name__ + "&" + \
    #        "steps=" + ",".join([step.name for step in self.step_range]) + "&" + \
    #        "samples=" + str(len(self.samples)) + "&" + \
    #        "AnonymizedList=" + ",".join([md5.md5(self.__class__.__name__ + "." + value).hexdigest() for key, value in listName.iteritems()])

#        print("""
#{separator_line}
## Call home with pipeline statistics
#{separator_line}
#wget "{server}?{request}" --quiet --output-document=/dev/null
#""".format(separator_line = "#" + "-" * 79, server=server, request=request))


    def submit_jobs(self):
        super(CAFPipeline, self).scheduler.submit(self)
    #    #if self.jobs and self.args.job_scheduler in ["pbs", "batch"]:
    #    #    self.mugqic_log()


# Abstract pipeline gathering common features of all Illumina sequencing pipelines (trimming, etc.)
# Specific steps must be defined in Illumina children pipelines.
#class Illumina(CAFPipeline):
#
#    def __init__(self):
#        self.argparser.add_argument("-r", "--readsets", help="readset file", type=file)
#        super(Illumina, self).__init__()
#
#    @property
#    def readsets(self):
#        if not hasattr(self, "_readsets"):
#            if self.args.readsets:
#                self._readsets = parse_illumina_readset_file(self.args.readsets.name)
#            else:
#                self.argparser.error("argument -r/--readsets is required!")
#        return self._readsets
#
#    @property
#    def run_type(self):
#        run_types = [readset.run_type for readset in self.readsets]
#        if len(set(run_types)) == 1 and re.search("^(PAIRED|SINGLE)_END$", run_types[0]):
#            return run_types[0]
#        else:
#            raise Exception("Error: readset run types " + ",".join(["\"" + run_type + "\"" for run_type in run_types]) +
#            " are invalid (should be all PAIRED_END or all SINGLE_END)!")
#
