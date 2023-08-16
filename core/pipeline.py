#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import hashlib
import logging
import os
import re
from collections import defaultdict
import json
import sys

# MECO Modules
from core.config import *
from core.scheduler import *
from core.step import *

log = logging.getLogger(__name__)

class Pipeline(object):
    
    def __init__(self):
        #sys.stderr.write("object: " + str(object) + "\n")
        def tree(): return defaultdict(tree)
        self._args = self.argparser.parse_args()

        logging.basicConfig(level=getattr(logging, self.args.log.upper()))
        config.parse_files(self.args.config)

        self._output_dir = os.path.abspath(self.args.output_dir)
        self._scheduler = create_scheduler(self.args.job_scheduler)
        self._force_jobs = self.args.force
        self._output_files = []

        # For JSON pipeline path.
        self._jobs_dict = dict()
        self._curr_el = 0
        self._tree = tree()
        self._nodes_tree = tree()

        step_counter = collections.Counter(self.steps)
        duplicated_steps = [step.__name__ for step in step_counter if step_counter[step] > 1]
        if duplicated_steps:
            raise Exception("Error: pipeline contains duplicated steps: " + ", ".join(duplicated_steps) + "!")
        else:
            # Will loop in Step class for each step of the pipeline selected or not. So yes __init__ is executed each time.
            self._step_list = [Step(step) for step in self.steps]

        if re.search("^\d+([,-]\d+)*$", self.args.steps):
            self._step_range = [self.step_list[i - 1] for i in parse_range(self.args.steps)]
        else:
            raise Exception("Error: step range \"" + self.args.steps +
                "\" is invalid (should match \d+([,-]\d+)*)!")
        #sys.stderr.write("self.create_jobs in Pipeline class: " + str(self.create_jobs) + "\n")
        self.create_jobs()
        self.generate_json()

    # Pipeline command line arguments parser

    @property
    def argparser(self):
        if not hasattr(self, "_argparser"):
            # Create ArgumentParser with numbered step list as epilog
            self._argparser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, epilog="Steps:\n" + "\n".join([str(idx + 1) + "- " + step.__name__ for idx, step in enumerate(self.steps)]))

            # Common options for all pipelines
            self._argparser.add_argument("-c", "--config", help="config INI-style file", nargs="+", type=argparse.FileType('r'), required=True)
            self._argparser.add_argument("-s", "--steps", help="step range e.g. '1-5', '3,6,7', '2,4-8'", required=True)
            self._argparser.add_argument("-o", "--output-dir", help="output directory (default: current)", default=os.getcwd())
            self._argparser.add_argument("-j", "--job-scheduler", help="job scheduler type (default: slurm)", choices=["slurm", "torque", "batch", "daemon", "sge"], default="slurm")
            self._argparser.add_argument("-f", "--force", help="force creation of jobs even if up to date (default: false)", action="store_true")
            self._argparser.add_argument("-l", "--log", help="log level (default: info)", choices=["debug", "info", "warning", "error", "critical"], default="info")
            self._argparser.add_argument("-z", "--json", help="generate pipeline path in json format", default=sys.stdout, type=argparse.FileType('w'), required=False)

        return self._argparser

    # Pipeline command line arguments
    @property
    def args(self):
        return self._args

    @property
    def output_dir(self):
        return self._output_dir

    @property
    def scheduler(self):
        return self._scheduler

    @property
    def force_jobs(self):
        return self._force_jobs

    @property
    def steps(self):
        # Needs to be defined in pipeline child class
        raise NotImplementedError

    @property
    def step_list(self):
        return self._step_list

    @property
    def step_range(self):
        return self._step_range
    
    @property
    def get_output_files(self):
        return self._output_files

    @property
    def jobs(self):
        jobs = []
        for step in self.step_range:
            jobs.extend(step.jobs)
        return jobs

    def dependency_jobs(self, current_job):
        dependency_jobs = []
        dependency_input_files = set()
        for step in self.step_range:
            for step_job in step.jobs:
                # If current job input files intersect with step job output files, step job is a dependency
                #sys.stderr.write("step_job.output_files:\n " + ",".join(step_job.output_files))
                self._output_files.extend(step_job.output_files)
                shared_files = set(current_job.input_files).intersection(set(step_job.output_files))
                if shared_files:
                    dependency_jobs.append(step_job)
                    dependency_input_files.update(shared_files)

        # Check if job input files not found in dependencies are on file system
        missing_input_files = set()
        # Add current_job.output_files in case of "... && ..." command
        # where first command output becomes second command input
        for remaining_input_file in set(current_job.input_files).difference(dependency_input_files).difference(set(current_job.output_files)):
            if not os.path.isfile(current_job.abspath(remaining_input_file)):
                missing_input_files.add(remaining_input_file)
        if missing_input_files:
            raise Exception("Error: missing input files for job " + current_job.name + ": " +
                ", ".join(missing_input_files) + " neither found in dependencies nor on file system!")

        return dependency_jobs

    def create_jobs(self):
        for step in self.step_range:
            log.info("Create jobs for step " + step.name + "...")
            jobs = step.create_jobs()
            for job in jobs:
                # Job name is mandatory to create job .done file name
                if not job.name:
                    raise Exception("Error: job \"" + job.command + "\" has no name!")

                # populate job dict.
                self._jobs_dict[job.name] = self._curr_el
                
                # Is that a job array type of job?
                if(job.job_array_num_task > 0):
                    # Job .done file name contains the command checksum.
                    # Thus, if the command is modified, the job is not up-to-date anymore
                    missing_tasks_num = []
                    job.done_job_array = os.path.join("job_output", step.name, job.name + "." + hashlib.md5(job.command.encode('utf-8')).hexdigest() + ".nrc.done")
                    for curr_task_num in range(0,job.job_array_num_task):
                        #sys.stderr.write("[DEBUG] curr_task_num (in pipeline.py): " + str(curr_task_num) + "\n")
                        job.done = os.path.join("job_output", step.name, job.name + "." + hashlib.md5(job.command.encode('utf-8')).hexdigest() + ".nrc.done")
                        job.output_dir = self.output_dir
                        job.dependency_jobs = self.dependency_jobs(job)
                        if not self.force_jobs and job.is_up2date(curr_task_num):
                            log.info("Job array " + job.name + "_{:07d}".format(curr_task_num) + " up to date... skipping")
                        else:
                            missing_tasks_num.append(curr_task_num)

                    missing_tasks_string = ','.join(str(e) for e in missing_tasks_num)
                    job.set_job_array_num_task_final(str(missing_tasks_string))
                    if len(missing_tasks_num) > 0:
                        step.add_job(job)
                
                # else if normal job
                else:
                    # Job .done file name contains the command checksum.
                    # Thus, if the command is modified, the job is not up-to-date anymore.
                    job.done = os.path.join("job_output", step.name, job.name + "." + hashlib.md5(job.command.encode('utf-8')).hexdigest() + ".nrc.done")
                    job.output_dir = self.output_dir
                    job.dependency_jobs = self.dependency_jobs(job)
                    if not self.force_jobs and job.is_up2date():
                        log.info("Job " + job.name + " up to date... skipping")
                    else:
                        step.add_job(job)
    
                        # Search for CPU and PPN
                        cpu_param = config.param(job.subname, 'cluster_cpu')
                        if re.match("-N (\d+) --ntasks-per-node=(\d+)", cpu_param):
                            m = re.search("-N (\d+) --ntasks-per-node=(\d+)", cpu_param)
                            cpu = int(m.group(1)) * int(m.group(2))
                        elif re.match("-l nodes=(\d+):ppn=(\d+)", cpu_param):
                            m = re.search("-l nodes=(\d+):ppn=(\d+)", cpu_param)
                            cpu = int(m.group(1)) * int(m.group(2))
                        elif re.match("-N \d+ -n \d+", cpu_param):
                            m = re.search("-N (\d+) -n (\d+)", cpu_param)
                            cpu = int(m.group(1)) * int(m.group(2))
                        elif re.search("-pe \S+ (\d+) -l res_cpus=(\d+)", cpu_param):
                            #-pe dev 1 -l res_cpus=6
                            m = re.search("-pe \S+ (\d+) -l res_cpus=(\d+)", cpu_param)
                            cpu = int(m.group(1)) * int(m.group(2))
                        else:
                            sys.stderr.write("[DEBUG]: Could not parse cpu_param parameter: " + cpu_param + "\n")
                            exit(1);

                        # Nodes
                        self._nodes_tree[self._curr_el]['name'] = job.name
                        self._nodes_tree[self._curr_el]['group'] = step.name
                        self._nodes_tree[self._curr_el]['value'] = str(cpu)
    
                        # Links
                        self._tree[self._curr_el]['job.name'] = job.name
                        self._tree[self._curr_el]['job.number'] = self._curr_el
                        i = 0
                        for curr_dep_job in job.dependency_jobs:
                            self._tree[self._curr_el]['job.dep.name'][str(i)] = curr_dep_job.name
                            self._tree[self._curr_el]['job.dep.number'][str(i)] = self._jobs_dict[curr_dep_job.name]
                            i += 1

                        self._curr_el += 1


            log.info("Step " + step.name + ": " + str(len(step.jobs)) + " job" + ("s" if len(step.jobs) > 1 else "") + " created" + ("" if step.jobs else "... skipping") + "\n")
        log.info("TOTAL: " + str(len(self.jobs)) + " job" + ("s" if len(self.jobs) > 1 else "") + " created" + ("" if self.jobs else "... skipping") + "\n")

    def generate_json(self):
        #self.scheduler.submit(self)
        
        if self.args.json:
            # Print pipeline path to JSON file
            self._json_file = self.args.json
            self._json_file = os.path.abspath(self._json_file.name)
            sys.stderr.write("JSON FILE: " + self._json_file + "\n")
            self._json_fh = open(self._json_file, 'w')

            # Print tree JSONesque style
            #self._json_fh.write(json.dumps(self._tree))
                    
            self._json_fh.write('{\n   "nodes":[\n')
            x = 0
            for k in self._nodes_tree.keys():
                node_name = ""
                node_group = ""
                node_value = ""

                for cat1 in self._nodes_tree[k]:

                    #sys.stderr.write("cat1:" + cat1 + "\n")

                    if(cat1 == "name"):
                        node_name = str(self._nodes_tree[k][cat1])
                    if(cat1 == "group"):
                        node_group = str(self._nodes_tree[k][cat1])
                    if(cat1 == "value"):
                        node_value = str(self._nodes_tree[k][cat1])

                if(x > self._curr_el - 2):  
                    self._json_fh.write('      {"name":"' + node_name + '", "group":"' + node_group + '", "value":"' + node_value + '"}\n' )
                else:
                    self._json_fh.write('      {"name":"' + node_name + '", "group":"' + node_group + '", "value":"' + node_value + '"},\n' )

                    x += 1
            self._json_fh.write('   ],\n   "links":[\n')
            
            #""" print a tree """
            y = 0
            for k in self._tree.keys():
                for cat1 in self._tree[k]:
                    if(cat1 == "job.dep.number"):
                        for number in self._tree[k][cat1]:
                            y += 1

            x = 0
            for k in self._tree.keys():
                job_number = ""
                job_name = ""
                job_dep_numbers = []
                job_dep_names = []

                for cat1 in self._tree[k]:

                    if(cat1 == "job.dep.number"):
                        for number in self._tree[k][cat1]:
                            job_dep_numbers.append(str(self._tree[k][cat1][number]))
                
                    if(cat1 == "job.dep.name"):
                        for name in self._tree[k][cat1]:
                            job_dep_names.append(self._tree[k][cat1][name])

                    if(cat1 == "job.number"):
                        job_number = str(self._tree[k][cat1])
                    
                    if(cat1 == "job.name"):
                        job_name = self._tree[k][cat1]

                for job_dep_number in job_dep_numbers:
                    job_dep_name = job_dep_names.pop(0)
                    if(x > y - 2):  
                        self._json_fh.write('      {"source":' + job_dep_number + ', "source_name":"' + job_dep_name + '", "target":' + job_number + ', "target_name":"' + job_name + '"}\n' )
                    else:
                        self._json_fh.write('      {"source":' + job_dep_number + ', "source_name":"' + job_dep_name + '", "target":' + job_number + ', "target_name":"' + job_name + '"},\n' )
                    x += 1
            
            self._json_fh.write('   ]\n}\n')
            self._json_fh.close()


# Return a range list given a string.
# e.g. parse_range('1,3,5-12') returns [1, 3, 5, 6, 7, 8, 9, 10, 11, 12]
def parse_range(astr):
    result = set()
    for part in astr.split(','):
        x = part.split('-')
        result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)
