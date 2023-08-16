#!/usr/bin/env python

# Python Standard Modules
import collections
import logging
import os

# MECO Modules
from core.config import *

log = logging.getLogger(__name__)

class Job:

    def __init__(self, input_files=[], output_files=[], module_entries=[], name="", command="", subname="", job_array_num_task=0):
        # Remove undefined input/output files if any
        self._input_files = list(filter(None, input_files))
        self._output_files = list(filter(None, output_files))

        # Retrieve modules from config, removing duplicates but keeping the order
        self._modules = list(collections.OrderedDict.fromkeys([config.param(section, option) for section, option in module_entries]))

        self._name = name
        self._subname = subname
        self._job_array_num_task = job_array_num_task
        self._command = command
        self._job_array_num_task_final = "unset"

    def show(self):
        print("Job: input_files: " + \
            ", ".join(self.input_files))

    @property
    def id(self):
        return self._id
    
    @id.setter
    def id(self, value):
        self._id = value

    # Added by JT. subname is to defined what cluster paramters must be picked up from the ini file for a particular job.
    @property
    def subname(self):
        return self._subname
    
    @subname.setter
    def subname(self, value):
        self._subname = value
    
    # Added by JT for job array support.
    @property
    def job_array_num_task(self):
        return self._job_array_num_task
    
    @job_array_num_task.setter
    def job_array_num_task(self, value):
        self._job_array_num_task = value
    
    # Added by JT for job array support.
    @property
    def job_array_num_task_final(self):
        return self._job_array_num_task_final
    
    @job_array_num_task_final.setter
    def job_array_num_task_final(self, value):
        self._job_array_num_task_final = value
    
    # Added by JT for job array support.
    #@job_array_num_task_final.setter
    def set_job_array_num_task_final(self, job_array_num_task_final):
        #sys.stderr.write("job_array_num_task_final" + job_array_num_task_final + "\n")
        self._job_array_num_task_final = job_array_num_task_final
    

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, value):
        self._name = value

    @property
    def output_dir(self):
        return self._output_dir
    
    @output_dir.setter
    def output_dir(self, value):
        self._output_dir = value

    @property
    def input_files(self):
        return self._input_files
    
    @input_files.setter
    def input_files(self, value):
        self._input_files = value

    @property
    def output_files(self):
        return self._output_files
    
    @output_files.setter
    def output_files(self, value):
        self._output_files = value

    @property
    def done(self):
        return self._done
    
    @done.setter
    def done(self, value):
        self._done = value
    
    @property
    def done_job_array(self):
        return self._done_job_array
    
    @done_job_array.setter
    def done_job_array(self, value):
        self._done_job_array = value

    @property
    def dependency_jobs(self):
        return self._dependency_jobs
    
    @dependency_jobs.setter
    def dependency_jobs(self, value):
        self._dependency_jobs = value

    @property
    def modules(self):
        return self._modules
    
    @modules.setter
    def modules(self, value):
        self._modules = value

    @property
    def command(self):
        return self._command

    @command.setter
    def command(self, value):
        self._command = value

    @property
    def command_with_modules(self):
        command = self.command
        if self.modules:
            command = "module load " + " ".join(self.modules) + " && \\\n" + command
        return command
    
    @command_with_modules.setter
    def command_with_modules(self, value):
        self._done = value
    
    @property
    def command_with_modules_and_unload(self):
        command = self.command
        if self.modules:
            command = "module load " + " ".join(self.modules) + " && \\\n" + command + " ; \\\n module unload " + " ".join(self.modules)
        return command

    def abspath(self, file):
        tmp_file = os.path.expandvars(file)
        if not os.path.isabs(tmp_file):
            # File path is relative to the job output directory
            tmp_file = os.path.normpath(os.path.join(self.output_dir, tmp_file))
        return tmp_file

    def is_up2date(self, job_array_num_task=None):
        # If job has dependencies, job is not up to date
        if self.dependency_jobs:
            return False

        # Retrieve absolute paths for .done, input and output files to avoid redundant OS function calls
        abspath_done = self.abspath(self.done)
        abspath_input_files = [self.abspath(input_file) for input_file in self.input_files]
        abspath_output_files = [self.abspath(output_file) for output_file in self.output_files]
    
        # Here job array jobs will have to set of outfiles: 1 tsv(real files) + 1 done (dummy files).
        # Here in is_up2date, only .done file will be considered as it is touched once the 
        # array[i] blast job has completed.
        if(job_array_num_task is not None):
            abspath_done_job_array = self.abspath(self.done_job_array)
            abspath_input_files = [abspath_input_files[job_array_num_task]]
            abspath_output_files = [abspath_output_files[job_array_num_task]]
        
            #sys.stderr.write("abspath_done: " + str(abspath_done) + "\n")
            #sys.stderr.write("abspath_input_files: " + str(abspath_input_files) + "\n")
            #sys.stderr.write("abspath_output_files: " + str(abspath_output_files) + "\n")
           
            # Because job array: just check in input and output files in which a dummy .done should have been included
            # If any .done, input or output file is missing, job is not up to date
            #for file in [abspath_done_job_array] + abspath_input_files + abspath_output_files:
            for file in abspath_input_files + abspath_output_files:
                if not os.path.isfile(file):
                    #sys.stderr.write("absent: " + str(file) + "\n")
                    return False
                #else:
                #    sys.stderr.write("present: " + str(file) + "\n")

            # Retrieve latest input file modification time i.e. maximum stat mtime
            latest_input_time = max([os.stat(input_file).st_mtime for input_file in abspath_input_files])
    
            # Same with earliest output file modification time
            earliest_output_time = min([os.stat(output_file).st_mtime for output_file in abspath_output_files])
    
            # If any input file is strictly more recent than all output files, job is not up to date
            if latest_input_time > earliest_output_time:
                return False
    
            # If all previous tests passed, job is up to date
            return True

        else:

            # If any .done, input or output file is missing, job is not up to date
            for file in [abspath_done] + abspath_input_files + abspath_output_files:
                if not os.path.isfile(file):
                    return False

            # Preprocess to remove empty elements (python 3 req).
            # Retrieve latest input file modification time i.e. maximum stat mtime
            latest_input_time = max([os.stat(input_file).st_mtime for input_file in abspath_input_files])
    
            # Same with earliest output file modification time
            earliest_output_time = min([os.stat(output_file).st_mtime for output_file in abspath_output_files])
    
            # If any input file is strictly more recent than all output files, job is not up to date
            if latest_input_time > earliest_output_time:
                return False
    
            # If all previous tests passed, job is up to date
            return True


# Create a new job by concatenating a list of jobs together
def concat_jobs(jobs, name="", subname=""):

    # Merge all input/output files and modules
    input_files = []
    output_files = []
    modules = []
    for job_item in jobs:
        input_files.extend([input_file for input_file in job_item.input_files if input_file not in input_files and input_file not in output_files])
        output_files.extend([output_file for output_file in job_item.output_files if output_file not in output_files])
        modules.extend([module for module in job_item.modules if module not in modules])

    job = Job(input_files, output_files, name=name, subname=subname)
    job.modules = modules

    # Merge commands
    job.command = " && \\\n".join([job_item.command for job_item in jobs])

    return job

# Create a new job by piping a list of jobs together
def pipe_jobs(jobs, name=""):

    job = Job(jobs[0].input_files, jobs[-1].output_files, name=name)

   # Merge all modules
    modules = []
    for job_item in jobs:
        modules.extend(job_item.modules)

    # Remove duplicates if any, keeping the order
    modules = list(collections.OrderedDict.fromkeys([module for module in modules]))
    job.modules = modules

    # Merge commands
    job.command = " | \\\n".join([job_item.command for job_item in jobs])

    return job

