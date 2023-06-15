#!/usr/bin/env python

# Python Standard Modules
import re
import sys

class Step:
    def __init__(self, create_jobs):
        #sys.stderr.write("Inside __init__ of class Step : " + str(self) + "\n")
        #sys.stderr.write("Inside __init__ of class Step : " + str(create_jobs) + "\n")
        #sys.stderr.write("Inside __init__ of class Step : " + str(create_jobs.__name__) + "\n")
        # Step name is used in Bash $JOB_ID variable, hence only alphanumeric and "_" characters are allowed
        step_name = create_jobs.__name__
        # here __name__ is equal to the step variable name (i.e. in the array where we defined all the steps).
        if re.search("^[a-zA-Z]\w+$", step_name):
            self._name = step_name
        else:
            raise Exception("Error: step name \"" + step_name +
                "\" is invalid (should match [a-zA-Z][a-zA-Z0-9_]+)!")

        self._name = step_name
        self._create_jobs = create_jobs
        self._jobs = []

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, value):
        self._name = value

    @property
    def create_jobs(self):
        return self._create_jobs
    
    @create_jobs.setter
    def create_jobs(self, value):
        self._creat_jobs = value

    @property
    def jobs(self):
        return self._jobs
    
    @jobs.setter
    def jobs(self, value):
        self._jobs = value

    #@property 
    def add_job(self, job):
        self.jobs.append(job)
        job.id = self.name + "_" + str(len(self.jobs)) + "_JOB_ID"
    
    #@add_job.setter
    def set_add_job(self, value):
        self._add_job = value
