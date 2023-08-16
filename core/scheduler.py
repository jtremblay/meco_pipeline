#!/usr/bin/env python

# Python Standard Modules
import datetime
import os

# MECO Modules
from core.config import *

# Output comment separator line
separator_line = "#" + "-" * 79

def create_scheduler(type):
    if type == "slurm":
        return SlurmScheduler()
    elif type == "torque":
        return TorqueScheduler()
    elif type == "sge":
        return SGEScheduler()
    elif type == "batch":
        return BatchScheduler()
    else:
        raise Exception("Error: scheduler type \"" + type + "\" is invalid!")

class Scheduler:
    def submit(self, pipeline):
        # Needs to be defined in scheduler child class
        raise NotImplementedError

    def print_header(self, pipeline):
        print(
"""#!/bin/bash

{separator_line}
# {pipeline.__class__.__name__} {scheduler.__class__.__name__} Job Submission Bash script
# Created on: {datetime}
# Steps:
{steps}
{separator_line}

OUTPUT_DIR={pipeline.output_dir}
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/{pipeline.__class__.__name__}_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR"""
            .format(
                separator_line=separator_line,
                pipeline=pipeline,
                scheduler=self,
                steps="\n".join(["#   " + step.name + ": " + str(len(step.jobs)) + " job" + ("s" if len(step.jobs) > 1 else "" if step.jobs else "... skipping") for step in pipeline.step_range]) + \
                "\n#   TOTAL: " + str(len(pipeline.jobs)) + " job" + ("s" if len(pipeline.jobs) > 1 else "" if pipeline.jobs else "... skipping"),
                datetime=datetime.datetime.now()
            )
        )

    def print_step(self, step):
        print("""
{separator_line}
# STEP: {step.name}
{separator_line}
STEP={step.name}
mkdir -p $JOB_OUTPUT_DIR/$STEP""".format(separator_line=separator_line, step=step)
        )

class SGEScheduler(Scheduler):

    def gen_ranges(self, lst):
        s = e = None
        for i in sorted(lst):
            if s is None:
                s = e = i
            elif i == e or i == e + 1:
                e = i
            else:
                yield (s, e)
                s = e = i
    
        if s is not None:
            yield (s, e)

    def submit(self, pipeline):
        self.print_header(pipeline)
        # Here also print job scripts dir
        print("""
SCRIPTS_DIR=$OUTPUT_DIR/job_scripts
mkdir -p $SCRIPTS_DIR
        """
        )
        for step in pipeline.step_range:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    #sys.stderr.write("job.job_array_num_task_final: " + job.job_array_num_task_final + "\n")
                    if job.job_array_num_task > 0:
                        #sys.stderr.write("[DEBUG] job.job_array_num_task_final:" + job.job_array_num_task_final + "\n")
                        # format num_task output using the gen_ranges() function defined above.
                        job_array_list = list(job.job_array_num_task_final.split(","))
                        job_array_list = map(int, job_array_list)
                        # For SGE, increment all 
                        job_array_list = [x+1 for x in job_array_list]
                        #sys.stderr.write("[DEBUG] job_array_list:" + str(job_array_list) + "\n")
                        job_array_list_formatted = repr(','.join(['%d' % s if s == e else '%d-%d' % (s, e) for (s, e) in self.gen_ranges(job_array_list)]))
                        cluster_job_array = "-t " + str(job_array_list_formatted)
                    else:
                        cluster_job_array = ""

                    if job.dependency_jobs:
                        # Chunk JOB_DEPENDENCIES on multiple lines to avoid lines too long
                        max_dependencies_per_line = 50
                        dependency_chunks = [job.dependency_jobs[i:i + max_dependencies_per_line] for i in range(0, len(job.dependency_jobs), max_dependencies_per_line)]
                        job_dependencies = "JOB_DEPENDENCIES=" + config.param('default', 'cluster_dependency_sep', 1, 'string').join(["$" + dependency_job.id for dependency_job in dependency_chunks[0]])
                        for dependency_chunk in dependency_chunks[1:]:
                            job_dependencies += "\nJOB_DEPENDENCIES=$JOB_DEPENDENCIES" + config.param('default', 'cluster_dependency_sep', 1, 'string') + config.param('default', 'cluster_dependency_sep', 1, 'string').join(["$" + dependency_job.id for dependency_job in dependency_chunk])
                    else:
                        job_dependencies = "JOB_DEPENDENCIES="

                    print("""
{separator_line}
# JOB: {job.id}: {job.name}
{separator_line}
JOB_NAME={job.name}
{job_dependencies}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
rm -f $SCRIPTS_DIR/{job.name}.sh
                          """.format(
                            job=job,
                            job_dependencies=job_dependencies,
                            separator_line=separator_line
                        )
                    )

                    cluster_other_arg      = config.param(job.name, 'cluster_other_arg')      
                    cluster_work_dir_arg   = config.param(job.name, 'cluster_work_dir_arg')   
                    cluster_output_dir_arg = config.param(job.name, 'cluster_output_dir_arg') 
                    cluster_job_name_arg   = config.param(job.name, 'cluster_job_name_arg')   
                    cluster_walltime       = config.param(job.subname, 'cluster_walltime')
                    cluster_queue          = config.param(job.subname, 'cluster_queue')
                    cluster_cpu            = config.param(job.subname, 'cluster_cpu')
                    cluster_pmem           = config.param(job.subname, 'cluster_pmem')
                    #cluster_job_array      = config.param(job.subname, 'cluster_job_array')

                    ##$ {cluster_submit_cmd} + " " + \
                    cmd = """echo "#!/bin/bash
#$ {cluster_other_arg}
#$ {cluster_work_dir_arg} $OUTPUT_DIR
#$ {cluster_output_dir_arg} $JOB_OUTPUT
#$ {cluster_job_name_arg} $JOB_NAME
#$ {cluster_walltime}
#$ {cluster_queue}
#$ {cluster_pmem}
#$ {cluster_cpu}
#$ {cluster_job_array}""".format(
                            #cluster_submit_cmd     = cluster_submit_cmd, 
                            cluster_other_arg      = cluster_other_arg,
                            cluster_work_dir_arg   = cluster_work_dir_arg,
                            cluster_output_dir_arg = cluster_output_dir_arg,
                            cluster_job_name_arg   = cluster_job_name_arg,
                            cluster_walltime       = cluster_walltime,
                            cluster_queue          = cluster_queue,
                            cluster_pmem           = cluster_pmem,
                            cluster_cpu            = cluster_cpu,
                            cluster_job_array      = cluster_job_array
                        )
                    
                    if job.dependency_jobs:
                        cmd += "\n#$ " + config.param(step.name, 'cluster_dependency_arg') + " $JOB_DEPENDENCIES" + "\n"

                    cmd += \
"""
rm -f $JOB_DONE
umask 0002
set -o pipefail
{install_home_and_modules}
{job.command_with_modules}
MECO_STATE=\$?
echo MECO_exitStatus:\$MECO_STATE
if [ \$MECO_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MECO_STATE" > $SCRIPTS_DIR/{job.name}.sh && 
{cluster_submit_cmd} $SCRIPTS_DIR/{job.name}.sh""".format(
        job = job,
        install_home_and_modules = "export TMPDIR=" + config.param('default', 'tmpdir', 1, 'string') + " && export INSTALL_HOME=" + config.param('default', 'install_home', 1, 'string') + " && source " + config.param('default', 'env_modules', 1, 'filepath') + " && module use \$INSTALL_HOME/modulefiles",
        cluster_submit_cmd = config.param(job.name, 'cluster_submit_cmd')
   )
                    
                    cmd += " " + config.param(step.name, 'cluster_submit_cmd_suffix')

                    if config.param(step.name, 'cluster_cmd_produces_job_id'):
                        cmd = job.id + "=$(" + cmd + ")"
                    else:
                        cmd += "\n" + job.id + "=" + job.name

                    # Write job parameters in job list file
                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST"
                    cmd += "\nsleep 1"
                    print(cmd)


class TorqueScheduler(Scheduler):

    def gen_ranges(self, lst):
        s = e = None
        for i in sorted(lst):
            if s is None:
                s = e = i
            elif i == e or i == e + 1:
                e = i
            else:
                yield (s, e)
                s = e = i
    
        if s is not None:
            yield (s, e)

    def submit(self, pipeline):
        self.print_header(pipeline)
        for step in pipeline.step_range:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    #sys.stderr.write("job.job_array_num_task_final: " + job.job_array_num_task_final + "\n")
                    if job.job_array_num_task > 0:
                        #sys.stderr.write("job.job_array_num_task_final:" + job.job_array_num_task_final + "\n")
                        # format num_task output using the gen_ranges() function defined above.
                        job_array_list = list(job.job_array_num_task_final.split(","))
                        job_array_list = map(int, job_array_list)
                        job_array_list_formatted = repr(','.join(['%d' % s if s == e else '%d-%d' % (s, e) for (s, e) in self.gen_ranges(job_array_list)]))
                        cluster_job_array = "--array=" + str(job_array_list_formatted)
                    else:
                        cluster_job_array = ""

                    if job.dependency_jobs:
                        # Chunk JOB_DEPENDENCIES on multiple lines to avoid lines too long
                        max_dependencies_per_line = 50
                        dependency_chunks = [job.dependency_jobs[i:i + max_dependencies_per_line] for i in range(0, len(job.dependency_jobs), max_dependencies_per_line)]
                        job_dependencies = "JOB_DEPENDENCIES=" + ":".join(["$" + dependency_job.id for dependency_job in dependency_chunks[0]])
                        for dependency_chunk in dependency_chunks[1:]:
                            job_dependencies += "\nJOB_DEPENDENCIES=$JOB_DEPENDENCIES:" + ":".join(["$" + dependency_job.id for dependency_job in dependency_chunk])
                    else:
                        job_dependencies = "JOB_DEPENDENCIES="

                    print("""
{separator_line}
# JOB: {job.id}: {job.name}
{separator_line}
JOB_NAME={job.name}
{job_dependencies}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH""".format(
                            job=job,
                            job_dependencies=job_dependencies,
                            separator_line=separator_line
                        )
                    )

                    cmd = \
"""echo "#!/bin/sh
rm -f $JOB_DONE && \\
set -o pipefail
{job.command_with_modules}
MECO_STATE=\$PIPESTATUS
echo MECO_exitStatus:\$MECO_STATE
if [ \$MECO_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MECO_STATE" | \\
""".format(job=job)
                    
                    cluster_submit_cmd     = config.param(job.name, 'cluster_submit_cmd')     
                    cluster_other_arg      = config.param(job.name, 'cluster_other_arg')      
                    cluster_work_dir_arg   = config.param(job.name, 'cluster_work_dir_arg')   
                    cluster_output_dir_arg = config.param(job.name, 'cluster_output_dir_arg') 
                    cluster_job_name_arg   = config.param(job.name, 'cluster_job_name_arg')   
                    cluster_walltime       = config.param(job.subname, 'cluster_walltime')
                    cluster_queue          = config.param(job.subname, 'cluster_queue')
                    cluster_cpu            = config.param(job.subname, 'cluster_cpu')
                    cluster_pmem           = config.param(job.subname, 'cluster_pmem')
                    #cluster_job_array      = config.param(job.subname, 'cluster_job_array')

                    cmd += \
                        cluster_submit_cmd + " " + \
                        cluster_other_arg + " " + \
                        cluster_work_dir_arg + " $OUTPUT_DIR " + \
                        cluster_output_dir_arg + " $JOB_OUTPUT " + \
                        cluster_job_name_arg + " $JOB_NAME " + \
                        cluster_walltime + " " + \
                        cluster_queue + " " + \
                        cluster_pmem + " " + \
                        cluster_cpu + " " + \
                        cluster_job_array
                    
                    if job.dependency_jobs:
                        cmd += " " + config.param(step.name, 'cluster_dependency_arg') + "$JOB_DEPENDENCIES"
                    cmd += " " + config.param(step.name, 'cluster_submit_cmd_suffix')

                    if config.param(step.name, 'cluster_cmd_produces_job_id'):
                        cmd = job.id + "=$(" + cmd + ")"
                    else:
                        cmd += "\n" + job.id + "=" + job.name

                    # Write job parameters in job list file
                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST"

                    print(cmd)

class SlurmScheduler(Scheduler):

    def gen_ranges(self, lst):
        s = e = None
        for i in sorted(lst):
            if s is None:
                s = e = i
            elif i == e or i == e + 1:
                e = i
            else:
                yield (s, e)
                s = e = i
    
        if s is not None:
            yield (s, e)

    def submit(self, pipeline):
        self.print_header(pipeline)
        for step in pipeline.step_range:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    #sys.stderr.write("job.job_array_num_task_final: " + job.job_array_num_task_final + "\n")
                    if job.job_array_num_task > 0:
                        #sys.stderr.write("job.job_array_num_task_final:" + job.job_array_num_task_final + "\n")
                        # format num_task output using the gen_ranges() function defined above.
                        job_array_list = list(job.job_array_num_task_final.split(","))
                        job_array_list = map(int, job_array_list)
                        job_array_list_formatted = repr(','.join(['%d' % s if s == e else '%d-%d' % (s, e) for (s, e) in self.gen_ranges(job_array_list)]))
                        cluster_job_array = "--array=" + str(job_array_list_formatted)
                    else:
                        cluster_job_array = ""

                    if job.dependency_jobs:
                        # Chunk JOB_DEPENDENCIES on multiple lines to avoid lines too long
                        max_dependencies_per_line = 50
                        dependency_chunks = [job.dependency_jobs[i:i + max_dependencies_per_line] for i in range(0, len(job.dependency_jobs), max_dependencies_per_line)]
                        job_dependencies = "JOB_DEPENDENCIES=" + ":".join(["$" + dependency_job.id for dependency_job in dependency_chunks[0]])
                        for dependency_chunk in dependency_chunks[1:]:
                            job_dependencies += "\nJOB_DEPENDENCIES=$JOB_DEPENDENCIES:" + ":".join(["$" + dependency_job.id for dependency_job in dependency_chunk])
                    else:
                        job_dependencies = "JOB_DEPENDENCIES="

                    print("""
{separator_line}
# JOB: {job.id}: {job.name}
{separator_line}
JOB_NAME={job.name}
{job_dependencies}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH""".format(
                            job=job,
                            job_dependencies=job_dependencies,
                            separator_line=separator_line
                        )
                    )

                    cmd = \
"""echo "#!/bin/sh
rm -f $JOB_DONE && \\
set -o pipefail
{job.command_with_modules}
MECO_STATE=\$PIPE_STATUS
echo MECO_exitStatus:\$MECO_STATE
if [ \$MECO_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MECO_STATE" | \\
""".format(job=job)
                    
                    cluster_submit_cmd     = config.param(job.name, 'cluster_submit_cmd')     
                    cluster_other_arg      = config.param(job.name, 'cluster_other_arg')      
                    cluster_work_dir_arg   = config.param(job.name, 'cluster_work_dir_arg')   
                    cluster_output_dir_arg = config.param(job.name, 'cluster_output_dir_arg') 
                    cluster_job_name_arg   = config.param(job.name, 'cluster_job_name_arg')   
                    cluster_walltime       = config.param(job.subname, 'cluster_walltime')
                    cluster_queue          = config.param(job.subname, 'cluster_queue')
                    cluster_cpu            = config.param(job.subname, 'cluster_cpu')
                    cluster_pmem           = config.param(job.subname, 'cluster_pmem')
                    #cluster_job_array      = config.param(job.subname, 'cluster_job_array')

                    cmd += \
                        cluster_submit_cmd + " " + \
                        cluster_other_arg + " " + \
                        cluster_work_dir_arg + " $OUTPUT_DIR " + \
                        cluster_output_dir_arg + " $JOB_OUTPUT " + \
                        cluster_job_name_arg + " $JOB_NAME " + \
                        cluster_walltime + " " + \
                        cluster_queue + " " + \
                        cluster_pmem + " " + \
                        cluster_cpu + " " + \
                        cluster_job_array
                    
                    if job.dependency_jobs:
                        cmd += " " + config.param(step.name, 'cluster_dependency_arg') + "$JOB_DEPENDENCIES"
                    cmd += " " + config.param(step.name, 'cluster_submit_cmd_suffix')

                    if config.param(step.name, 'cluster_cmd_produces_job_id'):
                        cmd = job.id + "=$(" + cmd + ")"
                    else:
                        cmd += "\n" + job.id + "=" + job.name

                    # Write job parameters in job list file
                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST"

                    print(cmd)

class BatchScheduler(Scheduler):
    def submit(self, pipeline):
        self.print_header(pipeline)
        print("SEPARATOR_LINE=`seq -s - 80 | sed 's/[0-9]//g'`")
        for step in pipeline.step_range:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    print("""
{separator_line}
# JOB: {job.name}
{separator_line}
JOB_NAME={job.name}
JOB_DONE={job.done}
printf "\\n$SEPARATOR_LINE\\n"
echo "Begin MECO Job $JOB_NAME at `date +%FT%H.%M.%S`" && \\
rm -f $JOB_DONE && \\
{command_with_modules_and_unload}
MECO_STATE=$PIPESTATUS
echo "End MECO Job $JOB_NAME at `date +%FT%H.%M.%S`"
echo MECO_exitStatus:$MECO_STATE
if [ $MECO_STATE -eq 0 ] ; then touch $JOB_DONE ; else return 0 ; fi""".format(
                            job=job,
                            separator_line=separator_line,
                            #command_with_modules=re.sub(r"\\(.)", r"\1", job.command_with_modules)
                            command_with_modules_and_unload=re.sub(r"\\(.)", r"\1", job.command_with_modules_and_unload)
                        )
                    )
