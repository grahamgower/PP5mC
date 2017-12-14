#!/usr/bin/env python
# Copyright (c) 2015,2017 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

from __future__ import print_function
import sys
import os
import tempfile
import subprocess
import time
import math
from Cheetah.Template import Template

script_template = """#!/bin/sh

### Job name.
#SBATCH -J $name

#if $queue
### Queue that job is submitted to.
#SBATCH -p $queue
#end if

#if $output
### Redirect stdout/stderr to this file.
#SBATCH -o $output
#end if

#if $email
### Send email to user when job fails.
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email
#end if

### Request resources.
#SBATCH --nodes $nodes
#SBATCH --ntasks $cpus
#SBATCH --time=${days}-${hours}:${mins}:${secs}
#SBATCH --mem=${mem}

### Inherit environment variables, so modules may be preloaded.
#SBATCH --export=ALL

cleanup() {
    ### Empty functions are syntax errors, so use a null command
    ### in case there are no cleanup commands.
    :

    #for $cmd in $cleanup_commands
    $cmd
    #end for
}

die() {
	ret=\$1
	shift
	echo "Failure with exit code \$ret, command was: \$@"
	cleanup
	exit \$ret
}

# make exit status propagate through pipes
set -o pipefail

#for $cmd in $run_commands
$cmd
ret=\$?
[ \$ret = 0 ] || die \$ret "$cmd"

#end for

cleanup
"""

class Resource():
    def __init__(self, rd, scale):
        assert isinstance(rd,dict), "Unexpected resource dict: %s" % type(rd)
        self.scale = scale

        self.nodes = rd.get("nodes", 1)

        if "cpus" in rd:
            self.cpus = int(rd["cpus"])*scale
        elif "_cpus" in rd:
            self.cpus = int(rd["_cpus"])
        else:
            self.cpus = 1

        if "mem" in rd:
            mem = self.siunit(rd["mem"])*scale
        elif "_mem" in rd:
            mem = self.siunit(rd["_mem"])
        else:
            mem = 16*2**20

        if mem < 2**30:
            mem_mb = int(math.ceil(mem/float(2**20)))
            self.mem = str(mem_mb) + "M"
        else:
            mem_gb = int(math.ceil(mem/float(2**30)))
            self.mem = str(mem_gb) + "G"

        t = 0 # time in seconds

        if "secs" in rd:
            t += float(rd["secs"])*scale
        elif "_secs" in rd:
            t += float(rd["_secs"])

        if "mins" in rd:
            t += 60 * float(rd["mins"])*scale
        elif "_mins" in rd:
            t += 60 * float(rd["_mins"])

        if "hours" in rd:
            t += 60*60 * float(rd["hours"])*scale
        elif "_hours" in rd:
            t += 60*60 * float(rd["_hours"])

        # round up
        t += math.ceil(t)

        self.days = int(t // (24*60*60))
        t -= self.days * (24*60*60)
        self.hours = int(t // (60*60))
        t -= self.hours * (60*60)
        self.mins = int(t // 60)
        t -= self.mins * 60
        self.secs = int(t)

    def siunit(self, s):
        # XXX: Does Slurm use 1000*1000 for megabyte, or 1024*1024?
        x = 1
        if s[-1] in "Bb":
            s = s[:-1]
        if s[-1] in "Kk":
            # K is kinda dubious for memory but slurm accepts it
            s = s[:-1]
            x = 2**10
        elif s[-1] in "Mm":
            s = s[:-1]
            x = 2**20
        elif s[-1] in "Gg":
            s = s[:-1]
            x = 2**30
        return float(s)*x

class JobScript:
    """
    A class for creating and submitting Slurm job scripts.

    Each command added to the job script has the exit status checked and the
    job will fail if it is non zero.  Cleanup commands are executed regardless
    of whether the job succeeds or fails.  Jobs may depend upon the completion
    of other jobs.

    See test suite below for usage examples.
    """

    ST_failed = set(["BOOT_FAIL", "CANCELLED", "TIMEOUT", "DEADLINE", "FAILED", "NODE_FAIL", "PREEMPTED", "TIMEOUT"])
    ST_success = set(["COMPLETED"])
    ST_incomplete = set(["CONFIGURING", "COMPLETING", "PENDING", "RUNNING", "RESIZING", "SUSPENDED"])

    def __init__(self, name="mr_jobby", email=None, queue=None,
            nodes=1, cpus=1, mem=16,
            days=0, hours=0, mins=0, secs=0,
            res=None):
        self.filename = None
        self.jobid = None

        self.tmpl = Template(script_template)
        self.tmpl.run_commands = []
        self.tmpl.cleanup_commands = []

        self.tmpl.name = name
        self.tmpl.output = None
        self.tmpl.email = email
        self.tmpl.queue = queue

        if res is not None:
            assert isinstance(res, Resource), "Unexpected resource type: %s" % type(Resource)
            self.tmpl.nodes = res.nodes
            self.tmpl.cpus = res.cpus
            self.tmpl.mem = res.mem
            self.tmpl.days = res.days
            self.tmpl.hours = res.hours
            self.tmpl.mins = res.mins
            self.tmpl.secs = res.secs
        else:
            self.tmpl.nodes = nodes
            self.tmpl.cpus = cpus
            self.tmpl.mem = mem
            self.tmpl.days = days
            self.tmpl.hours = hours
            if days == 0 and hours == 0 and mins == 0:
                    self.tmpl.mins = 5
            else:
                    self.tmpl.mins = mins
            self.tmpl.secs = secs

    def add_cmd(self, cmd):
        """Add a command to be run"""

        self.tmpl.run_commands.append(cmd)

    def add_cleanup_cmd(self, cmd):
        """
        Add a command which will be run at the end the job script.

        This command will execute regardless of whether the job succeeds or fails.
        """

        self.tmpl.cleanup_commands.append(cmd)

    def write(self, filename):
        """Write a job script to the specified file."""

        self.filename = filename
        with open(filename, "w") as f:
            f.write(str(self.tmpl))
    
    def wait(self):
        """Wait for the submitted job to complete."""

        if self.jobid == None:
            raise RuntimeError("Job has not been submitted yet.")

        cmdlist = ["sacct", "-n", "--format=State", "-j", self.jobid]
        while 1:
            statusstr = subprocess.check_output(cmdlist)
            status = statusstr.split("\n")[0].strip()
            if status in self.ST_failed:
                #print("job", self.jobid, "failed", file=sys.stderr)
                return
            if status in self.ST_success:
                return
            time.sleep(1)


    def sub(self, afterok=None):
        """
        Submit job to Slurm using sbatch(1).

        afterok -- job dependencies, a job id or list of job ids after which
                   this job should run.

        Returns: the job id.
        """

        cmdlist = ["sbatch"]

        if afterok:
            if not isinstance(afterok, list):
                afterok = [afterok]
            dep = "--dependency=afterok:{}".format(":".join(afterok))
            cmdlist.append(dep)

        if self.filename:
            cmdlist.append(self.filename)
        else:
            # Use a temp file for the script.
            # Piping to with the correct parameters is more complex.
            fd, tmp_script = tempfile.mkstemp()
            with os.fdopen(fd, "w") as f:
                f.write(str(self.tmpl))
            cmdlist.append(tmp_script)

        try:
            jobstr = subprocess.check_output(cmdlist)
        finally:
            if self.filename is None:
                os.unlink(tmp_script)
                #pass

        self.jobid = jobstr.split()[3]
        #print(tmp_script, "submitted, with jobid", self.jobid, file=sys.stderr)
        return self.jobid


if __name__ == "__main__":
    """Tests for the JobScript class"""

    import os.path

    _, tmp = tempfile.mkstemp(dir=os.environ["HOME"])
    os.unlink(tmp)

    # test multiple commands
    j1 = JobScript()
    j1.add_cmd("echo blah")
    j1.add_cmd("touch {}".format(tmp))
    j1.add_cmd("echo blah2")
    j1.sub()
    j1.wait()

    if not os.path.exists(tmp):
        raise RuntimeError("job 1 failed to create {}".format(tmp))

    os.unlink(tmp)

    # test exit after failed command
    j2 = JobScript()
    j2.add_cmd("/bin/false")
    j2.add_cmd("touch {}".format(tmp))
    j2.sub()
    j2.wait()

    if os.path.exists(tmp):
        os.unlink(tmp)
        raise RuntimeError("job 2 failed to exit after failed command")

    # test cleanup commands
    j3 = JobScript()
    j3.add_cmd("touch {}".format(tmp))
    j3.add_cleanup_cmd("echo rm -f {}".format(tmp))
    j3.add_cleanup_cmd("rm -f {}".format(tmp))
    j3.sub()
    j3.wait()

    if os.path.exists(tmp):
        os.unlink(tmp)
        raise RuntimeError("job 3 failed to remove {}".format(tmp))

    # test cleanup after failed command
    j4 = JobScript()
    j4.add_cmd("touch {}".format(tmp))
    j4.add_cmd("/bin/false")
    j4.add_cleanup_cmd("echo rm -f {}".format(tmp))
    j4.add_cleanup_cmd("rm -f {}".format(tmp))
    j4.sub()
    j4.wait()

    if os.path.exists(tmp):
        os.unlink(tmp)
        raise RuntimeError("job 4 failed to remove {} following failed command".format(tmp))

    # test "afterok" job dependencies
    j5 = JobScript()
    j5.add_cmd("touch {}.5".format(tmp))
    j5id = j5.sub()

    j6 = JobScript()
    j6.add_cmd("test -f {}.5".format(tmp))
    j6.add_cmd("touch {}.6".format(tmp))
    j6id = j6.sub(afterok=j5id)

    j7 = JobScript()
    j7.add_cmd("test -f {}.5".format(tmp))
    j7.add_cmd("test -f {}.6".format(tmp))
    j7.add_cmd("touch {}.7".format(tmp))
    j7.sub([j5id, j6id])
    j7.wait()

    failed = []
    for n in (5, 6, 7):
        filename = "{}.{}".format(tmp, n)
        if os.path.exists(filename):
            os.unlink(filename)
        else:
            failed.append(filename)

    if len(failed) != 0:
        raise RuntimeError("job 5/6/7 dependency problem? Missing file(s): {}".format(failed))


    # test that redirections work as expected
    j8 = JobScript()
    j8.add_cmd("sh -c 'echo foo; /bin/false' > {}".format(tmp))
    j8.add_cmd("echo > {}".format(tmp)) # truncate file, should not execute
    j8.sub()
    j8.wait()

    if not os.path.exists(tmp):
        raise RuntimeError("job 8 failed to create {} with file redirection".format(tmp))

    try:
        with open(tmp) as f:
            s = f.readline()
    finally:
        os.unlink(tmp)

    if s != "foo":
        RuntimeError("job 8 failed to redirect 'foo' into {}".format(tmp))

