#!/usr/bin/env python

from __future__ import print_function
import sys
import math
import json
import os
import collections
import subprocess

def parse_slurm_out(fn):
    completed = False
    with open(fn) as f:
        for line in f:
            if line.startswith("Job Name"):
                jobname_str = line.split()[3]
            elif line.startswith("State"):
                if line.split()[2] == "COMPLETED":
                    completed = True
            elif line.startswith("Walltime elapsed"):
                walltime_str = line.split()[3]
            elif line.startswith("% CPU used (Total)"):
                cpu_str = line.split()[5].strip("%")
            elif line.startswith("% Mem used (Max)"):
                mem_str = line.split()[6].strip("()")[:-5]

    if not completed:
        return None

    h, m, s = map(int, walltime_str.split(":"))
    secs = h*60*60 + m*60 + s
    
    x = 1
    if mem_str[-1] == "K":
        x = 1024
    elif mem_str[-1] == "M":
        x = 1024*1024
    elif mem_str[-1] == "G":
        x = 1024*1024*1024
    mem_f = x*float(mem_str[:-1])

    mem_mb = int(math.ceil(mem_f/float(2**20)))
    mem = mem_mb #str(mem_mb) + "M"

    cpus = int(math.ceil(float(cpu_str)/100))

    jfields = jobname_str.split(":")
    jobname = ":".join(jfields[:2])
    slr = ":".join(jfields[2:])

    return jobname, slr, mem, cpus, secs


def sacct_job_info(jobid_list):
    cmd = [ "sacct",
            "-n",
            "--format=JobID%-50,JobName%-100,NCPUS,Elapsed,MaxRSS,State",
            "--jobs={}".format(",".join(map(str, jobid_list)))
            ]
    out = subprocess.check_output(cmd)
    for line in out.split("\n"):
        fields = line.split()
        if len(fields) == 0:
            continue

        # ughhh... sacct output is nasty
        if len(fields) == 5:
            jobid = fields[0]
            jobname = fields[1]
            #cpus = int(fields[2])
            #elapsed = fields[3]
            #state = fields[4]
            continue
        elif len(fields) == 6:
            jobid2 = fields[0]
            jobid2 = jobid2[:-len(".batch")]
            if jobid != jobid2:
                continue
            #jobname = fields[1]
            cpus = int(fields[2])
            elapsed = fields[3]
            maxrss = fields[4]
            state = fields[5]
        else:
            raise Exception("sacct gave unexpected number of fields: {}".format(len(fields)))

        if state != "COMPLETED":
            continue

        jfields = jobname.split(":", 2)
        job = ":".join(jfields[:2])
        slr = jfields[2]

        x = 1
        if maxrss[-1] == "K":
            x = 1024
        elif maxrss[-1] == "M":
            x = 1024*1024
        elif maxrss[-1] == "G":
            x = 1024*1024*1024
        mem_f = x*float(maxrss[:-1])

        mem_mb = int(math.ceil(mem_f/float(2**20)))
        mem = mem_mb #str(mem_mb) + "M"

        # sacct(1) says this has format [DD-[HH:]]MM:SS
        # But I've never seen this in the output, its always HH:MM:SS.
        if "-" in elapsed:
            days, elapsed2 = elapsed.split("-", 1)
            days = int(days)
        else:
            elapsed2 = elapsed
            days = 0

        efields = map(int, elapsed2.split(":"))
        if len(efields) == 2:
            h = 0
            m, s = efields
        elif len(efields) == 3:
            h, m, s = efields
        else:
            raise Exception("sacct gave weird elapsed time `{}' for {}".format(elapsed, jobname))

        secs = days*24*60*60 + h*60*60 + m*60 + s

        yield job, slr, mem, cpus, secs

def parse_dependency_graph(fn):
    """
    Extract nodes from the dependency graph.

    NOTE: This is embarrassingly naiive and will only work for our dot file
          (and will easily break if python-graphviz changes things).
    """
    nodes = {}
    edges = []

    with open(fn) as f:
        for line in f:
            if line.startswith("digraph {") or line.startswith("}"):
                continue
            if line.startswith("\t\t"):
                # edge
                n1, n2 = map(int, line.split("->"))
                edges.append((n1,n2))
                if n1 not in nodes:
                    nodes[n1] = None
                if n2 not in nodes:
                    nodes[n2] = None
            elif line.startswith("\t"):
                # node
                node, rest = line.split(None, 1)
                node = int(node)
                _, jobname = rest.split("\"", 1)
                jobname = jobname.strip()
                nodes[node] = jobname

    return nodes, edges

def parse_config(fn):
    with open(fn) as f:
        d = json.load(f)

    prefix = d.get("prefix", "")
    if len(prefix) > 0 and prefix[-1] != '/':
        prefix = prefix + '/'

    size = collections.Counter()

    for sample, libs in d["samples"].iteritems():
        for lib, runs in libs.iteritems():
            for runid, files in runs.iteritems():
                r1, r2 = [prefix+f for f in files]
                filesize = sum((os.stat(f).st_size for f in (r1,r2))) / float(2**20)
                size["_".join((sample, lib, runid))] = filesize
                size["_".join((sample,lib))] += filesize
                size[sample] += filesize

    refs = d["refs"]

    return size, refs


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} config.json dependency_graph".format(sys.argv[0]),
                file=sys.stderr)
        exit(1)

    config_fn = sys.argv[1]
    dg_fn= sys.argv[2]

    size, refs = parse_config(config_fn)
    nodes, _ = parse_dependency_graph(dg_fn)
    jobid_list = nodes.keys()
    if max(jobid_list) == len(jobid_list):
        print("{} is from a --dryrun, cannot obtain slurm accounting info!".format(dg_fn), file=sys.stderr)
        exit(1)
    sacct = sacct_job_info(jobid_list)

    rmax = dict()

    for jobname, slr, mem, cpus, secs in sacct:
        for ref in refs.iterkeys():
            if slr.endswith(ref):
                slr = slr[:-(len(ref)+1)]
                break
        else:
            ref = None
        filesize = size[slr]
        #print(jobname, slr, ref, mem, cpus, secs, filesize)

        if jobname not in rmax:
            rmax[jobname] = (mem, cpus, secs, filesize)
        else:
            m, c, s, f = rmax[jobname]
            m = max(mem, m)
            c = max(cpus, c)
            if secs/filesize > s/f:
                s = secs
                f = filesize
            rmax[jobname] = (m, c, s, f)

    rd = dict()
    b = 1.2 # buffer of 20%
    for jobname in sorted(rmax.iterkeys()):
        m,c,s,f = rmax[jobname]
        m *= b
        s = b * s/f
        rd[jobname] = {"_mem": "{:.2f}".format(m)+"M",
                        "_cpus": c,
                        "secs": "{:.2f}".format(s)}

    json.dump({"resources": rd}, sys.stdout,
            sort_keys=True, indent=4, separators=(',', ': '))
