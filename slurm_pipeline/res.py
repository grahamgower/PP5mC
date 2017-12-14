#!/usr/bin/env python

from __future__ import print_function
import sys
import math
import json
import os
import collections

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
    if len(sys.argv) < 3:
        print("usage: {} config.json s1.out [... sN.out]".format(sys.argv[0]),
                file=sys.stderr)
        exit(1)

    config_fn = sys.argv[1]
    slurm_files = sys.argv[2:]

    size, refs = parse_config(config_fn)

    rmax = dict()

    for fn in slurm_files:
        slurm_out = parse_slurm_out(fn)
        if slurm_out is None:
            continue
        jobname, slr, mem, cpus, secs = slurm_out
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
