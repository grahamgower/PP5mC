#!/usr/bin/env python
# Copyright (c) 2017 Graham Gower <graham.gower@gmail.com>
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
from pipeline import Pipeline
import os, os.path
import errno
import re
import json
import collections

def mkdir_p(path):
    """
    https://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def filecheck(fn):
    if not os.path.exists(fn):
        raise OSError(os.errno.ENOENT, os.strerror(os.errno.ENOENT), fn)
    if not os.access(fn, os.R_OK):
        raise OSError(os.errno.EACCES, os.strerror(os.errno.EACCES), fn)
    return True

def do_fold(p, sample, lib, runid, r1, r2):
    odir = "{}/{}".format(sample, lib)
    pfx = "{}_{}_{}".format(sample, lib, runid)
    opfx = "{}/{}".format(odir, pfx)

    metrics = "{}.metrics".format(opfx)
    fastq = "{}.folded.fastq.gz".format(opfx)

    mkdir_p(odir)

    if os.path.exists(metrics) and os.path.exists(fastq):
        # nothing to do
        return None

    filecheck(r1)
    filecheck(r2)

    #p.new_job("01fold:{}".format(pfx), mem=128, cpus=3, hrs=3)
    j = p.new_job("01fold:{}".format(pfx), mem=128, cpus=3, mins=5)
    j.add_cmd("PATH=/data/acad/programs/pp5mc:$PATH")
    j.add_cmd("module load pigz/2.3.3-foss-2016b")

    fold_cmd = ["foldreads",
            "-m", metrics,
            "-1 {}".format(r1),
            "-2 {}".format(r2),
            "| pigz -p 2 -c -",
            "> {}".format(fastq)]
    j.add_cmd(" \\\n\t".join(fold_cmd))

    return j.sub()

def do_map(p, sample, lib, runid, refid, ref, deps):
    odir = "{}/{}".format(sample, lib)
    pfx = "{}_{}_{}".format(sample, lib, runid)
    opfx = "{}/{}".format(odir, pfx)

    fastq = "{}.folded.fastq.gz".format(opfx)
    bam = "{}.{}.bam".format(opfx, refid)
    bai = "{}.bai".format(bam)

    if os.path.exists(bam) and os.path.exists(bai):
        j = p.new_job("02check:{}".format(pfx), mem=32, mins=5)
        j.add_cmd("module load SAMtools/1.3.1-foss-2016b")
        j.add_cmd("samtools quickcheck {}".format(bam))
        return j.sub()

    j = p.new_job("02map:{}".format(pfx), mem=12*1024, mins=10, cpus=3)
    j.add_cmd("module load BWA/0.7.15-foss-2016b")
    j.add_cmd("module load SAMtools/1.3.1-foss-2016b")

    bwa_cmd = ["bwa mem",
            "-t 1",
            "-R \"@RG\tID:{}\tSM:{}\tLB:{}\tPL:ILLUMINA\"".format(pfx, sample, lib),
            "-C", # Carry through the fastq comments
            ref,
            fastq,
            "|samtools view",
            "-q 25",
            "-Sbu",
            "-",
            "|samtools sort",
            "-O bam",
            "-m 1G",
            "-@ 1",
            "-T tmp.sort.{}.{}".format(pfx, refid),
            "-",
            ">{}".format(bam)]
    j.add_cmd(" \\\n\t".join(bwa_cmd))

    j.add_cmd("samtools index {}".format(bam))

    return j.sub(afterok=deps)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Pre-process hairpin-ligated bisulfit-treated sequencing data")
    parser.add_argument("--dryrun", action="store_true", help="don't submit jobs to the queue")
    parser.add_argument("--email", type=str, help="email address for Slurm to send failure messages")
    parser.add_argument("config_fn", metavar="config.json", help="config file for pipeline run")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    with open(args.config_fn) as f:
        d = json.load(f)

    if "refs" not in d:
        print("Error: {}: no `refs' specified".format(args.config_fn), file=sys.stderr)
        exit(1)

    if "samples" not in d:
        print("Error: {}: no `samples' specified".format(args.config_fn), file=sys.stderr)
        exit(1)

    email = d.get("email", None)
    if args.email:
        email = args.email

    acctdir = "{}/acct".format(os.getcwd())
    mkdir_p(acctdir)

    prefix = d.get("prefix", "")
    if len(prefix) > 0 and prefix[-1] != '/':
        prefix = prefix + '/'

    queue = d.get("queue", "batch")

    p = Pipeline(queue, email, acctdir, args.dryrun)

    map_jobs = collections.defaultdict(list)

    for sample, libs in d["samples"].iteritems():
        for lib, runs in libs.iteritems():
            for runid, files in runs.iteritems():

                #
                # Sanity checks.
                #
                if type(files) != list and len(files) != 2:
                    print("Error: {}:{}:{}: expected two files (R1 and R2)".format(sample, lib, runid), file=sys.stderr)
                    exit(1)

                r1, r2 = [prefix+f for f in files]

                if r1 == r2:
                    print("Error: {}:{}:{}: has two files with same filename".format(sample, lib, runid), file=sys.stderr)
                    exit(1)

                # XXX: relax this requirement?
                if not re.search("R1", r1):
                    print("Error: {}: doesn't have `R1' in the filename, bailing".format(r1))
                    exit(1)
                if not re.search("R2", r2):
                    print("Error: {}: doesn't have `R2' in the filename, bailing".format(r2))
                    exit(1)


                #
                fold_jobid = do_fold(p, sample, lib, runid, r1, r2)

                for refid, ref in d["refs"].iteritems():
                    map_jobid = do_map(p, sample, lib, runid, refid, ref, fold_jobid)
                    map_jobs[(sample,lib,refid)].append(map_jobid)

    print("waiting for...")
    for k, map_jobid in map_jobs.iteritems():
        print(k, map_jobid)
