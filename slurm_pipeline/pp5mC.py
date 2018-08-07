#!/usr/bin/env python
# Copyright (c) 2017,2018 Graham Gower <graham.gower@gmail.com>
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
import os.path
import stat
import errno
import re
import json
import collections
import subprocess

from pipeline import Pipeline
from jobscript import Resource

def mkdir_p(fn):
    try:
        os.makedirs(fn)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(fn):
            pass
        else:
            raise

def rm_f(files):
    if not isinstance(files, list):
        files = [files]
    for fn in files:
        try:
            os.unlink(fn)
        except OSError as e:
            if e.errno != errno.ENOENT:
                raise

def filecheck(fn):
    if not os.path.exists(fn):
        raise OSError(os.errno.ENOENT, os.strerror(os.errno.ENOENT), fn)
    if not os.access(fn, os.R_OK):
        raise OSError(os.errno.EACCES, os.strerror(os.errno.EACCES), fn)
    return True


def do_fold(p, rd, jobsize, sample, lib, runid, r1, r2):
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

    res = Resource(rd["01:fold"], jobsize)
    j = p.new_job("01:fold:{}".format(pfx), res, force=True)

    fold_cmd = ["foldreads",
            "-m", metrics,
            "-1 {}".format(r1),
            "-2 {}".format(r2),
            "| gzip -c -",
            "> {}".format(fastq)]
    j.add_cmd(" \\\n\t".join(fold_cmd))

    return j.sub()

def do_map(p, rd, jobsize, sample, lib, runid, refid, ref, deps):
    odir = "{}/{}".format(sample, lib)
    pfx = "{}_{}_{}".format(sample, lib, runid)
    opfx = "{}/{}".format(odir, pfx)

    fastq = "{}.folded.fastq.gz".format(opfx)
    bam = "{}.{}.bam".format(opfx, refid)
    bai = "{}.bai".format(bam)

    if os.path.exists(bam) and os.path.exists(bai):
        ret = subprocess.call(["samtools", "quickcheck", bam])
        if ret == 0:
            # nothing to do
            return None

    # remove files if they exist, we will remap
    if os.path.lexists(bai):
        os.unlink(bai)
    if os.path.lexists(bam):
        os.unlink(bam)

    res = Resource(rd["02:map"], jobsize)
    j = p.new_job("02:map:{}.{}".format(pfx, refid), res, force=True)

    bwa_cmd = ["bwa mem",
            "-t {}".format(res.cpus),
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

def bam2rgids(bam):
    ids = set()
    hdr = subprocess.check_output(["samtools", "view", "-H", bam])
    for hline in hdr.split("\n"):
        if not hline.startswith("@RG"):
            continue
        for field in hline.split("\t")[1:]:
            tag, val = field.split(":", 1)
            if tag == "ID":
                ids.add(val)
                break
    return ids

def do_merge_runs(p, rd, jobsize, sample, lib, runid_list, refid, deps):
    odir = "{}/{}".format(sample, lib)
    pfx = "{}_{}.{}".format(sample, lib, refid)
    opfx = "{}/{}".format(odir, pfx)
    obam = "{}.bam".format(opfx)
    obai = "{}.bai".format(obam)

    if os.path.exists(obam) and os.path.exists(obai):
        ret = subprocess.call(["samtools", "quickcheck", obam])
        if ret == 0:
            # check header for the expected readgroups
            exp_rgs = {"{}_{}_{}".format(sample,lib,r) for r in runid_list}
            hdr_rgs = bam2rgids(obam)
            if hdr_rgs == exp_rgs:
                # nothing to do
                return None, obam

    if len(runid_list) == 1:
        # nothing to merge, just symlink
        runid = runid_list[0]
        rbam = "{}_{}_{}.{}.bam".format(sample, lib, runid, refid)
        rbai = "{}.bai".format(rbam)
	rm_f([obam, obai])
        os.symlink(rbam, obam)
        os.symlink(rbai, obai)
        # pass through dependencies
        return deps, obam

    # merge list
    ibamlist = []
    for runid in runid_list:
        rpfx = "{}_{}_{}".format(sample, lib, runid)
        rbam = "{}/{}.{}.bam".format(odir, rpfx, refid)
        ibamlist.append(rbam)

    res = Resource(rd["03:mergeruns"], jobsize)
    j = p.new_job("03:mergeruns:{}".format(pfx), res, force=True)

    j.add_cmd("rm -f {}".format(obai))
    j.add_cmd("rm -f {}".format(obam))
    merge_cmd = ["samtools merge", obam]
    merge_cmd.extend(ibamlist)
    j.add_cmd(" \\\n\t".join(merge_cmd))

    j.add_cmd("samtools index {}".format(obam))
    jobid = j.sub(afterok=deps)

    return jobid, obam

def do_dedup(p, rd, jobsize, sample, lib, refid, deps):
    odir = "{}/{}".format(sample, lib)
    pfx = "{}_{}.{}".format(sample, lib, refid)
    opfx = "{}/{}".format(odir, pfx)
    ibam = "{}.bam".format(opfx)
    obam = "{}.dedup.bam".format(opfx)
    obai = "{}.bai".format(obam)

    if os.path.exists(obam) and os.path.exists(obai):
        ret = subprocess.call(["samtools", "quickcheck", obam])
        if ret == 0:
            # nothing to do
            return None, obam

    res = Resource(rd["04:dedup"], jobsize)
    j = p.new_job("04:dedup:{}".format(pfx), res, force=True)

    j.add_cmd("rm -f {}".format(obai))
    j.add_cmd("rm -f {}".format(obam))

    dedup_cmd = ["rmdup_collapsed.py --remove-duplicates",
                "< {}".format(ibam),
                "> {}".format(obam)]
    j.add_cmd(" \\\n\t".join(dedup_cmd))

    j.add_cmd("samtools index {}".format(obam))
    jobid =  j.sub(afterok=deps)

    return jobid, obam


def do_merge_libs(p, rd, jobsize, sample, lib_list, sl_info, refid, deps):
    odir = sample
    pfx = "{}.{}".format(sample, refid)
    opfx = "{}/{}".format(odir, pfx)
    obam = "{}.bam".format(opfx)
    obai = "{}.bai".format(obam)

    if os.path.exists(obam) and os.path.exists(obai):
        ret = subprocess.call(["samtools", "quickcheck", obam])
        if ret == 0:
            # check header for the expected readgroups
            exp_rgs = set()
            for lib in lib_list:
                for runid in sl_info[(sample,lib)]:
                    exp_rgs.add("{}_{}_{}".format(sample, lib, runid))
            hdr_rgs = bam2rgids(obam)
            if hdr_rgs == exp_rgs:
                # nothing to do
                return None, obam

    if len(lib_list) == 1:
        # nothing to merge, just symlink
        lib = lib_list[0]
        lbam = "{}/{}_{}.{}.dedup.bam".format(lib, sample, lib, refid)
        lbai = "{}.bai".format(lbam)
	rm_f([obam, obai])
        os.symlink(lbam, obam)
        os.symlink(lbai, obai)
        # pass through dependencies
        return deps, obam

    # merge list
    ibamlist = []
    for lib in lib_list:
        lbam = "{}/{}/{}_{}.{}.dedup.bam".format(sample, lib, sample, lib, refid)
        ibamlist.append(lbam)

    res = Resource(rd["05:mergelibs"], jobsize)
    j = p.new_job("05:mergelibs:{}".format(pfx), res, force=True)

    j.add_cmd("rm -f {}".format(obai))
    j.add_cmd("rm -f {}".format(obam))

    merge_cmd = ["samtools merge", obam]
    merge_cmd.extend(ibamlist)
    j.add_cmd(" \\\n\t".join(merge_cmd))

    j.add_cmd("samtools index {}".format(obam))

    jobid = j.sub(afterok=deps)

    return jobid, obam

def do_indel_realign(p, rd, jobsize, sample, refid, ref, deps):
    odir = sample
    pfx = "{}.{}".format(sample, refid)
    opfx = "{}/{}".format(odir, pfx)
    ibam = "{}.bam".format(opfx)
    obam1 = "{}.realigned.bam".format(opfx)
    intervals = "{}.intervals".format(opfx)

    obam2 = "{}.{}.bam".format(sample, refid)
    obai2 = "{}.bai".format(obam2)

    if os.path.exists(obam2) and os.path.exists(obai2):
        ret = subprocess.call(["samtools", "quickcheck", obam2])
        if ret == 0:
            # nothing to do
            return None, obam2

    # train realigner
    res1 = Resource(rd["06:realign1"], jobsize)
    j1 = p.new_job("06:realign1:{}".format(pfx), res1, force=True)
    j1.add_cmd("rm -f {}".format(intervals))
    gatk_cmd1 = ["$gatk", # must have this defined somewhere
                "-T RealignerTargetCreator",
                "-nt {}".format(res1.cpus),
                "-R {}".format(ref),
                "-I {}".format(ibam),
                "-o {}".format(intervals)]
    j1.add_cmd(" \\\n\t".join(gatk_cmd1))
    stage1 = j1.sub(afterok=deps)

    # do indel realignment
    res2 = Resource(rd["06:realign2"], jobsize)
    j2 = p.new_job("06:realign2:{}".format(pfx), res2, force=True)
    j2.add_cmd("rm -f {}".format(obam1))
    gatk_cmd2 = ["$gatk",
                "-T IndelRealigner",
                "-R {}".format(ref),
                "-I {}".format(ibam),
                "--bam_compression 0",
                "--disable_bam_indexing",
                "-targetIntervals {}".format(intervals),
                "-o {}".format(obam1)]
    j2.add_cmd(" \\\n\t".join(gatk_cmd2))
    j2.add_cmd("rm {}".format(intervals))
    stage2 = j2.sub(afterok=stage1)

    # regenerate MD tags to match realignments
    res3 = Resource(rd["06:realign3"], jobsize)
    j3 = p.new_job("06:realign3:{}".format(pfx), res2, force=True)
    j3.add_cmd("rm -f {}".format(obam2))
    j3.add_cmd("rm -f {}".format(obai2))
    calmd_cmd = ["samtools calmd",
                "-b",
                obam1,
                ref,
                "> {}".format(obam2)]
    j3.add_cmd(" \\\n\t".join(calmd_cmd))

    j3.add_cmd("samtools index {}".format(obam2))
    stage3 = j3.sub(afterok=stage2)

    return stage3, obam2

def do_mark_5mC(p, rd, jobsize, sample, refid, ref, deps):
    pfx = "{}.{}".format(sample, refid)
    ibam = "{}.bam".format(pfx)
    methlist = "{}.methlist.txt.gz".format(pfx)

    skip = True
    for fmt in ("methylkit", "pileOmeth"):
        for ctx in ("CpG", "CHG", "CHH"):
            if not os.path.exists("{}.{}.{}.txt.gz".format(pfx, fmt, ctx)):
                skip = False
                break
        if skip == False:
            break
    if os.path.exists(methlist) and skip:
        # all files present, nothing to do
        return None

    res = Resource(rd["07:mark5mC"], jobsize)
    j = p.new_job("07:mark5mC:{}".format(pfx), res, force=True)

    j.add_cmd("rm -f {}".format(methlist))
    cmd1 = ["mark5mC",
            "-5 10",
            "-3 10",
            ibam,
            ref,
            "| gzip -c -",
            "> {}".format(methlist)]
    j.add_cmd(" \\\n\t".join(cmd1))

    j.add_cmd("rm -f {}.methylkit.{{CpG,CHG,CHH}}.txt.gz".format(pfx))
    j.add_cmd("rm -f {}.pileOmeth.{{CpG,CHG,CHH}}.txt.gz".format(pfx))
    cmd2 = ["frobmethlist.py",
            "--all",
            "--gzip",
            methlist,
            pfx]
    j.add_cmd(" \\\n\t".join(cmd2))

    return j.sub(afterok=deps)

def do_scanbp(p, rd, jobsize, sample, refid, ref, deps):
    pfx = "{}.{}".format(sample, refid)
    ibam = "{}.bam".format(pfx)
    pairs_txt = "{}.pairs.txt".format(pfx)
    pairs_pdf = "{}.pairs.pdf".format(pfx)

    if os.path.exists(pairs_txt) and os.path.exists(pairs_pdf):
        # all files present, nothing to do
        return None

    res = Resource(rd["08:scanbp"], jobsize)
    j = p.new_job("08:scanbp:{}".format(pfx), res, force=True)

    j.add_cmd("rm -f {} {}".format(pairs_txt, pairs_pdf))
    cmd1 = ["scanbp",
            ibam,
            "> {}".format(pairs_txt)]
    j.add_cmd(" \\\n\t".join(cmd1))

    cmd2 = ["plot_nt_pairing.py",
            "--title {}".format(pfx),
            pairs_txt,
            pairs_pdf]
    j.add_cmd(" \\\n\t".join(cmd2))

    return j.sub(afterok=deps)

def do_samtools_flagstat(p, rd, jobsize, pfx, bam, deps):
    flagstat = "{}.flagstat.txt".format(bam)

    if os.path.exists(flagstat):
        return None

    res = Resource(rd["09:flagstat"], jobsize)
    j = p.new_job("09:flagstat:{}".format(pfx), res, force=True)

    j.add_cmd("samtools flagstat {} > {}".format(bam, flagstat))

    return j.sub(afterok=deps)

def do_samtools_stats(p, rd, jobsize, pfx, bam, deps):
    stats_mt = "{}.stats.MT.txt".format(bam)
    stats_aut = "{}.stats.Aut.txt".format(bam)
    stats_X = "{}.stats.chrX.txt".format(bam)

    if os.path.exists(stats_mt) and os.path.exists(stats_aut) and os.path.exists(stats_X):
        return None

    res = Resource(rd["10:stats"], jobsize)
    j = p.new_job("10:stats:{}".format(pfx), res, force=True)

    chrmax = 50
    autlist = [str(c) for c in range(1,chrmax)] \
                + ["chr{}".format(c) for c in range(1,chrmax)] \
                + ["Chr{}".format(c) for c in range(1,chrmax)]
    autstr = " ".join(autlist)

    # 'samtools stats' doesn't fail if you specify a non-existant chr name,
    # so we try multiple common names in each case.
    j.add_cmd("samtools stats {} MT Mt M > {}".format(bam, stats_mt))
    j.add_cmd("samtools stats {} {} > {}".format(bam, autstr, stats_aut))
    j.add_cmd("samtools stats {} X chrX ChrX > {}".format(bam, stats_X))

    return j.sub(afterok=deps)

def do_samtools_bedcov(p, rd, jobsize, pfx, ref_fa, bam, deps):
    fai = "{}.fai".format(ref_fa)
    bed = "{}.bed.tmp".format(bam)
    bedcov = "{}.bedcov.txt".format(bam)

    if os.path.exists(bedcov):
        return None

    res = Resource(rd["11:bedcov"], jobsize)
    j = p.new_job("11:bedcov:{}".format(pfx), res, force=True)

    j.add_cmd("awk 'BEGIN {OFS=\"\t\"} {print $1, 0, $2}' {} > {}".format(fai, bed))
    j.add_cmd("samtools bedcov {} {} > {}".format(bed, bam, bedcov))
    j.add_cleanup_cmd = "rm -f {}".format(bed)

    return j.sub(afterok=deps)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Pre-process hairpin-ligated bisulfite-treated sequencing data")
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

    if "resources" not in d:
        print("Error: {}: no `resources' specified".format(args.config_fn), file=sys.stderr)
        exit(1)

    resources = ["01:fold", "02:map", "03:mergeruns", "04:dedup", "05:mergelibs",
            "06:realign1", "06:realign2", "06:realign3", "07:mark5mC"]
    for r in resources:
        if r not in d["resources"]:
            print("Error: {}: missing resources for `{}'".format(args.config), file=sys.stderr)
            exit(1)

    if "gatk" not in os.environ:
        print("Error: $gatk not set. E.g first run:\nexport gatk='java -Xmx2g -XX:ParallelGCThreads=2 -jar /path/to/GenomeAnalysisTK.jar'", file=sys.stderr)
        exit(1)

    email = d.get("email", None)
    if args.email:
        email = args.email

    acctdir = "{}/acct".format(os.getcwd())
    mkdir_p(acctdir)

    prefix = d.get("prefix", "")
    if len(prefix) > 0 and prefix[-1] != '/':
        prefix = prefix + '/'

    queue = d.get("queue", None)

    p = Pipeline(queue, email, acctdir, args.dryrun)


    slr_info = dict()
    sl_info = collections.defaultdict(list)
    sl_size = collections.defaultdict(lambda:0)
    s_info = collections.defaultdict(list)
    s_size = collections.defaultdict(lambda:0)

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

                # XXX: R1/R2 might be elsewhere in the filename; foolish requirement?
                if not re.search("R1", r1):
                    print("Error: {}: doesn't have `R1' in the filename, bailing".format(r1))
                    exit(1)
                if not re.search("R2", r2):
                    print("Error: {}: doesn't have `R2' in the filename, bailing".format(r2))
                    exit(1)

                # A measure of the size of the job, in Mb
                filesize = sum((os.stat(f).st_size for f in (r1,r2))) / float(2**20)

                slr_info[(sample, lib, runid)] = (r1, r2, filesize)
                sl_info[(sample, lib)].append(runid)
                sl_size[(sample,lib)] += filesize
                s_size[sample] += filesize

            s_info[sample].append(lib)

    rd = d["resources"]

    # tasks per runid
    map_jobs = collections.defaultdict(list)
    for (sample, lib, runid), (r1,r2,filesize) in slr_info.iteritems():
                fold_jobid = do_fold(p, rd, filesize, sample, lib, runid, r1, r2)

                for refid, ref in d["refs"].iteritems():
                    map_jobid = do_map(p, rd, filesize, sample, lib, runid, refid, ref, fold_jobid)
                    map_jobs[(sample,lib,refid)].append(map_jobid)

    # tasks per library
    dedup_jobs = collections.defaultdict(list)
    for (sample, lib), runid_list in sl_info.iteritems():
        pfx = "{}_{}.{}".format(sample, lib, refid)
        filesize = sl_size[(sample,lib)]
        deplist = map_jobs[(sample,lib,refid)]
        merge_runs_jobid, raw_bam = do_merge_runs(p, rd, filesize, sample, lib, runid_list, refid, deplist)
        flagstat1_jobid = do_samtools_flagstat(p, rd, filesize, pfx+":raw", raw_bam, merge_runs_jobid)
        dedup_jobid, dedup_bam = do_dedup(p, rd, filesize, sample, lib, refid, merge_runs_jobid)
        flagstat2_jobid = do_samtools_flagstat(p, rd, filesize, pfx+":dedup", dedup_bam, dedup_jobid)
        dedup_jobs[sample].append(dedup_jobid)

    # tasks per sample
    for sample, lib_list in s_info.iteritems():
        pfx = "{}.{}".format(sample, refid)
        filesize = s_size[sample]
        deplist = dedup_jobs[sample]
        merge_libs_jobid, _ = do_merge_libs(p, rd, filesize, sample, lib_list, sl_info, refid, deplist)
        realign_jobid, realigned_bam = do_indel_realign(p, rd, filesize, sample, refid, ref, merge_libs_jobid)

        stats_jobid = do_samtools_stats(p, rd, filesize, pfx, realigned_bam, realign_jobid)
        mark5mC_jobid = do_mark_5mC(p, rd, filesize, sample, refid, ref, realign_jobid)
        scanbp_jobid = do_scanbp(p, rd, filesize, sample, refid, ref, realign_jobid)

    p.print_graph()
