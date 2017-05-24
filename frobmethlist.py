#!/usr/bin/env python

from __future__ import print_function
import sys

class gzopen:
    """
    Pipe file through a gzip subprocess.

    This is substantially faster than just using the gzip library,
    particularly when dealing with multiple files at a time.
    """
    def __init__(self, fn, mode="r"):
        from subprocess import Popen, PIPE
        import os

        try:
            if "r" in mode:
                self.pipe = Popen(["gzip", "-dc", fn], bufsize=-1, stdout=PIPE)
                self.f = self.pipe.stdout
            elif "w" in mode:
                self.pipe = Popen(["gzip", "-c"], bufsize=-1, stdin=PIPE, stdout=open(fn, mode))
                self.f = self.pipe.stdin
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                # no gzip
                # ...
                raise
            else:
                raise

        self.read = self.f.read
        self.write = self.f.write

    def __enter__(self):
        return self.f

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self.pipe.stdin is not None:
            self.pipe.stdin.close()
        if self.pipe.stdout is not None:
            self.pipe.stdout.close()


def parse_methlist(filename):
    with gzopen(filename) as f:
        next(f) # skip header
        # rintf("chrom\tpos-0\tpos-1\tstrand\tdepth\tC\tmC\tcontext\n");
        for line in f:
            line = line.rstrip()
            fields = line.split("\t")
            chrom = fields[0]
            start, end = map(int, fields[1:3])
            strand = fields[3]
            depth, C, mC = map(int, fields[4:7])
            context = fields[7]
            yield chrom, start, end, strand, depth, C, mC, context

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Convert mark_5mC output to methylkit/pileOmeth formats")
    parser.add_argument("-m", "--methylkit", action="store_true", default=False)
    parser.add_argument("-p", "--pileOmeth", action="store_true", default=False)
    parser.add_argument("--cpg", action="store_true", default=False)
    parser.add_argument("--chg", action="store_true", default=False)
    parser.add_argument("--chh", action="store_true", default=False)
    parser.add_argument("--all", action="store_true", default=False)
    parser.add_argument("--gzip", action="store_true", default=False)
    parser.add_argument("infile")
    parser.add_argument("oprefix")
    args = parser.parse_args()

    if args.all:
        args.methylkit = args.pileOmeth = True
        args.cpg = args.chg = args.chh = True

    if not args.methylkit and not args.pileOmeth:
        print("Must set one or more of --methylkit or --pileOmeth",
                file=sys.stderr)
        exit(1)

    if not args.cpg and not args.chg and not args.chh:
        print("Must set one or more of --cpg, --chg, or --chh",
                file=sys.stderr)
        exit(1)

    return args

def print_methylkit(f, chrom, start, end, strand, C, mC):
    chrbase = "{}.{}".format(chrom, start)
    strand = "FR"[strand == '-']
    print(chrbase, chrom, start, strand, C+mC, C, mC, file=f, sep="\t")

def print_pileOmeth(f, chrom, start, end, C, mC):
    print(chrom, start, end, int(100*mC/(C+mC)), mC, C, file=f, sep="\t")

if __name__ == "__main__":
    args = parse_args()

    mk_files = [None, None, None]
    pm_files = [None, None, None]

    if args.gzip:
        suffix = ".gz"
    else:
        suffix = ""

    if args.methylkit:
        if args.cpg:
            mk_files[0] = gzopen("{}.methylkit.CpG.txt{}".format(args.oprefix, suffix), "w")
        if args.chg:
            mk_files[1] = gzopen("{}.methylkit.CHG.txt{}".format(args.oprefix, suffix), "w")
        if args.chh:
            mk_files[2] = gzopen("{}.methylkit.CHH.txt{}".format(args.oprefix, suffix), "w")

    if args.pileOmeth:
        if args.cpg:
            pm_files[0] = gzopen("{}.pileOmeth.CpG.txt{}".format(args.oprefix, suffix), "w")
        if args.chg:
            pm_files[1] = gzopen("{}.pileOmeth.CHG.txt{}".format(args.oprefix, suffix), "w")
        if args.chh:
            pm_files[2] = gzopen("{}.pileOmeth.CHH.txt{}".format(args.oprefix, suffix), "w")


    try:
        for line in parse_methlist(args.infile):
            chrom, start, end, strand, depth, C, mC, context = line

            if C+mC == 0 or context[0] != 'C':
                continue

            if context[1] == 'G':
                # CpG
                ctx = 0
            elif context[1] in "ACT":
                if context[2] == 'G':
                    # CHG
                    ctx = 1
                elif context[2] in "ACT":
                    # CHH
                    ctx = 2

            if args.methylkit and mk_files[ctx] is not None:
                print_methylkit(mk_files[ctx], chrom, start, end, strand, C, mC)
            if args.pileOmeth and pm_files[ctx] is not None:
                print_pileOmeth(pm_files[ctx], chrom, start, end, C, mC)

    finally:
        for f in mk_files:
            if f is not None:
                f.close()
        for f in pm_files:
            if f is not None:
                f.close()
