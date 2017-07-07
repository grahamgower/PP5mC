#!/usr/bin/env python

from __future__ import print_function
from gzopen import gzopen
import sys

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
    # coordinate is 1-based
    chrbase = "{}.{}".format(chrom, start+1)
    strand = "FR"[strand == '-']
    print(chrbase, chrom, start+1, strand, C+mC,
                "{:.2f}".format(100*mC/(C+mC)),
                "{:.2f}".format(100*C/(C+mC)),
                file=f, sep="\t")

def print_pileOmeth(f, chrom, start, end, C, mC):
    # coordinates are 0-based, half open
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

        for f in mk_files:
            if f is not None:
                print("chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT", file=f)

    if args.pileOmeth:
        if args.cpg:
            pm_files[0] = gzopen("{}.pileOmeth.CpG.txt{}".format(args.oprefix, suffix), "w")
            print("track type=\"bedGraph\" description=\"CpG methylation levels\"", file=pm_files[0])
        if args.chg:
            pm_files[1] = gzopen("{}.pileOmeth.CHG.txt{}".format(args.oprefix, suffix), "w")
            print("track type=\"bedGraph\" description=\"CHG methylation levels\"", file=pm_files[1])
        if args.chh:
            pm_files[2] = gzopen("{}.pileOmeth.CHH.txt{}".format(args.oprefix, suffix), "w")
            print("track type=\"bedGraph\" description=\"CHH methylation levels\"", file=pm_files[2])


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
