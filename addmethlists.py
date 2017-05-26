#!/usr/bin/env python

from __future__ import print_function
from gzopen import gzopen
import sys
import itertools

def parse_col0(fn):
    with gzopen(fn) as f:
        for line in f:
            yield line.split(None, 1)[0]

def parse_tsv(fn, ichrom, ipos, chrmap, skip=1):
    with gzopen(fn) as f:
        while skip:
            skip -= 1
            next(f)
        for line in f:
            fields = line.split()
            key = (chrmap[fields[ichrom]], int(fields[ipos]))
            yield key, fields

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Sum methylation counts from two files")
    parser.add_argument("-m", "--methylkit", action="store_true", default=False)
    parser.add_argument("-p", "--pileOmeth", action="store_true", default=False)
    parser.add_argument("-z", "--gzip", action="store_true", default=False, help="output is gzipped")
    parser.add_argument("-l", "--chr-list", required=True, help="chromosome list for sort order (e.g. faidx file)")
    parser.add_argument("in1")
    parser.add_argument("in2")
    args = parser.parse_args()

    if not args.methylkit and not args.pileOmeth:
        print("Must set one of --methylkit or --pileOmeth",
                file=sys.stderr)
        exit(1)

    if args.methylkit and args.pileOmeth:
        print("Parameters --methylkit and --pileOmeth are mutually exclusive",
                file=sys.stderr)
        exit(1)

    return args

if __name__ == "__main__":
    args = parse_args()

    if args.methylkit:
        ichrom = 1
        ipos = 2
    elif args.pileOmeth:
        ichrom = 0
        ipos = 1
    else:
        raise Exception("Are we doing methylkit or pileOmeth?")

    chrmap = {s:i for i,s in enumerate(parse_col0(args.chr_list))}
    f1 = parse_tsv(args.in1, ichrom, ipos, chrmap)
    f2 = parse_tsv(args.in2, ichrom, ipos, chrmap)

    l1 = l2 = oline = None


    with gzopen("/dev/stdout", "w", args.gzip) as fout:
        try:
            k1, l1 = next(f1)
            k2, l2 = next(f2)
            while True:
                if k1<k2:
                    oline = l1
                    l1 = None
                    k1, l1 = next(f1)
                elif k2<k1:
                    oline = l2
                    l2 = None
                    k2, l2 = next(f2)
                else:
                    if args.methylkit:
                        oline = l1[:4]
                        C = int(l1[5]) + int(l2[5])
                        mC = int(l1[6]) + int(l2[6])
                        oline.extend([C+mC, C, mC])
                    elif args.pileOmeth:
                        oline = l1[:3]
                        mC = int(l1[4]) + int(l2[4])
                        C = int(l1[5]) + int(l2[5])
                        oline.extend([100*mC/(C+mC), mC, C])

                    l1 = l2 = None
                    k1, l1 = next(f1)
                    k2, l2 = next(f2)

                print(*oline, file=fout, sep="\t")

        except StopIteration:
            pass

        if oline is not None:
            print(*oline, file=fout, sep="\t")
        if l1 is not None:
            print(*l1, file=fout, sep="\t")
        if l2 is not None:
            print(*l2, file=fout, sep="\t")

        for _, lx in itertools.chain(f1, f2):
            print(*lx, file=fout, sep="\t")
