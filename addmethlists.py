#!/usr/bin/env python

from __future__ import print_function
from gzopen import gzopen
import sys

def parse_tsv(fn, istart, iend):
    with gzopen(fn) as f:
        for line in f:
            line = line.rstrip()
            fields = line.split("\t")
            yield fields[istart:iend], fields

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Sum methylation counts from two files")
    parser.add_argument("-m", "--methylkit", action="store_true", default=False)
    parser.add_argument("-p", "--pileOmeth", action="store_true", default=False)
    parser.add_argument("-z", "--gzip", action="store_true", default=False, help="output is gzipped")
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
        istart = 1
        iend = 3
    elif args.pileOmeth:
        istart = 0
        iend = 2
    else:
        raise Exception("Are we doing methylkit or pileOmeth?")

    f1 = parse_tsv(args.in1, istart, iend)
    f2 = parse_tsv(args.in2, istart, iend)

    k1, l1 = next(f1)
    k2, l2 = next(f2)
    oline = None

    with gzopen("/dev/stdout", "w", args.gzip) as fout:
        try:
            while True:
                if k1<k2:
                    oline = l1
                    k1, l1 = next(f1)
                if k2<k1:
                    oline = l2
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
                    k1, l1 = next(f1)
                    k2, l2 = next(f2)

                print(*oline, file=fout, sep="\t")

        except StopIteration:
            if oline == l1:
                oline = l2
                remainder = f2
            else:
                oline = l1
                remainder = f1

            print(*oline, file=fout, sep="\t")
            for _, oline in remainder:
                print(*oline, file=fout, sep="\t")
