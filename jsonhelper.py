#!/usr/bin/env python

from __future__ import print_function
import sys
import errno
import json

def parse_tab(fn):
    prefix = None
    samples = {}
    with open(fn) as f:
        for lineno, line in enumerate(f,1):
            fields = line.split()
            assert len(fields) == 5, "expected 5 columns on input, got %d: %s" % (len(fields), fields)
            sample,lib,runid,r1,r2 = fields
            if sample not in samples:
                samples[sample] = {}
            if lib not in samples[sample]:
                samples[sample][lib] = {}

            if runid in samples[sample][lib]:
                print("%s:%d: duplicate entry for %s:%s:%s" % (fn, lineno, sample, lib, runid), file=sys.stderr)
                exit(1)

            samples[sample][lib][runid] = [r1,r2]

            # build a prefix
            if prefix is None:
                prefix = r2
            for i, (c,d) in enumerate(zip(prefix, r1)):
                if c!=d:
                    prefix = prefix[:i]
                    while not prefix.endswith("/"):
                        i-=1
                        prefix = prefix[:i]
                    break

    # trim the prefix
    for sample, libs in samples.iteritems():
        for lib, runs in libs.iteritems():
            for runid, files in runs.iteritems():
                r1 = files[0][len(prefix):]
                r2 = files[1][len(prefix):]
                samples[sample][lib][runid] = [r1,r2]

    d = {"prefix": prefix, "samples": samples}
    return d

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Parse JSON so shell scripts don't have to")
    parser.add_argument("-i", "--import-tsv", action="store_true", help="import tab separated data, output JSON")
    parser.add_argument("input", help="input file")
    parser.add_argument("keys", nargs="*", help="json object keys")

    return parser.parse_args()

def main():
    args = parse_args()

    if args.import_tsv:
        d = parse_tab(args.input)
        json.dump(d, sys.stdout,
                sort_keys=True, indent=4, separators=(',', ': '))
        exit(0)

    with open(args.input) as f:
        d = json.load(f)


    if len(args.keys)>0 and args.keys[0] == "refs":
        for refid, fasta in d["refs"].iteritems():
            print(refid, fasta, sep="\t")
    elif len(args.keys)==0 or args.keys[0] == "samples":
        prefix = d.get("prefix", "")
        if len(prefix) > 0 and prefix[-1] != '/':
            prefix = prefix + '/'

        for sample, libs in d["samples"].iteritems():
            if len(args.keys)>1 and sample != args.keys[1]:
                continue
            for lib, runs in libs.iteritems():
                if len(args.keys)>2 and lib != args.keys[2]:
                    continue
                for runid, files in runs.iteritems():
                    if type(files) != list and len(files) != 2:
                        print("%s:%s:%s: expected two files (R1 and R2)" % (sample, lib, runid), file=sys.stderr)
                        exit(1)
                    r1, r2 = files
                    if r1 == r2:
                        print("%s:%s:%s: has two files with same filename" % (sample, lib, runid), file=sys.stderr)
                        exit(1)

                    print(sample, lib, runid, prefix+r1, prefix+r2, sep="\t")
    elif len(args.keys)>0:
        k = args.keys[0]
        if k in d:
            print(d[k])

if __name__ == "__main__":
    try:
        main()
    except IOError as e:
        if e.errno == errno.EPIPE:
            exit(0)
