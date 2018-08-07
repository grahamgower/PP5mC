#!/usr/bin/env python
#
# Measure the fidelity of sequence reconstruction.
#
# Copyright (c) 2018 Graham Gower <graham.gower@gmail.com>
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

def parse_fq(filename):
    """
    Fastq parser. Doesn't do fasta.
    """

    label = None
    comment = None
    qual = None

    if filename.endswith(".gz"):
        import gzip
        xopen = gzip.open
    else:
        xopen = open

    with xopen(filename) as f:
        for lineno, line in enumerate(f):
            line = line.rstrip("\r\n")
            if len(line) == 0:
                continue

            state = lineno%4

            if state == 0:
                assert line[0] == "@", "invalid fastq file"

                if label is not None:
                    yield label, comment, seq, qual
                lfields = line.split(None, 1)
                label = lfields[0]
                if len(lfields) > 1:
                    comment = lfields[1]
                else:
                    comment = None

            elif state == 1:
                seq = line
            elif state == 3:
                qual = line

    if label is not None:
        yield label, comment, seq, qual

def revcomp(seq):
    """
    Reverse complement of sequence @seq.
    """
    revmap = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
    return "".join((revmap[s] for s in reversed(seq)))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} sim.fq recover.fq".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    readlen = 150

    sims = {}
    for label, comment, _, _ in parse_fq(sys.argv[1]):
        tags = {}
        for f in comment.split():
            tag, type, val = f.split(":")
            if type == "i":
                val = int(val)
            tags[tag] = val

        seq = tags["om"]
        #pos = tags["ps"]
        sims[label] = seq.upper()

    nbases = 0
    matches = 0
    for label, _, seq, _ in parse_fq(sys.argv[2]):
        s1 = sims[label][:readlen]
        s2 = seq.upper()
        nbases += len(s1)

        if s1 == s2:
            matches += len(s1)
            continue

        rs1 = revcomp(sims[label][-readlen:])
        if rs1 == s2:
            matches += len(s1)
            continue
        else:
            m = 0 # fwd matches
            rm = 0 # rev matches
            for ss1,rss1,ss2 in zip(s1,rs1,s2):
                if ss1 == ss2:
                    m += 1
                if rss1 == ss2:
                    rm += 1
            matches += max(m,rm)

    if nbases > 0:
        ratio = float(matches)/nbases
    else:
        ratio = 0
    print(matches, nbases, ratio)
