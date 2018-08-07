#!/usr/bin/env python
#
# Measure the proportion of reads that map to/near the correct location.
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
import collections

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
                    yield label[1:], comment, seq, qual
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
        yield label[1:], comment, seq, qual

# SAM flag field
F_PAIRED        = 0x001 # the read is paired in sequencing
F_PAIR_MAPPED   = 0x002 # the read is mapped in a proper pair
F_UNMAPPED      = 0x004 # the query sequence itself is unmapped
F_MATE_UNMAPPED = 0x008 # the mate is unmapped
F_STRAND        = 0x010 # strand of the query (1 for reverse)
F_MATE_STRAND   = 0x020 # strand of the mate
F_FIRST_READ    = 0x040 # the read is the first read in a pair
F_SECOND_READ   = 0x080 # the read is the second read in a pair
F_SECONDARY     = 0x100 # the alignment is not primary
F_QCFAIL        = 0x200 # QC failure
F_DUP           = 0x400 # optical or PCR duplicate
F_SUPP          = 0x800 # supplementary alignment

def parse_sam(fn):
    f_filter = F_UNMAPPED|F_SECONDARY|F_QCFAIL|F_DUP|F_SUPP
    with open(fn) as f:
        for line in f:
            if line[0] == "@":
                continue
            line = line.rstrip()
            fields = line.split("\t")
            name = fields[0]
            flag = int(fields[1])
            pos = int(fields[3])
            maq = int(fields[4])
            seq = fields[9]

            if flag & f_filter != 0:
                continue

            if maq < 10:
                continue

            rev = (flag & F_STRAND) == F_STRAND

            yield name, pos, seq, rev

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} sim.fq recover.sam".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    sims = {}
    for label, comment, _, _ in parse_fq(sys.argv[1]):
        tags = {}
        for f in comment.split():
            tag, type, val = f.split(":")
            if type == "i":
                val = int(val)
            tags[tag] = val

        #seq = tags["om"]
        pos = tags["ps"]
        sims[label] = pos

    nseqs = 0
    correctpos = 0
    closepos = 0
    disthist = collections.Counter()

    for label, pos2, seq, rev in parse_sam(sys.argv[2]):
        #if not rev:
        #    continue

        if label.endswith("#/1"):
            label = label[:-3]
        elif label.endswith("#/2"):
            continue

        nseqs += 1
        pos1 = sims[label]

        dist = pos2-pos1
        if dist == 0:
            correctpos += 1
        if abs(dist) < 10:
            closepos += 1
        disthist[dist] += 1

    if nseqs > 0:
        correct = float(correctpos) / nseqs
        close = float(closepos) / nseqs
    else:
        correct = 0
        close = 0

    print(nseqs, correctpos, closepos, correct, close)
    #print(disthist.most_common(20))
