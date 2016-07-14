#!/usr/bin/env python
# Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
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
"""
Simulate damaged and/or bisulfite treated DNA fragments, sequenced from
a library with a hairpin on one end.

 p5  5'   (-) strand   3'   /----\
-----|------------------|---      \
-----|------------------|---      /
 p7  3'   (+) strand   5'   \----/

When denatured, the hairpin unfolds:

 p5  5'   (-) strand   3'  hairpin  5'   (+) strand   3'  p7
-----|------------------|-----------|------------------|------

This program simulates the two possible outcomes following adapter removal
and merging of paired end reads (e.g. by AdapterRemoval or leeHom).
1) Reads overlap and are collapsed into a single sequence.
2) Reads do not overlap, remaining as uncollapsed paired sequences.
"""

import sys
import random
import numpy as np
from scipy.stats import geom, randint, poisson

def randstr(size, alpha="ACGT"):
    return "".join(alpha[i] for i in randint.rvs(0, len(alpha), size=size))

complement = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N", "c": "g", "g": "c"}
def revcomp(seq):
    return "".join(complement[s] for s in reversed(seq))

def damage(seq,
        aa, bb,
        sigma_ss, # Cytosine deamination rate for single stranded DNA
        sigma_ds, # Cytosine deamination rate for double stranded DNA
        ):
    """
    Damage a DNA sequence according to Briggs et al. (2007).

    aa       -- Overhanging bases on the left.
    bb       -- Overhanging bases on the right.
    sigma_ss -- Cytosine deamination rate for single stranded DNA
    sigma_ds -- Cytosine deamination rate for double stranded DNA

    Assume 25% chance for each of:
    (1) 5' overhang on both strands.
    (2) 3' overhang on both strands.
    (3) 5' and 3' overhang on (-) strand.
    (4) 5' and 3' overhang on (+) strand.
    """

    x1 = random.random()
    if x1 < 0.25:
        # 5' overhang on both strands.
        # End repair ensures no observed damage in terminal regions.
        l_overhang, r_overhang = aa, bb
        l_chomp, r_chomp = 0, 0
    elif x1 < 0.5:
        # 3' overhang on both strands.
        # Ends removed by T4 DNA polymerase, no overhangs remain.
        l_overhang, r_overhang = 0, 0
        l_chomp, r_chomp = aa, bb

    # 5' and 3' overhang on the same strand.
    # No observed damage on one end, no overhang on the other.
    elif x1 < 0.75:
        l_overhang, r_overhang = aa, 0
        l_chomp, r_chomp = 0, bb
    else:
        l_overhang, r_overhang = 0, bb
        l_chomp, r_chomp = aa, 0

    r_overhang = len(seq)-1 - r_overhang

    seq = seq[l_chomp:len(seq)-r_chomp]

    # (-), (+) sequences
    dseq1, dseq2 = [], []

    for i, n1 in enumerate(seq):
        n2 = complement[n1]
        if n1 in "Cc":
            x2 = random.random()
            if i < l_overhang:
                if x2 < sigma_ss:
                    # deamination at 5' end of (-) strand, (+) strand repaired
                    n1 = "T"
                    n2 = "A"
            elif i <= r_overhang:
                # deamination, without repair on the other strand
                if x2 < sigma_ds:
                    n1 = "T"
        elif n2 in "Cc":
            x2 = random.random()
            if i > r_overhang:
                if x2 < sigma_ss:
                    # deamination at 5' end of (+) strand, (-) strand repaired
                    n1 = "A"
                    n2 = "T"
            elif i >= l_overhang:
                # deamination, without repair on the other strand
                if x2 < sigma_ds:
                    n2 = "T"
        dseq1.append(n1)
        dseq2.append(n2)

    return "".join(dseq1), "".join(reversed(dseq2))

def methylate(seq, cg_rate=0.2):
    """
    Methylate cytosines in CpG dinucleotides.

    Methylated sites are denoted by a lowercase 'c', other strand methylations
    by a lowercase 'g'.
    Assumes both both strands are methylated, or not, in tandem.

    cg_rate is plucked from thin air and does not reflect reality.
    """

    if len(seq) < 2:
        # ignore short sequences
        return seq

    mseq = []
    iseq = iter(seq)
    s1 = next(iseq)

    for s2 in iseq:
        if (s1, s2) == ("C", "G") and random.random() < cg_rate:
            s1 = "c"
            s2 = "g"
        mseq.append(s1)
        s1 = s2

    mseq.append(s1)

    return "".join(mseq)

def bisulfite(r1, r2, bis_rate=0.99):
    """
    Convert unmethylated C's to T's at rate bis_rate.
    """

    bseq1, bseq2 = [], []

    for n1, n2 in zip(r1, r2):
        if n1 == "C":
            if random.random() < bis_rate:
                n1 = "T"
        if n2 == "C":
            if random.random() < bis_rate:
                n2 = "T"
        bseq1.append(n1)
        bseq2.append(n2)

    return "".join(bseq1), "".join(bseq2)

def mutate(seq, p=0.02, alpha="ACGT"):
    """
    Naive mutations, for e.g. sequencing error.
    """
    mseq = []
    for s in seq:
        if random.random() <= p:
            s = random.choice(alpha)
        mseq.append(s)
    return "".join(mseq)

def sample_sequence(ref1, bis_treat, hairpin, fraglen, l_overhang, r_overhang, sigma_ss, sigma_ds):
    """
    Sample a hairpin library sequence.

    bis_treat   -- Is library bisulfite treated?
    hairpin     -- Hairpin sequence.
    fraglen     -- Length of the initial DNA fragment.
    l_overhang  -- Overhanging bases on the left.
    r_overhang  -- Overhanging bases on the right.
    sigma_ss    -- Cytosine deamination rate for single stranded DNA
    sigma_ds    -- Cytosine deamination rate for double stranded DNA
    """

    pos = random.randint(0, len(ref1)-fraglen)
    seq = ref1[pos:pos+fraglen+1]

    if bis_treat:
        seq = methylate(seq)

    minus, plus = damage(seq, l_overhang, r_overhang, sigma_ss, sigma_ds)

    if bis_treat:
        minus, plus = bisulfite(minus, plus)
        # remove lowercase "methylated" sites
        minus = minus.upper()
        plus = plus.upper()

    return mutate(minus + hairpin + plus), len(minus)

def poisson_mix_rvs(shapes, mix, size):
    """
    Produce random variables from mixed poisson distributions.

    shapes  -- List of parameters for the poisson distributions.
    mix     -- List of mixture proportions for the first len(shapes)-1
               distributions.
    size    -- Number of random variables to produce.
    """

    if len(shapes)-1 != len(mix):
        raise Exception("len(shapes)-1 != len(mix)")

    if sum(mix) > 1:
        raise Exception("sum(mix) > 1")

    rvslist = []
    sizes = [size*m for m in mix]
    sizes.append(size-sum(sizes))

    for sh, sz in zip(shapes, sizes):
        rvs = poisson.rvs(sh, size=sz)
        rvslist.append(rvs)

    rvs = np.concatenate(rvslist)
    np.random.shuffle(rvs)
    return rvs

def fastq_write(f, name, seq, qual):
    f.write("@{}\n".format(name))
    f.write(seq)
    f.write("\n+\n")
    f.write(qual)
    f.write("\n")

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Simulate hairpin sequenced reads from a random reference.")
    parser.add_argument("-p", "--hairpin", default="ATCGTTTTTCGAT", help="Hairpin sequence [%(default)s]")
    parser.add_argument("-b", "--bisulfite", default=False, action="store_true", help="Deaminate unmethylated cytosines [%(default)s]")
    parser.add_argument("-n", "--numreads", default=10000, type=int, help="Number of reads [%(default)s]")
    parser.add_argument("-l", "--readlen", default=100, type=int, help="Length of sequenced reads [%(default)s].")
    parser.add_argument("-o", "--oprefix", default="", help="Prefix for fastq output files [%(default)s]")
    parser.add_argument("--mean-sslen", default=1.7, help="Average length of single stranded overhangs [%(default)s]")
    parser.add_argument("--sigma-ss", default=0.68, help="Cytosine deamination rate for single stranded DNA [%(default)s]")
    parser.add_argument("--sigma-ds", default=0.0097, help="Cytosine deamination rate for double stranded DNA [%(default)s]")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    args.hairpin = args.hairpin.upper()

    ref1 = randstr(1000000)

    # trimodal fragment length distribution
    fraglens = poisson_mix_rvs((60, 80, 100), (0.7, 0.05), args.numreads)

    # overhang distribution
    overhangs = geom.rvs(1.0/(1.0+args.mean_sslen), loc=-1, size=args.numreads*2)

    # output filenames
    collapsed = "{}collapsed.fastq".format(args.oprefix)
    uncollapsed1 = "{}uncollapsed_r1.fastq".format(args.oprefix)
    uncollapsed2 = "{}uncollapsed_r2.fastq".format(args.oprefix)

    hlen = len(args.hairpin)

    with open(collapsed, "w") as f_col, \
            open(uncollapsed1, "w") as f_unc1, \
            open(uncollapsed2, "w") as f_unc2:
        for rnum, fraglen in enumerate(fraglens):

            l_overhang, r_overhang = overhangs[2*rnum:2*rnum+2]
            read, fraglen2 = sample_sequence(ref1, args.bisulfite, args.hairpin, fraglen, l_overhang, r_overhang, args.sigma_ss, args.sigma_ds)

            if fraglen2 + hlen < args.readlen:
                # paired reads would be collapsed into a single sequence
                fastq_write(f_col, rnum, read, len(read)*'B')
            else:
                # uncollapsed reads
                r1 = read[:args.readlen]
                r2 = revcomp(read[-args.readlen:])
                qual_str = len(r1) * 'B'
                fastq_write(f_unc1, rnum, r1, qual_str)
                fastq_write(f_unc2, rnum, r2, qual_str)
