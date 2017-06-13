#!/usr/bin/env python

from __future__ import print_function
import sys
from collections import Counter

try:
    range = xrange
except NameError:
    pass

def label2ssqq(label):
    for field in label.split():
        if field.startswith("XF:Z:"):
            field = field[5:]
            break
    else:
        raise Exception("missing XF field for label: `{}'".format(label))
    
    return field.split("|")

def parse_fq(filename):
    """
    fastq parser
    """

    state = 0
    label = None
    qual = None

    if filename.endswith(".gz"):
        import gzip
        xopen = gzip.open
    else:
        xopen = open

    with xopen(filename) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if len(line) == 0:
                continue

            # fastq
            if line[0] == "@":
                if label is not None:
                    yield label, "".join(seq), "".join(qual)
                state = 1
                label = line
                seq = []
                qual = []
                continue
            elif line[0] == "+":
                state = 2
                continue

            if state == 1:
                seq.append(line)
            elif state == 2:
                qual.append(line)

    if label is not None:
        yield label, "".join(seq), "".join(qual)

def mmfind(needle, haystack, mm_allow=2):
    nlen = len(needle)
    for i in range(0, len(haystack)-nlen):
        mm = 0
        for h, n in zip(haystack[i:i+nlen], needle):
            if h != n:
                mm += 1
                if mm > mm_allow:
                    break
        if mm <= mm_allow:
            return i
    return -1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} r1.fq r2.fq".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    fn1 = sys.argv[1]
    fn2 = sys.argv[2]

    hairpin = "ACGCCGGCGGCAAGTGAAGCCGCCGGCGT"
    rhairpin = "ACGCCGGCGGCTTCACTTGCCGCCGGCGT"

    p5 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    p7 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"

    print("""
<!DOCTYPE html>
<html>
<head>
<style>
.bs {color:blue}
.nobs {color:red}
.mismatch {color:purple}
.hairpin {color:orange}
.yadapter {color:cyan}
</style>
</head>
<body>
<h2>
""")
    def f(s, cls):
        if cls is None:
            return s
        else:
            return "<font class={}>{}</font>".format(cls, s)

    fq1 = parse_fq(fn1)
    fq2 = parse_fq(fn2)

    while True:
    #for i in range(1000):
        try:
            _, ss1, qq1 = next(fq1)
            _, ss2, qq2 = next(fq2)
        except StopIteration:
            break

        hpi = mmfind(hairpin,ss1)
        rhpi = mmfind(rhairpin,ss2)
        p5i = mmfind(p5,ss1)
        p7i = mmfind(p7,ss2)
        s1_list = []
        s2_list = []

        print("<p><tt>")
        for i,(s1,s2,q1,q2) in enumerate(zip(ss1,ss2,qq1,qq2)):
            s1 = s1.upper()
            s2 = s2.upper()
            q1 = ord(q1)-33
            q2 = ord(q2)-33
            s1_cls = s2_cls = None
            if s1 == s2:
                if s1 == 'C' or s1 == 'G':
                    s1_cls = s2_cls = "nobs"
            elif (s1,s2) in [('T','C'), ('G','A')]:
                s1_cls = s2_cls = "bs"
            else:
                if q1<=q2:
                    s1_cls = "mismatch"
                if q1>=q2:
                    s2_cls = "mismatch"
            if hpi != -1 and i>=hpi and i<hpi+len(hairpin):
                s1_cls = "hairpin"
            if rhpi != -1 and i>=rhpi and i<rhpi+len(rhairpin):
                s2_cls = "hairpin"
            if p5i != -1 and i>=p5i:
                s1_cls = "yadapter"
            if p7i != -1 and i>=p7i:
                s2_cls = "yadapter"
            s1_list.append(f(s1, s1_cls))
            s2_list.append(f(s2, s2_cls))
        print("".join(s1_list))
        print("</br>")
        print("".join(s2_list))
        #print("</br>")
        #print(seq)
        print("</tt></p>")

    print("</h2></html></body>")
