#!/usr/bin/env python

from __future__ import print_function
import sys
from collections import Counter

try:
    range = xrange
except NameError:
    pass

def comment2ssqq(label):
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
    comment = None
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
                    yield label, comment, "".join(seq), "".join(qual)
                state = 1
                lfields = line.split(None, 1)
                label = lfields[0]
                if len(lfields) > 1:
                    comment = lfields[1]
                else:
                    comment = None
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
        yield label, comment, "".join(seq), "".join(qual)

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

def revcomp(seq):
    """
    Reverse complement of sequence @seq.
    """
    revmap = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
    return "".join((revmap[s] for s in reversed(seq)))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: {} folded.fq [hairpin]".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    folded_fn = sys.argv[1]

    if len(sys.argv) == 3:
        hairpin = sys.argv[2]
    else:
        hairpin = "ACGCCGGCGGCAAGTGAAGCCGCCGGCGT"

#    p5 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
#    p7 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
    p5 = p7 = "AGATCGGAAGAGC"

    rhairpin = revcomp(hairpin)

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

    def g(slist):
        sl = []
        out = []
        last = None
        for s, cls in slist:
            if cls != last and len(sl) != 0:
                out.append(f("".join(sl), last))
                last = cls
                sl = []
            sl.append(s)
        if len(sl) != 0:
            out.append(f("".join(sl), last))
        return "".join(out)

    for seqno, (label, comment, ss, qq) in enumerate(parse_fq(folded_fn), 1):

        if seqno > 1000:
            break

        ss1, ss2, qq1, qq2, hlen_str = comment2ssqq(comment)
        hlen = int(hlen_str)

        hpi = rhpi = len(ss)
        p5i = p7i = 2*len(ss)+hlen

        if p5i > len(ss1):
            p5i = mmfind(p5,ss1)
        if p7i > len(ss1):
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
            if hpi != -1 and i>=hpi and i<hpi+hlen:
                s1_cls = "hairpin"
            if rhpi != -1 and i>=rhpi and i<rhpi+hlen:
                s2_cls = "hairpin"
            if p5i != -1 and i>=p5i:
                s1_cls = "yadapter"
            if p7i != -1 and i>=p7i:
                s2_cls = "yadapter"
            s1_list.append((s1, s1_cls))
            s2_list.append((s2, s2_cls))
        print("R1", g(s1_list))
        print("</br>")
        print("R2", g(s2_list))
        print("</br>")
        print("FS ", ss)
        print("</br>")
        print("FQ ", qq)
        print("</tt></p>")

    print("</h2></html></body>")
