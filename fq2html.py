#!/usr/bin/env python

from __future__ import print_function
import sys
import math
from collections import Counter

try:
    range = xrange
except NameError:
    pass

def comment2ssqq(label):
    tags = {}
    for field in label.split():
        try:
            tag, type, val = field.split(":", 2)
        except ValueError:
            raise Exception("unknown comment field ``{}''".format(field))

        tags[tag] = val

    return tags["r1"], tags["r2"], tags["q1"], tags["q2"], tags.get("hp", "")

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

# Inverse poisson CDF, stolen from bwa: bwtaln.c
def bwa_cal_maxdiff(l, err=0.02, thres=0.01, maxlen=1000):
	esum = elambda = math.exp(-l * err)
	y = 1.0
	k = x = 1
        while (k < maxlen):
		y *= l * err
		x *= k
		esum += elambda * y / x
                if (1.0 - esum < thres):
                    return k
                k += 1
	return 2;

def mmfind(needle, haystack, mm_allow=1, start=0, min_matches=-1):
    nlen = len(needle)
    if min_matches == -1:
        min_matches = nlen
    for i in range(start, len(haystack)-min_matches):
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

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="HTMLify fastq files to troubleshoot hairpin-ligated bisulfite-treated files")
    parser.add_argument("-n", "--nseqs", type=int, default=1000, help="only output this many sequences [(%default)s]")
    parser.add_argument("-p", "--hairpin", action="append", help="hairpin sequence(s)")
    parser.add_argument("--latex", action="store_true", default=False, help="LaTeX output [%(default)s]")
    parser.add_argument("-i", "--interleaved", action="store_true", default=False, help="R1/R2 are interleaved in one fastq")
    parser.add_argument("-m", "--adapter-matchlen", type=int, default=9, help="min number of bases to match hairpin/adapter at end of read")
    parser.add_argument("fq1", metavar="f1.fq", help="fastq r1 or folded")
    parser.add_argument("fq2", metavar="f2.fq", help="fastq r2", nargs='?', default=None)
    args = parser.parse_args()

    if args.hairpin is None:
        args.hairpin = ["ACGCCGGCGGCAAGTGAAGCCGCCGGCGT",
                        "ACGCCGGCGGCAAGTAAGCCGCCGGCGT",
                        "ACGCCGGCGGCAAGTAGCCGCCGGCGT"]

    args.rhairpin = [revcomp(h) for h in args.hairpin]

    return args

if __name__ == "__main__":
    args = parse_args()

    p5 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA" #[:args.adapter_matchlen]
    p7 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" #[:args.adapter_matchlen]

    yy_maxdiff = bwa_cal_maxdiff(len(p5))
    hp_maxdiff = [bwa_cal_maxdiff(len(s)) for s in args.hairpin]
    #print(yy_maxdiff, hp_maxdiff, file=sys.stderr)

    if args.latex:
        print("""
\\documentclass[a4paper,10pt]{article}
\\usepackage[landscape,margin=5mm]{geometry}
\\usepackage{xcolor}

\\definecolor{bs}{RGB}{0, 0, 255} % blue
\\definecolor{nobs}{RGB}{255, 0, 0} % red
\\definecolor{mismatch}{RGB}{255, 0, 255} % magenta
\\definecolor{hairpin}{RGB}{255, 255, 0} % yellow
\\definecolor{yadapter}{RGB}{0, 255, 255} % cyan

\\setlength{\parindent}{0pt}

\\begin{document}
""")
        def f(s, cls):
            if cls is None:
                return s
            elif cls in ("hairpin", "yadapter"):
                # set as highlight colour
                return "{{\\begingroup\\setlength{{\\fboxsep}}{{0pt}}\\colorbox{{{}}}{{{}\\/}}\\endgroup}}".format(cls, s)
            else:
                return "{{\\color{{{}}}{}}}".format(cls, s)
    else:
        print("""
<!DOCTYPE html>
<html>
<head>
<style>
.bs {color:blue}
.nobs {color:red}
.mismatch {color:magenta}
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

    ss = qq = None
    fq1 = parse_fq(args.fq1)
    fq2 = None
    if args.interleaved:
        fq2 = fq1
    elif args.fq2:
        fq2 = parse_fq(args.fq2)

    seqno = 0

    while True:
        seqno += 1

        if fq2 is None:
            # folded reads
            try:
                label, comment, ss, qq = next(fq1)
            except StopIteration:
                break
            ss1, ss2, qq1, qq2, hp = comment2ssqq(comment)
            hlen = len(hp)

            hpi = rhpi = len(ss)
            p5i = p7i = 2*len(ss)+hlen
        else:
            # unfolded paired-end reads
            try:
                label1, comment1, ss1, qq1 = next(fq1)
                label2, comment2, ss2, qq2 = next(fq2)
            except StopIteration:
                break

            if label1[-2:] in ("/1", ".1"):
                label1 = label1[:-2]
                label2 = label2[:-2]
            if label1 != label2:
                print("read name mismatch for pair {}".format(seqno), file=sys.stderr)
                exit(1)

            label = label1
            hpi = rhpi = -1
            hp = ""

            # look for a hairpin
            for i, (h,rh) in enumerate(zip(args.hairpin, args.rhairpin)):
                hpi = mmfind(h, ss1, hp_maxdiff[i], min_matches=args.adapter_matchlen)
                rhpi = mmfind(rh, ss2, hp_maxdiff[i], min_matches=args.adapter_matchlen)
                #print(label, h, hpi, file=sys.stderr)
                if hpi != -1:
                    hp = h
                    break
                if rhpi != -1:
                    hp = h
                    break

            hlen = len(hp)

            if hpi != -1:
                p7i = mmfind(p7, ss1, yy_maxdiff, hpi+hlen, min_matches=args.adapter_matchlen)
            else:
                p7i = mmfind(p7, ss1, yy_maxdiff, min_matches=args.adapter_matchlen)

            if rhpi != -1:
                p5i = mmfind(p5, ss2, yy_maxdiff, rhpi+hlen, min_matches=args.adapter_matchlen)
            else:
                p5i = mmfind(p5, ss2, yy_maxdiff, min_matches=args.adapter_matchlen)


        if seqno > args.nseqs:
            break

        s1_list = []
        s2_list = []

        if args.latex:
            print("\\begin{samepage}\n{\\tt")
        else:
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
                s1_cls = s2_cls = "mismatch"

            if hpi != -1 and i>=hpi and i<hpi+hlen:
                s1_cls = "hairpin"
            if rhpi != -1 and i>=rhpi and i<rhpi+hlen:
                s2_cls = "hairpin"
            if p7i != -1 and i>=p7i:
                s1_cls = "yadapter"
            if p5i != -1 and i>=p5i:
                s2_cls = "yadapter"

            if q1 <= 20:
                if args.latex:
                    s1 = "{{\\underline{{{}}}\\/}}".format(s1)
                else:
                    s1 = "<u>{}</u>".format(s1)
            if q2 <= 20:
                if args.latex:
                    s2 = "{{\\underline{{{}}}\\/}}".format(s2)
                else:
                    s2 = "<u>{}</u>".format(s2)
            s1_list.append((s1, s1_cls))
            s2_list.append((s2, s2_cls))

        if args.latex:
            print(label, "\\\\")
            print("R1", g(s1_list), "\\\\")
            if ss is None:
                print("R2", g(s2_list))
            else:
                print("R2", g(s2_list), "\\\\")
                print("FS ", ss, "\\\\")
                qq = qq.replace("&", "\\&")
                qq = qq.replace("$", "\\$")
                print("FQ ", qq)
            print("}\\vskip\\baselineskip \\end{samepage}\n")
        else:
            print(label, "</br>")
            print("R1", g(s1_list), "</br>")
            print("R2", g(s2_list), "</br>")
            if ss is not None:
                print("FS ", ss, "</br>")
                print("FQ ", qq)
            print("</tt></p>")

    if args.latex:
        print("\n\\end{document}")
    else:
        print("</h2></html></body>")
