#!/usr/bin/env python
#
# Plot nucleotide pairing info, as output from `scanbp'.
#
# Copyright (c) 2016-2018 Graham Gower <graham.gower@gmail.com>
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
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from pylab import setp
import numpy as np
import collections


class ParseError(Exception):
    pass

def parse_nt_pairing(filename, bpairs):

    ctx5p_pairs = None
    ctx5p_pos = []
    ctx5p = collections.defaultdict(list)

    ctx3p_pairs = None
    ctx3p_pos = []
    ctx3p = collections.defaultdict(list)

    with open(filename) as f1:

        for lineno, line in enumerate(f1,1):
            line = line.rstrip()
            if not line:
                continue

            if line.startswith("#CTX5p"):
                ctx5p_pairs = line.split("\t")[2:]
                if not bpairs.issubset(set(ctx5p_pairs)):
                    raise ParseError("{}: line {}: missing pairs.".format(filename, lineno))

                continue

            if line.startswith("CTX5p"):
                fields = [int(a) for a in line.split("\t")[1:]]
                ctx5p_pos.append(fields[0])
                x = np.array(fields[1:], dtype=float)
                sum_x = np.sum(x)
                if (sum_x):
                    x /= sum_x

                for i, xx in enumerate(x):
                    ctx5p[ctx5p_pairs[i]].append(xx)

                continue

            if line.startswith("#CTX3p"):
                ctx3p_pairs = line.split("\t")[2:]
                if not bpairs.issubset(set(ctx3p_pairs)):
                    raise ParseError("{}: line {}: missing pairs.".format(filename, lineno))
                continue

            if line.startswith("CTX3p"):
                fields = [int(a) for a in line.split("\t")[1:]]
                ctx3p_pos.append(fields[0])
                x = np.array(fields[1:], dtype=float)
                sum_x = np.sum(x)
                if (sum_x):
                    x /= sum_x

                for i, xx in enumerate(x):
                    ctx3p[ctx3p_pairs[i]].append(xx)

                continue

    return  ctx5p_pos, ctx5p, ctx3p_pos, ctx3p

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="plot nucleotide pairing info, as output from `scanbp'")
    parser.add_argument("--only5p", action="store_true", default=False, help="only plot 5' end")
    parser.add_argument("--scale", type=float, default=1, help="scale the plot")
    parser.add_argument("--title", help="text for title")
    parser.add_argument("--ratio4x3", action="store_true", default=False, help="plot 4x3 ratio")
    parser.add_argument("--allpairs", action="store_true", default=False, help="plot all base pairs, not just canonical pairs")
    parser.add_argument("infile", help="input pairs.txt")
    parser.add_argument("outfile", help="output pairs.pdf")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    pal16d = {"black":"000000",
            "darkgray":"#575757",
            "red":"#ad2323",
            "blue":"#2a4bd7",
            "green":"#1d6914",
            "brown":"#814a19",
            "purple":"#8126c0",
            "lightgray":"#a0a0a0",
            "lightgreen":"#81c57a",
            "lightblue":"#9dafff",
            "cyan":"#29d0d0",
            "orange":"#ff9233",
            "yellow":"#eedd22",
            "tan":"#e9debb",
            "pink":"#ffcdf3",
            "white":"#ffffff"}

    if args.allpairs:
        pairlist = ("A/T", "T/A", "G/T", "T/G", "G/C", "C/G",
                "A/C", "C/A", "A/G", "G/A", "C/T", "T/C",
                "A/A", "C/C", "G/G", "T/T")
        labels = pairlist
    else:
        pairlist = ("A/T", "T/A", "G/T", "T/G", "G/C", "C/G")
        labels = ("A/T", "T/A", "G/C", "C/G", "G/mC", "mC/G")

    pair2sym = {
            "A/T": ('o', "cyan", "cyan"),
            "T/A": ('v', "orange", "orange"),
            "G/T": ('*', "green", "green"),
            "T/G": ('p', "yellow", "blue"),
            "G/C": ('s', "tan", "brown"),
            "C/G": ('d', "pink", "red"),

            "A/C": ('v', "cyan", "red"),
            "C/A": ('v', "yellow", "green"),
            "A/G": ('^', "blue", "tan"),
            "G/A": ('^', "orange", "blue"),
            "C/T": ('>', "red", "red"),
            "T/C": ('>', "cyan", "brown"),

            "A/A": ('<', "tan", "orange"),
            "C/C": ('x', "black", "black"),
            "G/G": ('h', "pink", "blue"),
            "T/T": ('p', "tan", "green")
            }

    ctx5p_pos, ctx5p, ctx3p_pos, ctx3p = parse_nt_pairing(args.infile, set(pairlist))
    ctx3p_pos = np.array(ctx3p_pos)+1

    plot_file = args.outfile
    if plot_file.endswith("pdf"):
        pdf = PdfPages(plot_file)

    if args.ratio4x3:
        fig_w, fig_h = plt.figaspect(3.0/4.0)
    else:
        fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig1 = plt.figure(figsize=(args.scale*fig_w, args.scale*fig_h))

    gs1 = gridspec.GridSpec(10, 8)
    if args.only5p:
        ax1 = fig1.add_subplot(gs1[:9,:8])
    else:
        ax1 = fig1.add_subplot(gs1[:9,:4])
        ax2 = fig1.add_subplot(gs1[:9,-4:], sharey=ax1)
        ax2.invert_xaxis()
        ax2.yaxis.tick_right()
        ax2.tick_params(right="on", left="on")

    ax1.tick_params(right="on", left="on")

    alpha = 1.0
    linestyle = ":"
    linewidth = args.scale*1.0
    markeredgewidth = args.scale*0.5
    markersize = args.scale*5
    vlinestyle = "--"
    vlinecolour = "black"
    vlinealpha = 0.5

    if min(ctx5p_pos) < 0:
        ax1.axvline(0, color=vlinecolour, markeredgecolor=vlinecolour, linestyle=vlinestyle, alpha=vlinealpha)
    if not args.only5p and min(ctx3p_pos) < 0:
        ax2.axvline(0, color=vlinecolour, markeredgecolor=vlinecolour, linestyle=vlinestyle, alpha=vlinealpha)

    for i, (pair, label) in enumerate(zip(pairlist,labels)):
        pairs_5p = ctx5p[pair]
        pairs_3p = ctx3p[pair]
        c1 = pal16d[pair2sym[pair][1]]
        c2 = pal16d[pair2sym[pair][2]]
        m = pair2sym[pair][0]

        ax1.plot(ctx5p_pos, pairs_5p, color=c2, markerfacecolor=c1, markeredgecolor=c2,
                marker=m, markersize=markersize, markeredgewidth=markeredgewidth,
                linestyle=linestyle, linewidth=linewidth,
                alpha=alpha, label=label)
        if not args.only5p:
            ax2.plot(ctx3p_pos, pairs_3p, color=c2, markerfacecolor=c1, markeredgecolor=c2,
                    marker=m, markersize=markersize, markeredgewidth=markeredgewidth,
                    linestyle=linestyle, linewidth=linewidth,
                    alpha=alpha, label=label)

    ax1.set_xlim(min(ctx5p_pos)-0.5, max(ctx5p_pos)+0.5)
    ax2.set_xlim(min(ctx3p_pos)-0.5, max(ctx3p_pos)+0.5)

    if not args.only5p:
        ax2.set_xlim(min(ctx3p_pos)-0.5, max(ctx3p_pos)+0.5)
        ax2.yaxis.set_label_position("right")

    ax1.set_xlabel("Distance from $5'$ end of $+$ strand", labelpad=10)
    if not args.only5p:
        ax2.set_xlabel("Distance from $3'$ end of $+$ strand", labelpad=10)
    ax1.set_ylabel("Frequency", labelpad=10)

    ax1_handles,ax1_labels = ax1.get_legend_handles_labels()
    empty = mpatches.Patch(color='white')
    handles = [empty] + ax1_handles
    labels = ["+/-"] + ax1_labels

    def flip(items, ncol):
        import itertools
        return list(itertools.chain(*[items[i::ncol] for i in range(ncol)]))

    if args.allpairs:
        handles = [empty,empty,empty] + flip(ax1_handles, 6)
        labels = ["+/-", "", ""] + flip(ax1_labels, 6)
    else:
        handles = [empty] + ax1_handles
        labels = ["+/-"] + ax1_labels

    leg = ax1.legend(handles, labels, numpoints=1, frameon=False, ncol=7, loc='upper center', bbox_to_anchor=(1.00, -0.16))
    pm_text = leg.get_texts()[0]
    pm_text.set_ha('center')

    if args.title:
        title = "Base pair frequencies ({})".format(args.title)
    else:
        title = "Base pair frequencies"
    fig1.suptitle(title, fontsize=int(np.sqrt(args.scale)*12))

    if args.allpairs:
        plt.tight_layout(rect=[0, 0.15, 0.95, 0.95])
    else:
        plt.tight_layout(rect=[0, 0.06, 0.95, 0.95])

    if plot_file.endswith("pdf"):
        pdf.savefig()
        pdf.close()
    else:
        plt.savefig(plot_file)
