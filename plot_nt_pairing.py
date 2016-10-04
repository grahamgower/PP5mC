#!/usr/bin/env python

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
                if bpairs != set(ctx5p_pairs):
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
                if bpairs != set(ctx3p_pairs):
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
    parser = argparse.ArgumentParser(description="plot nucleotide pairing info")
    parser.add_argument("--only5p", action="store_true", default=False, help="only plot 5' end")
    parser.add_argument("--scale", type=int, default=1, help="scale the plot")
    parser.add_argument("--title", help="text for title")
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

    pairlist = ("A/T", "T/A", "G/T", "T/G", "G/C", "C/G",
            "A/C", "C/A", "A/G", "G/A", "C/T", "T/C",
            "A/A", "C/C", "G/G", "T/T")
    pair2sym = {
            "A/T": ('o', "cyan", "cyan"),
            "T/A": ('o', "orange", "orange"),
            "G/T": ('*', "blue", "blue"),
            "T/G": ('*', "yellow", "blue"),
            "G/C": ('s', "red", "red"),
            "C/G": ('d', "green", "green"),

            "A/C": ('v', "cyan", "yellow"),
            "C/A": ('v', "yellow", "green"),
            "A/G": ('^', "blue", "tan"),
            "G/A": ('^', "orange", "blue"),
            "C/T": ('>', "red", "lightgreen"),
            "T/C": ('>', "cyan", "brown"),

            "A/A": ('<', "tan", "brown"),
            "C/C": ('x', "pink", "red"),
            "G/G": ('h', "pink", "blue"),
            "T/T": ('p', "tan", "green")
            }

    ctx5p_pos, ctx5p, ctx3p_pos, ctx3p = parse_nt_pairing(args.infile, set(pairlist))
    #ctx3p_pos = list(range(29,-31,-1))

    plot_file = args.outfile
    pdf = PdfPages(plot_file)

    scale = 1
    #fig_w, fig_h = plt.figaspect(9.0/32.0)
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(scale*fig_w, scale*fig_h))
    #fig1 = plt.figure(figsize=(fig_w, fig_h))
    #fig1 = plt.figure(figsize=(scale*11.69, scale*8.27), dpi=100)

    gs1 = gridspec.GridSpec(2, 9)
    if args.only5p:
        ax1 = fig1.add_subplot(gs1[0,:8])
        ax3 = fig1.add_subplot(gs1[1,:8], sharex=ax1)
    else:
        ax1 = fig1.add_subplot(gs1[0,:4])
        ax2 = fig1.add_subplot(gs1[0,-4:], sharey=ax1)
        ax3 = fig1.add_subplot(gs1[1,:4], sharex=ax1)
        ax4 = fig1.add_subplot(gs1[1,-4:], sharey=ax3, sharex=ax2)
        ax2.invert_xaxis()
    #    ax4.invert_xaxis()
        ax2.yaxis.tick_right()
        ax2.tick_params(right="on", left="on")
        ax4.yaxis.tick_right()
        ax4.tick_params(right="on", left="on")

    alpha = 1.0
    linestyle = ":"
    linewidth = scale*1.0
    markeredgewidth = scale*0.5
    markersize = scale*5
    vlinestyle = "--"
    vlinecolour = "black"
    vlinealpha = 0.5

    if min(ctx5p_pos) < 0:
        ax1.axvline(0, color=vlinecolour, markeredgecolor=vlinecolour, linestyle=vlinestyle, alpha=vlinealpha)
        ax3.axvline(0, color=vlinecolour, markeredgecolor=vlinecolour, linestyle=vlinestyle, alpha=vlinealpha)
    if not args.only5p and min(ctx3p_pos) < 0:
        ax2.axvline(-1, color=vlinecolour, markeredgecolor=vlinecolour, linestyle=vlinestyle, alpha=vlinealpha)
        ax4.axvline(-1, color=vlinecolour, markeredgecolor=vlinecolour, linestyle=vlinestyle, alpha=vlinealpha)
    
    for i, label in enumerate(pairlist):
        pairs_5p = ctx5p[label]
        pairs_3p = ctx3p[label]
        c1 = pal16d[pair2sym[label][1]]
        c2 = pal16d[pair2sym[label][2]]
        m = pair2sym[label][0]

        if np.mean(pairs_5p) > 0.1:
            ax5p = ax1
        else:
            ax5p = ax3

        if not args.only5p:
            if np.mean(pairs_3p) > 0.1:
                ax3p = ax2
            else:
                ax3p = ax4

        ax5p.plot(ctx5p_pos, pairs_5p, color=c1, markeredgecolor=c2,
                marker=m, markersize=markersize, markeredgewidth=markeredgewidth,
                linestyle=linestyle, linewidth=linewidth,
                alpha=alpha, label=label)
        if not args.only5p:
            ax3p.plot(ctx3p_pos, pairs_3p, color=c1, markeredgecolor=c2,
                    marker=m, markersize=markersize, markeredgewidth=markeredgewidth,
                    linestyle=linestyle, linewidth=linewidth,
                    alpha=alpha, label=label)

    ax3.set_xlim(min(ctx5p_pos)-0.5, max(ctx5p_pos)+0.5)
    #ax3.set_xlim(-20.5, 25.5)
    setp(ax1.get_xticklabels(), visible=False)
    #setp(ax2.get_yticklabels(), visible=False)

    if not args.only5p:
        ax4.set_xlim(min(ctx3p_pos)-0.5, max(ctx3p_pos)+0.5)
        #ax4.set_xlim(-20.5, 25.5)

        setp(ax2.get_xticklabels(), visible=False)
        #setp(ax4.get_yticklabels(), visible=False)
        ax4.yaxis.set_label_position("right")

    ax1_handles,ax1_labels = ax1.get_legend_handles_labels()
    empty = mpatches.Patch(color='white')
    handles = [empty] + ax1_handles
    labels = ['+/-'] + ax1_labels

    ax1.legend(handles, labels, numpoints=1, frameon=False, loc='center left', bbox_to_anchor=(1.01, 0.75))
    ax3.legend(numpoints=1, frameon=False, loc='center left', bbox_to_anchor=(1.01, 0.75))

    ax3.set_xlabel("Distance from p5/p7")
    if not args.only5p:
        ax4.set_xlabel("Distance from hairpin")
    ax1.set_ylabel("Frequency", labelpad=20, y=-0.15)

    if args.title:
        title = "Observed pairing frequencies ({})".format(args.title)
    else:
        title = "Observed pairing frequencies"
    fig1.suptitle(title, fontsize=int(np.sqrt(scale)*12))

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    pdf.savefig()
#    pdf.savefig(papertype="a4", orientation="landscape")
    pdf.close()
