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

class ParseError(Exception):
    pass

def parse_nt_pairing(filename):

    ctx5p_pairs = None
    ctx5p_pos = []
    ctx5p = []

    ctx3p_pairs = None
    ctx3p_pos = []
    ctx3p = []

    with open(filename) as f1:

        for line in f1:
            line = line.rstrip()
            if not line:
                continue

            if line.startswith("#CTX5p"):
                ctx5p_pairs = line.split("\t")[2:]
                continue

            if line.startswith("CTX5p"):
                fields = [int(a) for a in line.split("\t")[1:]]
                ctx5p_pos.append(fields[0])
                x = np.array(fields[1:], dtype=float)
                sum_x = np.sum(x)
                if (sum_x):
                    x /= sum_x
                ctx5p.append(x)
                continue

            if line.startswith("#CTX3p"):
                ctx3p_pairs = line.split("\t")[2:]
                continue

            if line.startswith("CTX3p"):
                fields = [int(a) for a in line.split("\t")[1:]]
                ctx3p_pos.append(fields[0])
                x = np.array(fields[1:], dtype=float)
                sum_x = np.sum(x)
                if (sum_x):
                    x /= sum_x
                ctx3p.append(x)
                continue

    return ctx5p_pairs, ctx5p_pos, ctx5p, ctx3p_pairs, ctx3p_pos, ctx3p

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} pairs.txt plot.pdf".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    ctx5p_pairs, ctx5p_pos, ctx5p, ctx3p_pairs, ctx3p_pos, ctx3p = parse_nt_pairing(sys.argv[1])

    pal12 = ("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffdd5f", "#b15928", "black", "gray", "c", "m")
    markers = ('o', 'v', '^', '<', '>', '1', '2', '3',
                '4', 'd', 'p', 'h', '*', '+', 'x', '|')

    plot_file = sys.argv[2]
    pdf = PdfPages(plot_file)

    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))

    gs1 = gridspec.GridSpec(2, 7)
    ax1 = fig1.add_subplot(gs1[0,:3])
    ax2 = fig1.add_subplot(gs1[0,-3:], sharey=ax1)
    ax3 = fig1.add_subplot(gs1[1,:3], sharex=ax1)
    ax4 = fig1.add_subplot(gs1[1,-3:], sharey=ax3, sharex=ax2)
    ax2.invert_xaxis()
#    ax4.invert_xaxis()
    ax2.yaxis.tick_right()
    ax2.tick_params(right="on", left="on")
    ax4.yaxis.tick_right()
    ax4.tick_params(right="on", left="on")

    bpairs = set(ctx5p_pairs)
    if bpairs != set(ctx3p_pairs):
        print("Error: 5' and 3' pairings do not match.", file=sys.stderr)
        exit(1)

    alpha = 1.0
    plotlinestyle = ":"
    vlinestyle = "--"
    vlinecolour = "black"
    vlinealpha = 0.5

    if min(ctx5p_pos) < 0:
        ax1.axvline(0, color=vlinecolour, markeredgecolor=vlinecolour, linestyle=vlinestyle, alpha=vlinealpha)
        ax3.axvline(0, color=vlinecolour, markeredgecolor=vlinecolour, linestyle=vlinestyle, alpha=vlinealpha)
    if min(ctx3p_pos) < 0:
        ax2.axvline(-1, color=vlinecolour, markeredgecolor=vlinecolour, linestyle=vlinestyle, alpha=vlinealpha)
        ax4.axvline(-1, color=vlinecolour, markeredgecolor=vlinecolour, linestyle=vlinestyle, alpha=vlinealpha)

    for i, (c, m, label_5p, label_3p) in enumerate(zip(pal12, markers, ctx5p_pairs, ctx3p_pairs)):
        pairs_5p = [ctx[i] for ctx in ctx5p]
        pairs_3p = [ctx[i] for ctx in ctx3p]

        if np.mean(pairs_5p) > 0.1:
            ax5p = ax1
        else:
            ax5p = ax3

        if np.mean(pairs_3p) > 0.1:
            ax3p = ax2
        else:
            ax3p = ax4

        ax5p.plot(ctx5p_pos, pairs_5p, color=c, markeredgecolor=c, marker=m, linestyle=plotlinestyle, alpha=alpha, label=label_5p)
        ax3p.plot(ctx3p_pos, pairs_3p, color=c, markeredgecolor=c, marker=m, linestyle=plotlinestyle, alpha=alpha, label=label_3p)

    #ax3.set_xlim(min(ctx5p_pos)-0.5, max(ctx5p_pos)+0.5)
    #ax4.set_xlim(min(ctx3p_pos)-0.5, max(ctx3p_pos)+0.5)
    ax3.set_xlim(-20.5, 25.5)
    ax4.set_xlim(-20.5, 25.5)

    setp(ax1.get_xticklabels(), visible=False)
    #setp(ax2.get_yticklabels(), visible=False)
    setp(ax2.get_xticklabels(), visible=False)
    #setp(ax4.get_yticklabels(), visible=False)
    ax4.yaxis.set_label_position("right")

    ax1_handles,ax1_labels = ax1.get_legend_handles_labels()
    empty = mpatches.Patch(color='white')
    handles = [empty] + ax1_handles
    labels = ['+/-'] + ax1_labels

    ax1.legend(handles, labels, numpoints=1, frameon=False, loc='center left', bbox_to_anchor=(0.955, 0.75))
    ax3.legend(numpoints=1, frameon=False, loc='center left', bbox_to_anchor=(0.955, 0.75))

    ax3.set_xlabel("Distance from p5/p7")
    ax4.set_xlabel("Distance from hairpin")
    ax1.set_ylabel("Frequency", labelpad=20, y=-0.15)
    fig1.suptitle("Observed pairing frequencies")

    #plt.tight_layout()
    pdf.savefig()
    pdf.close()
