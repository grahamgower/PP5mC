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

def parse_fold_log(filename):

    pos_ds, pos_l, pos_r = None, None, None
    pairs_ds, pairs_l, pairs_r = {}, {}, {}

    with open(filename) as f1:

        for line in f1:
            if line.startswith("Watson<->Crick [>"):
                pos_ds = True
                break

        if pos_ds is None:
            raise ParseError("Could not find Watson<->Crick [>Nbp from fragment end] line")

        for line in f1:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("Watson<->Crick [5']"):
                pos_l = [int(pp) for pp in line.split("\t")[1:]]
                break
            fields = line.split(": ")
            pairs_ds[fields[0]] = float(fields[1])

        if pos_l is None:
            raise ParseError("Could not find Watson<->Crick [5'] line")

        for line in f1:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("Watson<->Crick [3']"):
                pos_r = [int(pp) for pp in line.split("\t")[1:]]
                break
            fields = line.split("\t")
            pairs_l[fields[0]] = [float(ff) for ff in fields[1:]]

        if pos_r is None:
            raise ParseError("Could not find Watson<->Crick [3'] line")

        for line in f1:
            line = line.rstrip()
            if not line:
                continue
            fields = line.split("\t")
            pairs_r[fields[0]] = [float(ff) for ff in fields[1:]]

    return pos_l, pos_r, pairs_l, pairs_r, pairs_ds

if __name__ == "__main__":
    if len(sys.argv) not in (2,3):
        print("usage: {} fold_out.log [plot.pdf]".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    pos_l, pos_r, pairs_l, pairs_r, pairs_ds = parse_fold_log(sys.argv[1])

    pal12 = ("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffdd5f", "#b15928", "black", "gray", "c", "m")
    markers = ('o', 'v', '^', '<', '>', '1', '2', '3',
                '4', 'd', 'p', 'h', '*', '+', 'x', '|')

    plot_file = sys.argv[2] if len(sys.argv) == 3 else "plot_fold_log.pdf"
    pdf = PdfPages(plot_file)

    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig_w, fig_h = plt.figaspect(3.0/4.0)
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

    bpairs = set(pairs_l.keys())
    if bpairs != set(pairs_r.keys()):
        print("Error: left and right pairings do not match.", file=sys.stderr)
        exit(1)

#    xx = 25
#    pos_l = pos_l[:xx]
#    pos_r = pos_r[:xx]
#    for pair in bpairs:
#        pairs_l[pair] = pairs_l[pair][:xx]
#        pairs_r[pair] = pairs_r[pair][:xx]

    alpha = 1.0
    hlinestyle = "none"
    plotlinestyle = ":"

    for c, m, pair in zip(pal12, markers, bpairs):
        b0,b1 = pair.split("<->")
        label = b0 + "/" + b1

        if np.mean(pairs_l[pair]) > 0.1:
            ax1.plot(pos_l, pairs_l[pair], color=c, markeredgecolor=c, marker=m, linestyle=plotlinestyle, alpha=alpha, label=label)
            ax1.axhline(pairs_ds[pair], color=c, markeredgecolor=c, linestyle=hlinestyle)
        else:
            ax3.plot(pos_l, pairs_l[pair], color=c, markeredgecolor=c, marker=m, linestyle=plotlinestyle, alpha=alpha, label=label)
            ax3.axhline(pairs_ds[pair], color=c, markeredgecolor=c, linestyle=hlinestyle)

        if np.mean(pairs_r[pair]) > 0.1:
            ax2.plot(pos_r, pairs_r[pair], color=c, markeredgecolor=c, marker=m, linestyle=plotlinestyle, alpha=alpha, label=label)
            ax2.axhline(pairs_ds[pair], color=c, markeredgecolor=c, linestyle=hlinestyle)
        else:
            ax4.plot(pos_r, pairs_r[pair], color=c, markeredgecolor=c, marker=m, linestyle=plotlinestyle, alpha=alpha, label=label)
            ax4.axhline(pairs_ds[pair], color=c, markeredgecolor=c, linestyle=hlinestyle)

    ax3.set_xlim(-0.5, max(pos_l)+0.5)
    ax4.set_xlim(max(pos_r)+0.5, -0.5)

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
