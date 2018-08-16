#!/usr/bin/env python
#
# Plot quality profile(s), as output from `qualprofile'.
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
import os.path
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np

def parse_qualprofile(fn):
    V = []
    i = 0
    with open(fn) as f:
        for line in f:
            line = line.rstrip()
            if not line or line[0] == "#":
                continue

            if line.startswith("MEAN"):
                X = map(float, line.split()[1:])
                Sigma = np.zeros([len(X),len(X)])
            elif line.startswith("COV"):
                Si = map(float, line.split()[1:])
                Sigma[i,:len(Si)] = Si
                V.append(Si[-1])
                i += 1

    if len(X) != len(V):
        raise Exception("{}: MEAN and COV length mismatch".format(fn))

    return np.array(X), np.array(V), Sigma

def ribbonplot(ax, x, y, yerr, colour, label):
    ax.plot(x, y, label=label, color=colour, lw=1)
    ax.plot(x, y+yerr, color=colour, lw=0.5)
    ax.plot(x, y-yerr, color=colour, lw=0.5)
    ax.fill_between(x, y-yerr, y+yerr, color=colour, alpha=0.3)
    ax.set_ylim(0,np.max(y+yerr))

if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print("{}: out.pdf profile1 [profile2]".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    X1, V1, Sigma1 = parse_qualprofile(sys.argv[2])
    fn1 = os.path.basename(sys.argv[2])
    if len(sys.argv) == 4:
        X2, V2, Sigma2 = parse_qualprofile(sys.argv[3])
        fn2 = os.path.basename(sys.argv[3])
    else:
        X2 = V2 = Sigma2 = fn2 = None

    plot_file = sys.argv[1]
    if plot_file.endswith("pdf"):
        pdf = PdfPages(plot_file)
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    scale = 1
    fig1 = plt.figure(figsize=(scale*fig_w, scale*fig_h))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = fig1.add_subplot(gs1[0])

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    h = np.arange(1,len(X1)+1)
    ribbonplot(ax1, h, X1, np.sqrt(V1), colors[0], fn1)
    if X2 is not None:
        ribbonplot(ax1, h, X2, np.sqrt(V2), colors[1], fn2)

    ax1.set_xlabel("Position in read")
    ax1.set_ylabel("Qual score")
    ax1.legend()

    plt.tight_layout()
    if plot_file.endswith("pdf"):
        pdf.savefig(figure=fig1)
    else:
        plt.savefig(plot_file[:-4]+"-1"+plot_file[-4:])

    # next page
    fig2 = plt.figure(figsize=(scale*fig_w, scale*fig_h))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = fig2.add_subplot(gs1[0])

    im  = ax1.imshow(Sigma1, cmap='magma')
    ax1.figure.colorbar(im, ax=ax1)
    ax1.set_title("Covariance matrix ({})".format(fn1))

    plt.tight_layout()
    if plot_file.endswith("pdf"):
        pdf.savefig(figure=fig2)
    else:
        plt.savefig(plot_file[:-4]+"-2"+plot_file[-4:])

    if Sigma2 is not None:
        # next page
        fig3 = plt.figure(figsize=(scale*fig_w, scale*fig_h))
        gs1 = gridspec.GridSpec(1, 1)
        ax1 = fig3.add_subplot(gs1[0])

        im  = ax1.imshow(Sigma2, cmap='magma')
        ax1.figure.colorbar(im, ax=ax1)
        ax1.set_title("Covariance matrix ({})".format(fn2))

        plt.tight_layout()
        if plot_file.endswith("pdf"):
            pdf.savefig(figure=fig3)
        else:
            plt.savefig(plot_file[:-4]+"-3"+plot_file[-4:])

    if plot_file.endswith("pdf"):
        pdf.close()
