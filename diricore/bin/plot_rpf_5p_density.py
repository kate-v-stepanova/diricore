#!/usr/bin/env python
#

"""
"""

from __future__ import absolute_import
import os
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties
import seaborn
seaborn.reset_orig()
import numpy
from matplotlib.ticker import FormatStrFormatter

sys.path.append(os.path.join(os.path.dirname(__file__), '../lib/'))


from plottxcoords import get_means, read_sampleinfo


SUBPLOTSIZE = (6, 6 * 10 / 16.)
SUBPLOTPARS = {"left": 0.125, "right": 0.9, "top": 0.9, "bottom": 0.1, "wspace": 0.2, "hspace": 0.2}

XTICKFONT = FontProperties(weight="bold", size=18)
YTICKFONT = FontProperties(weight="bold", size=12)
LABELFONT = FontProperties(weight="bold", size=13)
TITLEFONT = FontProperties(weight="bold", size=22)


def calc_figsize(fignrc, subplotsize, subplotpars):
    # fignrc is (nrow, ncol)
    # but subplotsize is (x, y)
    x = fignrc[1] * subplotsize[0] / \
            (subplotpars.get("right", 1.0) - subplotpars.get("left", 0.0) - subplotpars.get("wspace", 0.0))
    y = fignrc[0] * subplotsize[1] / \
            (subplotpars.get("top", 1.0) - subplotpars.get("bottom", 0.0) - subplotpars.get("hspace", 0.0))

    return (x, y)


def plot_5p_rpf_density_difference(h5fn, pairs, codongroups, h5mapsfn, min_reads=100, ylim=None, prettynames=None, centered_at_x0=False, range=None):
    # assume samples is list of (cond, ref, color)
    # assume codondata is map of aa -> indexname
    if range is None:
        min_x = -30
        max_x = 31
    if centered_at_x0:
        X = numpy.arange(min_x, max_x)
    else:
        X = numpy.arange(min_x, max_x) + 1

    figsize = calc_figsize((len(pairs), len(codongroups)), SUBPLOTSIZE, SUBPLOTPARS)
    fig, axs = pyplot.subplots(len(pairs), len(codongroups), sharex=True, sharey=True, squeeze=False, figsize=figsize)
    for i_pair, pairdata in enumerate(pairs):
        refsampleid, condsampleid, linecolor = pairdata
        #samples = [(h5fn, condsampleid, None), (h5fn, refsampleid, None)] # this should be correct. because plottxcoords.py.
        samples = [(h5fn, refsampleid, None), (h5fn, condsampleid, None)]
        means, passing_cutoff = get_means(samples, codongroups, h5mapsfn, min_reads, intersect=True, smooth3nt=True)

        if len(passing_cutoff) != len(codongroups):
             print("passing_cutoff != len(codongroupgs). {} != {}".format(passing_cutoff, len(codongroups)))

        if means.shape[1] != len(codongroups):
             print("means.shape[1] != len(codongroups). {} != {}".format(means.shape[1], len(codongroups)))

        for j, codongroup in enumerate(codongroups):
            ws = passing_cutoff[j]
            assert len(ws) == 2
            assert (ws[1] == ws[0]).all(), "Assumed if intersect==True"

            ax = axs[i_pair][j]

            ax.plot(X + 0.5, means[1, j] - means[0, j], color=linecolor, linewidth=3.0, alpha=0.65)

            ax.set_xticks((X + 0.5)[X % 3 == 0])
            ax.set_xticklabels(map(str, X[X % 3 == 0]), fontproperties=XTICKFONT)
            ax.xaxis.set_ticks_position('bottom')
            ax.set_xlim(X.min(), X.max())

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()

            ax.set_title("%s" % (",".join(sorted(codongroup)), ), fontproperties=TITLEFONT)

            ylabel = "%s vs %s" % (refsampleid, condsampleid)

            if prettynames is not None:
                ylabel = "%s vs %s" % (prettynames[condsampleid], prettynames[refsampleid])

            ax.set_ylabel(ylabel, ha="center", va="center", labelpad=20, fontproperties=LABELFONT)

            ax.axhline(0.0, linewidth=1.25, linestyle="--", color="#3c3c3c")

            ax.text(0.985, 0.98, "$n = %d$" % ws[0].sum(), ha="right", va="top", transform=ax.transAxes, fontsize=13, weight="bold")

    if ylim is not None:
        axs[0][0].set_ylim(ylim)

    for i in xrange(len(pairs)):
        for j in xrange(len(codongroups)):
            ax = axs[i][j]
            pyplot.setp(ax.get_yticklabels(), visible=True)
            ax.set_yticklabels(ax.get_yticks(), fontproperties=YTICKFONT)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')) # 2 decimals

    for i in xrange(len(pairs)):
        for j in xrange(len(codongroups)):
            ax = axs[i][j]
            rectheight = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.025
            rectbottom = ax.get_ylim()[0] + 0.4 * rectheight
            ax.add_patch(Rectangle((0.5, rectbottom), 3, rectheight, facecolor="black", edgecolor="none"))
    return fig


def main():
    import argparse
    ap = argparse.ArgumentParser(description="")
    ap.add_argument('--y-limits', default=None, required=False)
    ap.add_argument('--centered-at-x0', action="store_true", default=False, required=False)
    ap.add_argument('-c', '--contrasts', required=False)
    ap.add_argument('-n', '--sample-names', default=None, required=False)
    ap.add_argument('-o', '--outfile', metavar="PDF", required=False)
    ap.add_argument('-m', '--min-reads', type=int, default=100, required=False)
    ap.add_argument(dest="h5file", metavar="HDF5FILE", nargs=1)
    ap.add_argument(dest="h5mapsfile", metavar="HDF5MAPSFILE", nargs=1)
    ap.add_argument(dest="codongroups", metavar="CODON,CODON,[...]", nargs="+")

    args = ap.parse_args()

    h5mapsfn = args.h5mapsfile[0]
    pairs = read_sampleinfo(args.contrasts)
    if args.sample_names is not None:
        sampleinfo = read_sampleinfo(args.sample_names, cols_num=2)
        sample_names = dict(sampleinfo)
        for (refsampleid, condsampleid, color) in pairs:
            print_samples = False
            if refsampleid not in sample_names:
                 print('WARNING: {} not in sample_names'.format(refsampleid))
                 print_samples = True
            if condsampleid not in sample_names:
                 print('WARNING: {} not in sample_names'.format(condsampleid))
                 print_samples = True
            if print_samples:
                 print(sample_names)
    else:
        sample_names = None

    ylim = args.y_limits
    if ylim is not None:
        ylim = map(float, ylim.split(","))[:2]

    codongroups = map(lambda s: set(s.split(",")), args.codongroups)
    print(codongroups)
    # Note (APU): condongroups must be a feature to be able to analyze two codons together
    # But in general, this will generate a list  containing single codon sets
    # If for example the input is "TTT,CCC ATG" , it will make: ( ('TTT','CCC') , (ATG) 

    fig = plot_5p_rpf_density_difference(args.h5file[0], pairs, codongroups, h5mapsfn, args.min_reads, ylim, sample_names, args.centered_at_x0)

    fig.savefig(args.outfile, format="pdf")

    return


if __name__ == "__main__":
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

    try:
        main()
    except IOError as e:
        if e.errno != 32:
            raise
    except KeyboardInterrupt as e:
        pass
