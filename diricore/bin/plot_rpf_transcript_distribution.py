#!/usr/bin/env python
#

import os
from itertools import izip
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')

import h5py
from numpy import array, zeros, empty, arange, diff, concatenate, interp, linspace
from numpy import convolve, apply_along_axis
import scipy.signal

import seaborn
from matplotlib import pyplot


DEFAULT_NGRID = 2000
DEFAULT_MINREADS = 100


def select_representative_isoforms(infofn):
    """
    Longest CDS per gene, if tied: least # of exons, if tied: longest transcript
    """

    info = h5py.File(infofn, 'r')

    tids = info["tid"][:]
    gids = info["gid"][:]
    cdspos = info["cdspos"][:]
    coordlut = info["coordlut"][:]
    tidlengths = diff(concatenate((info["seqlut"][:], array([info["seq"].shape[0]], dtype=info["seqlut"].dtype))))
    cdslengths = diff(cdspos, axis=1).squeeze()

    gi = defaultdict(lambda: list())
    for tid, gid, (nexons, exonsoffset), tidlength, cdslength in izip(tids, gids, coordlut, tidlengths, cdslengths):
        gi[gid].append((tid, cdslength, nexons, tidlength))

    chosentids = dict()
    for gid, vals in gi.iteritems():
        chosen = sorted(vals, key=lambda t: (t[1], -t[2], t[3]))[-1]
        if chosen[1] == 0:
            # cdslength == 0
            continue
        chosentids[gid] = chosen[0]

    return chosentids


def get_transcript_rpf_density(A, ngrid, infofn, chosentids):
    """
    WARNING: may allocate a lot of memory
    (Usage =~ ngrid * len(chosentids) * 3 * 4 bytes, which typically comes down to 0.5-1 Gb)
    """

    info = h5py.File(infofn, 'r')

    seqlut = info["seqlut"][:]
    tids = info["tid"][:]
    cdspos = info["cdspos"][:]

    tidlengths = concatenate((diff(seqlut), array([int(info["seq"].shape[0] - seqlut[-1])], dtype=seqlut.dtype)))
    assert tidlengths.sum() == info["seq"].shape[0]
    cdslengths = diff(cdspos).squeeze()
    assert (cdslengths <= tidlengths).all()

    datatid2idx = {x: i for (i, x) in enumerate(tids)}
    tid2itid = {x: i for (i, x) in enumerate(sorted(chosentids.itervalues()))}

    X = linspace(0, 1, ngrid + 1)
    interpcounts = zeros((len(chosentids), 3, ngrid), dtype="float32")
    counts = zeros(len(chosentids), dtype="uint32")
    for gid in chosentids:
        tid = chosentids[gid]
        tididx = datatid2idx[tid]
        offset = seqlut[tididx]
        l = tidlengths[tididx]
        cdsstart, cdsend = cdspos[tididx]
        assert cdsstart < l
        assert cdsstart < cdsend

        a = A[offset:(offset + l)]
        assert a.size == l
        t = a.sum()

        itid = tid2itid[tid]

        counts[itid] = t

        if t == 0:
            continue

        t = float(t)

        if (cdsstart > 0) and (a[:cdsstart].sum() > 0):
            interpcounts[itid, 0] = diff(interp(X, linspace(0, 1, cdsstart, endpoint=True), a[:cdsstart].cumsum() / t))

        if ((l - cdsend) > 0) and (a[cdsend:].sum() > 0):
            interpcounts[itid, 2] = diff(interp(X, linspace(0, 1, (l - cdsend), endpoint=True), a[cdsend:].cumsum() / t))

        if (a[cdsstart:cdsend].sum() > 0):
            interpcounts[itid, 1] = diff(interp(X, linspace(0, 1, (cdsend - cdsstart), endpoint=True), a[cdsstart:cdsend].cumsum() / t))

    return interpcounts, counts


def get_averaged_transcript_rpf_densities_mutliple_samples(sampleinfo, infofn, ngrid, minreads):
    chosentids = select_representative_isoforms(infofn)
    means = empty((len(sampleinfo), 3, ngrid), dtype="float32")
    ws = empty(len(sampleinfo), dtype=int)
    for isample, sample in enumerate(sampleinfo):
        h5data = h5py.File(sample[1], 'r')
        A = h5data["counts"][sample[0]][:]
        interpcounts, counts = get_transcript_rpf_density(A, ngrid, infofn, chosentids)
        w = counts >= minreads

        means[isample] = interpcounts[w].mean(axis=0) * 3 * ngrid
        ws[isample] = w.sum()

    W = scipy.signal.windows.gaussian(int(round(ngrid / 25.)), int(round(ngrid / 25.)) / 6.)
    W /= W.sum()

    v = concatenate((means[:, 0], means[:, 1], means[:, 2]), axis=1).T
    vc = apply_along_axis(lambda a: convolve(a, W, mode="same"), 0, v)

    return vc, ws


def plot_average_transcript_rpf_densities(sampleinfo, vc, ws):
    fig, ax = pyplot.subplots(1, 1, figsize=(10, 5))
    xmax = vc.shape[0]
    for isample, sample in enumerate(sampleinfo):
        ax.plot(vc[:, isample], color=sample[3], label="%s\n$n=%d$" % (sample[2], ws[isample]))
    ax.set_xticks(arange(0, xmax + 1, xmax / 3))
    ax.set_xticklabels(["5'end of transcript", "CDS START", "CDS STOP", "3'end of transcript"], ha="center", weight="bold")
    ax.get_xticklabels()[0].set_horizontalalignment("left")
    ax.get_xticklabels()[-1].set_horizontalalignment("right")
    ax.set_ylabel("Average transcript RPF density\n(intra-gene normalized)")
    legend = ax.legend()
    #pyplot.setp(legend.get_texts(), ha="left", va="bottom")
    #pyplot.setp(ax.get_yticklabels(), weight="bold")
    #fig.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    fig.tight_layout()
    return fig


def main():
    import argparse
    ap = argparse.ArgumentParser(description="")
    ap.add_argument('-m', '--minreads', required=False, default=DEFAULT_MINREADS, type=int)
    ap.add_argument('-o', '--outfile', required=True)
    ap.add_argument('-t', '--txinfofile', required=True)
    ap.add_argument(dest="sampleinfos", metavar="SAMPLEINFO", nargs="*")

    args = ap.parse_args()

    if not args.minreads > 0:
        ap.error("Invalid value for minreads: %d" % args.minreads)

    if not os.access(args.txinfofile, os.R_OK):
        ap.error("Unable to open file: %s" % args.txinfofile)

    #sampleinfo = []
    #for line in open(args.sampleinfofile[0]):
    #    fields = line.rstrip().split("\t")
    #    if len(fields) != 4:
    #        continue

    #    sampleinfo.append(fields)

    sampleinfo = []
    for s in args.sampleinfos:
        fields = s.split(",")
        if len(fields) != 4:
            ap.error("Invalid sampleinfo argument: %s" % s)
        sampleinfo.append(fields)

    vc, ws = get_averaged_transcript_rpf_densities_mutliple_samples(sampleinfo, args.txinfofile, DEFAULT_NGRID, args.minreads)
    fig = plot_average_transcript_rpf_densities(sampleinfo, vc, ws)
    fig.savefig(args.outfile, format="pdf", bbox_inches='tight')

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
