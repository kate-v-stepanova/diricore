#!/usr/bin/env python
#

"""
"""

from itertools import chain

import h5py
import numpy
from numpy.lib.stride_tricks import as_strided
from numpy import convolve, repeat, empty, concatenate, sort, array


def get_cutoff(cvs, min_reads):
    return reduce(lambda x, y: x & y, map(lambda cv: cv.sum(axis=1) >= min_reads, cvs))

#(APU) Normalization of windows (each value divided by the sum and multiplied by the width of the window
def get_normab(cv):
    return ((cv.T * cv.shape[1]) / cv.sum(axis=1).astype(float)).T


def get_smoothed_mean(normab, W):
    return convolve(normab.mean(axis=0), W, mode="same")


def get_means(samples, codongroups, h5mapsfn, min_reads, intersect=False, smooth3nt=True):
    h5maps = h5py.File(h5mapsfn, 'r')

    if smooth3nt:
        W = repeat(1. / 3, 3)
    else:
        W = array([0.0, 1.0, 0.0])

    h5countfiles = {fn: h5py.File(fn, 'r') for fn in set([h5fn for (h5fn, _, _) in samples])}

    # NOTE: here usually a good chunk of memory is needed.. (~700Mb per sample)
    cs = [h5countfiles[h5fn]["counts"][sampleid][:] for (h5fn, sampleid, _) in samples]

    widths = [h5maps["codonmap"][x].attrs["region_width"] for x in set(chain(*codongroups))]
    assert len(set(widths)) == 1
    width = widths[0]

    cvs = [as_strided(c, shape=((c.shape[0] - width) / 1 + 1, width), strides=(1 * c.itemsize, 1 * c.itemsize)) for c in cs]

    means = empty((len(samples), len(codongroups), width))
    passing_cutoff = []

    for icg, cg in enumerate(codongroups):
        m = sort(concatenate([h5maps["codonmap"][c] for c in cg]))

        if intersect:
            w = get_cutoff([cv[m] for cv in cvs], min_reads)
            ws = [w] * len(cvs)
        else:
            ws = [get_cutoff([cv[m]], min_reads) for cv in cvs]

        passing_cutoff.append(ws)

        nabs = [get_normab(cvs[i][m][ws[i]]) for i in xrange(len(cvs))]
        assert all(map(lambda x: numpy.isfinite(x).all(), nabs))
        smabs = [get_smoothed_mean(nab, W) for nab in nabs]
        means[:, icg] = array(smabs)

    return means, passing_cutoff


def read_sampleinfo(fn):
    fh = open(fn)
    samples = []
    for line in fh:
        line = line.rstrip()
        fields = line.split("\t")
        samples.append(tuple(fields))

    return samples
