#!/usr/bin/env python
#

import os
from itertools import izip
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')

import h5py
from numpy import array, zeros, empty, arange, diff, concatenate, split
from matplotlib import pyplot


DEFAULT_TXINFOFN = "./staticdata/human/transcript_data.hdf5"
DEFAULT_MINREADS = 100

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


def get_transcript_rpf_density(A, infofn, chosentids,minreads):

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

    counts = zeros((len(chosentids),5), dtype="uint32")
    counts_norm = zeros((len(chosentids),5), dtype="float32")
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

        if t < minreads:
            continue
	
	idxs = [cdsstart -30, cdsstart+30, cdsend -30, cdsend +30]	

	if cdsstart > 50 and ((l-cdsend) > 50):
		cs = split(a, idxs)
		cs = array(map(sum, cs))
		counts[itid] = cs * 100.0 / cs.sum()
		cs = cs * 1000.0 / diff([0] + idxs + [l])
		counts_norm[itid] = cs * 100.0 / cs.sum()

		
    return counts, counts_norm


def get_averaged_transcript_rpf_densities_mutliple_samples(sampleinfo, infofn, minreads,outf):
    chosentids = select_representative_isoforms(infofn)
    res_counts = zeros([len(sampleinfo), 5],dtype='float32')
    res_counts_norm = zeros([len(sampleinfo), 5],dtype='float32')
    for isample, sample in enumerate(sampleinfo):
        h5data = h5py.File(sample[1], 'r')
        A = h5data["counts"][sample[0]][:]
        counts, counts_norm = get_transcript_rpf_density(A, infofn, chosentids,minreads)
	counts = counts[counts.sum(axis=1) > 0,]
	counts_norm = counts_norm[counts_norm.sum(axis=1) > 0]
	res_counts[isample,] = counts.mean(axis=0)
	res_counts_norm[isample,] = counts_norm.mean(axis=0)
    outfo = open(outf + '.txt', 'w')
    outfon = open(outf + 'norm.txt', 'w')
    outfo.write('Sample\tUTR5\tStart\tCDS\tStop\tUTR3\n')
    outfon.write('Sample\tUTR5\tStart\tCDS\tStop\tUTR3\n')
    for isample, sample in enumerate(sampleinfo):
	outfo.write(sample[0] + '\t' + '\t'.join(map(str, res_counts[isample,])) + '\n') 
	outfon.write(sample[0] + '\t' + '\t'.join(map(str, res_counts_norm[isample,])) + '\n')    
    outfo.close()
    outfon.close()
    return res_counts, res_counts_norm

def main():
    import argparse
    ap = argparse.ArgumentParser(description="")
    ap.add_argument('-m', '--minreads', required=False, default=DEFAULT_MINREADS, type=int)
    ap.add_argument('-o', '--outfile', required=True)
    ap.add_argument('-t', '--txinfofile', required=False, default=DEFAULT_TXINFOFN)
    #ap.add_argument(dest="sampleinfofile", metavar="SAMPLEINFOFILE", nargs=1)
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


    counts, countsn = get_averaged_transcript_rpf_densities_mutliple_samples(sampleinfo, args.txinfofile,args.minreads,args.outfile)


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
