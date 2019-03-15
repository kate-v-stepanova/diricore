#!/bin/bash

set -e;
set -u;


defaultmeta="./data/input/metadata/rpf_transcript_distribution_sampleinfo.tsv"

project=$1
minreads=${2};
outname=${3};
metafile=${defaultmeta};

OUTDIR="./data/output/figures/rpf_transcript_distribution";
TXFILE="./data/output/rpf_5p_density/${project}.txcoord_counts.hdf5";
python_bin="./diricore/bin/plot_rpf_transcript_distribution.py";

mkdir -p ${OUTDIR} || true;

outfile="${OUTDIR}/${project}.rpf_transcript_distribution_plot.${outname}.pdf"

sampinfo=`cat $metafile | while read -a LINE
do
	echo -n "'${LINE[0]},${TXFILE},${LINE[1]},${LINE[2]}' "
done`

eval $(echo python ${python_bin} \
       	-m ${minreads} \
      	-o $outfile \
	"$sampinfo")
