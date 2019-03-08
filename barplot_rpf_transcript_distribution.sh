#!/bin/bash

. ./diricore_virtualenv/bin/activate

set -e;
set -u;



defaultmeta="./data/input/metadata/rpf_transcript_distribution_sampleinfo.tsv"


project=$1
minreads=${2:-100};
outname=${3:-1};
metafile=${4:-$defaultmeta};

OUTDIR="./data/output/rpf_5p_density/";
TXFILE="./data/output/rpf_5p_density/${project}.txcoord_counts.hdf5";
python_bin="./diricore/bin/barplot_rpf_transcript_distribution.py";


mkdir ${OUTDIR} || true;



outfile="${OUTDIR}/${project}.rpf_in_regs.${outname}"
sampinfo=`cat $metafile | while read -a LINE
do
	echo -n "'${LINE[0]},${TXFILE},${LINE[1]},${LINE[2]}' "
done`


eval $(echo python ${python_bin} \
       	-m ${minreads} \
       	-o $outfile \
	"$sampinfo")








