#!/bin/bash

set -e;
set -u;

metafile="./data/input/metadata/rpf_transcript_distribution_sampleinfo.tsv"

#project=$1
minreads=${1:-100};
#outname=${3:-1};

OUTDIR="./data/output/rpf_5p_density";
TXFILE="./data/output/subsequence_data/miseq_data.subsequence_data.hdf5"
#TXFILE=$(ls ${OUTDIR}/*.hdf5)
python_bin="/home/e984a/diricore/diricore/bin/barplot_rpf_transcript_distribution.py";
project=$(basename ${TXFILE})
project=${project%.*.*}
outname=$project

mkdir -p ${OUTDIR}

outfile="${OUTDIR}/${project}.rpf_in_regs.${outname}.txt"
echo $outfile
sampinfo=`cat $metafile | while read -a LINE
do
	echo -n "'${LINE[0]},${TXFILE},${LINE[1]},${LINE[2]}' "
done`
echo $sampinfo
echo $minreads

eval $(echo python ${python_bin} \
       	-m ${minreads} \
	-t ${TXFILE} \
       	-o $outfile \
	"$sampinfo")
