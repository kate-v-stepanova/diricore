#!/bin/bash

set -e;
set -u;

dataset_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
defaultmeta="$PROJECT_DIR/analysis/input/metadata/rpf_transcript_distribution_sampleinfo.tsv"


project=$dataset_id
minreads=${2};
outname=$dataset_id
metafile=${defaultmeta};

OUTDIR="$PROJECT_DIR/analysis/output/figures/rpf_transcript_distribution";
TXFILE="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.hdf5";
DIRICORE_DIR="/home/e984a/diricore"
python_bin="$DIRICORE_DIR/diricore/bin/plot_rpf_transcript_distribution.py";

mkdir -p ${OUTDIR}

outfile="${OUTDIR}/${project}.rpf_transcript_distribution_plot.${outname}.pdf"

sampinfo=`cat $metafile | while read -a LINE
do
	echo -n "'${LINE[0]},${TXFILE},${LINE[1]},${LINE[2]}' "
done`

eval $(echo python ${python_bin} \
       	-m ${minreads} \
      	-o $outfile \
	"$sampinfo")

echo "Done. File created: $outfile"
