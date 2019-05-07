#!/bin/bash

set -e;
set -u;

dataset_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
metafile="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv"

project=$dataset_id
minreads=${2};
genome=$3

OUTDIR="$PROJECT_DIR/analysis/output/figures/rpf_transcript_distribution";
TX_ALL="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.all.hdf5";
TX_ALL_DEDUP="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.all.dedup.hdf5";
TX_HQ="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.hq.hdf5";
TX_HQ_DEDUP="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.hq.dedup.hdf5";
TX_INFO_FILE="$BASE_DIR/static/$genome/transcript_data.hdf5"

DIRICORE_DIR="/home/e984a/diricore"
python_bin="$DIRICORE_DIR/diricore/bin/plot_rpf_transcript_distribution.py";

mkdir -p ${OUTDIR}

outfile="${OUTDIR}/all.${project}.rpf_transcript_distribution_plot.${dataset_id}.pdf"

sampinfo=`cat $metafile | while read -a LINE
do
	echo -n "'${LINE[0]},${TX_ALL},${LINE[1]},${LINE[2]}' "
done`

eval $(echo python ${python_bin} \
       	-m ${minreads} \
      	-o $outfile \
        -t $TX_INFO_FILE \
	"$sampinfo")

echo "Done. File created: $outfile"

outfile="${OUTDIR}/all.dedup.${project}.rpf_transcript_distribution_plot.${dataset_id}.pdf"

sampinfo=`cat $metafile | while read -a LINE
do
	echo -n "'${LINE[0]},${TX_ALL_DEDUP},${LINE[1]},${LINE[2]}' "
done`

eval $(echo python ${python_bin} \
       	-m ${minreads} \
      	-o $outfile \
        -t $TX_INFO_FILE \
	"$sampinfo")

echo "Done. File created: $outfile"


outfile="${OUTDIR}/hq.${project}.rpf_transcript_distribution_plot.${dataset_id}.pdf"

sampinfo=`cat $metafile | while read -a LINE
do
	echo -n "'${LINE[0]},${TX_HQ},${LINE[1]},${LINE[2]}' "
done`

eval $(echo python ${python_bin} \
       	-m ${minreads} \
      	-o $outfile \
        -t $TX_INFO_FILE \
	"$sampinfo")

echo "Done. File created: $outfile"


outfile="${OUTDIR}/hq.dedup.${project}.rpf_transcript_distribution_plot.${dataset_id}.pdf"

sampinfo=`cat $metafile | while read -a LINE
do
	echo -n "'${LINE[0]},${TX_HQ_DEDUP},${LINE[1]},${LINE[2]}' "
done`

eval $(echo python ${python_bin} \
       	-m ${minreads} \
      	-o $outfile \
        -t $TX_INFO_FILE \
	"$sampinfo")

echo "Done. File created: $outfile"
