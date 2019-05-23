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
TX_ALL="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.all.${minreads}.hdf5";
TX_ALL_DEDUP="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.all.dedup.${minreads}.hdf5";
TX_HQ="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.hq.${minreads}.hdf5";
TX_HQ_DEDUP="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.hq.dedup.${minreads}.hdf5";
TX_INFO_FILE="$BASE_DIR/static/$genome/transcript_data.hdf5"

DIRICORE_DIR="/home/e984a/diricore"
python_bin="$DIRICORE_DIR/diricore/bin/plot_rpf_transcript_distribution.py";

mkdir -p ${OUTDIR}

if [[ -f $TX_ALL ]]; then
  outfile="${OUTDIR}/all.${project}.${minreads}.rpf_transcript_distribution_plot.${dataset_id}.pdf"
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
else
  echo "$TX_ALL does not exist. Skipping"
fi

if [[ -f $TX_ALL_DEDUP ]]; then
  outfile="${OUTDIR}/all.dedup.${project}.${minreads}.rpf_transcript_distribution_plot.${dataset_id}.pdf"
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
else
  echo "$TX_ALL_DEDUP does not exist. Skipping"
fi

if [[ -f $TX_HQ ]]; then
  outfile="${OUTDIR}/hq.${project}.${minreads}.rpf_transcript_distribution_plot.${dataset_id}.pdf"
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
else
  echo "$TX_HQ does not exist. Skipping"
fi

if [[ -f $TX_HQ_DEDUP ]]; then
  outfile="${OUTDIR}/hq.dedup.${project}.${minreads}.rpf_transcript_distribution_plot.${dataset_id}.pdf"
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
else
  echo "$TX_HQ_DEDUP does not exist. Skipping"
fi
