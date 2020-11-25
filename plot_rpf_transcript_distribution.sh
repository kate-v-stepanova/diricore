#!/bin/bash

set -e;
set -u;

dataset_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
metafile="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv"

project=$dataset_id
minreads=$3;
genome=$2

OUTDIR="$PROJECT_DIR/analysis/output/figures/rpf_transcript_distribution";
TX_ALL="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.all.${minreads}.hdf5";
TX_ALL_DEDUP="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.all.dedup.${minreads}.hdf5";
TX_HQ="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.hq.${minreads}.hdf5";
TX_HQ_DEDUP="$PROJECT_DIR/analysis/output/rpf_5p_density/${project}.txcoord_counts.hq.dedup.${minreads}.hdf5";
TX_INFO_FILE="$BASE_DIR/static/$genome/transcript_data.hdf5"

DIRICORE_DIR="/icgc/dkfzlsdf/analysis/OE0532/software/diricore"
python_bin="$DIRICORE_DIR/diricore/bin/plot_rpf_transcript_distribution.py";


# $4 can be: all, all_unique, hq, hq_unique
infile=$TX_HQ_DEDUP
file_type="hq_unique"
if [[ $# -ge 4 ]]; then
    file_type=$4
    if [[ $4 == "hq" ]]; then
        infile=$TX_HQ
    elif [[ $4 == "all_unique" ]]; then
        infile=$TX_ALL_DEDUP
    elif [[ $4 == "all" ]]; then
        infile=$TX_ALL
    else 
        file_type="hq_unique"
    fi
fi
p_id=${dataset_id/"/"/"_"}
outfile="${OUTDIR}/${file_type}.${p_id}.${minreads}.rpf_transcript_distribution_plot.${p_id}.pdf"

if [[ $# -ge 5 ]]; then
   subset=$5
   metafile="$PROJECT_DIR/analysis/input/metadata/${subset}_rpf_density_samplenames.tsv"

   outfile="${OUTDIR}/${file_type}.${subset}.${p_id}.${minreads}.rpf_transcript_distribution_plot.${p_id}.pdf"
fi

mkdir -p ${OUTDIR}


if [[ -f $infile ]]; then
  echo "Processing $infile"
  sampinfo=`cat $metafile | while read -a LINE
  do
	echo -n "'${LINE[0]},${infile},${LINE[1]},${LINE[2]}' "
  done`

  eval $(echo python ${python_bin} \
       	-m ${minreads} \
      	-o $outfile \
        -t $TX_INFO_FILE \
	"$sampinfo")

  echo "Done. File created: $outfile"
else
  echo "$infile does not exist."
fi
