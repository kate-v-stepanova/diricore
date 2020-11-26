#!/bin/bash

set -e;
set -u;

dataset_id=$1
species=$2;
projectname=$dataset_id

if [[ $# -ge 3 ]]; then
  minreads=$3;
else
  minreads=15
fi

bam_type="hq_unique"
# can be: hq, hq_unique, all, all_unique
if [[ $# -ge 4 ]]; then
   bam_type=$4
fi
echo "$bam_type"

plots_only=0
if [[ $# -ge 5 ]]; then
    plots=$5
    if [[ $plots -eq "plots_only" ]]; then
        plots_only=1
    fi
fi
echo "$plots_only"

y_limits=1
if [[ $# -eq 6 ]]; then
    y_limits=$6
fi

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"

OUTDIR="$PROJECT_DIR/analysis/output/subsequence_data";
#INDIR="$PROJECT_DIR/analysis/output/tophat_out";
INDIR="$PROJECT_DIR/analysis/output/alignments/toGenome"
PLOTDIR="$PROJECT_DIR/analysis/output/figures";
CONTRASTS="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";
SAMPLENAMES="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv"
DIRICORE_DIR="/home/e984a/diricore"
INDEXDATAFN="$BASE_DIR/static/${species}/subseq_index_data.pkl.gz";

all_frame_file="${OUTDIR}/${projectname}.subsequence_data.all.${minreads}.hdf5"
all_dedup_frame_file="${OUTDIR}/${projectname}.subsequence_data.all.dedup.${minreads}.hdf5"
hq_frame_file="${OUTDIR}/${projectname}.subsequence_data.hq.${minreads}.hdf5"
hq_dedup_frame_file="${OUTDIR}/${projectname}.subsequence_data.hq.dedup.${minreads}.hdf5"

if [ $bam_type == "all" ]; then
    frame_file=$all_frame_file
    bam_pattern="_toGenome.bam"
elif [ $bam_type == "all_unique" ]; then
    frame_file=$all_dedup_frame_file
    bam_pattern="_toGenome_dedup.bam"
elif [ $bam_type == "hq" ]; then
    frame_file=$hq_frame_file
    bam_pattern="_toGenome.hqmapped.bam"
else
    frame_file=$hq_dedup_frame_file
    bam_pattern="_toGenome.hqmapped_dedup.bam"
fi

###
mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}/subsequence_shift_plots/"

if [[ -f $frame_file ]]; then
    echo "File exists: $frame_file"
    echo "Delete file or use option 'plots_only'"
    exit
fi

if [ $plots_only -eq 0 ]; then
    # rm -f $OUTDIR/*
    rm -f $frame_file;
    echo "Extracting subsequences ($bam_type)";
    for bamfn in $(ls $INDIR/*$bam_pattern); do
        b=$(basename $bamfn);
        b=${b%"$bam_pattern"};
        $DIRICORE_DIR/diricore/bin/extract_subsequences.py \
           -v \
           run \
           -o $frame_file \
           -f 0 \
           ${INDEXDATAFN} \
           ${b},${bamfn}
    done;
    echo "Done. Generated file: ${frame_file}";
fi


echo "Creating plots"
if [[ -f $frame_file ]]; then
    # create subsequence shift plots
    mkdir -p ${PLOTDIR}/subsequence_shift_plots/$bam_type
    $DIRICORE_DIR/diricore/bin/plot_subsequence_shifts.py \
    -o ${PLOTDIR}/subsequence_shift_plots/$bam_type/${projectname}.$bam_type.m${minreads}. \
    -m $minreads \
    --sample-names ${SAMPLENAMES} \
    --contrasts ${CONTRASTS} \
    --y-limits ${y_limits} \
    ${frame_file}
    echo "Created plots in ${PLOTDIR}/subsequence_shift_plots"
else
    echo "ERROR: file $frame_file not found!! Skipping"
fi
###

