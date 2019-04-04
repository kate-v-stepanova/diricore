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

plots_only=0
if [[ $# -eq 4 ]]; then
    plots=$4
    if [[ $plots -eq "plots_only" ]]; then
        plots_only=1
    fi
fi

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"

OUTDIR="$PROJECT_DIR/analysis/output/subsequence_data";
INDIR="$PROJECT_DIR/analysis/output/tophat_out";
PLOTDIR="$PROJECT_DIR/analysis/output/figures";
CONTRASTS="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";
SAMPLENAMES="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv"
DIRICORE_DIR="/home/e984a/diricore"
INDEXDATAFN="$DIRICORE_DIR/staticdata/${species}/subseq_index_data.pkl.gz";

frame_file="${OUTDIR}/${projectname}.subsequence_data.frame0.hdf5"
hq_frame_file="${OUTDIR}/${projectname}.subsequence_data.frame0.hq.hdf5"


#if [ ! -f ${CONTRAST_FILE} ]; then
#    cut -f1,2 ${RPF_CONTRASTS} > ${CONTRAST_FILE}
#fi

#if [ ! -f ${SAMPLENAME_FILE} ]; then
#    cp ${RPF_SAMPLENAME} ${SAMPLENAME_FILE}
#fi

###
mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}/subsequence_shift_plots/"

if [[ $plots_only -eq 0 ]]; then
    # rm -f $OUTDIR/*
    rm -f $frame_file
    rm -f $hq_frame_file
    echo "Extracting HQ subsequences"
    ls -1 ${INDIR}/*/accepted_hits.hqmapped_dedup.bam | sort -V | while read bamfn; do
        b=$(basename $(dirname "$bamfn"));
        b=${b%%.*};
        b=${b#"dem_"};
        b=${b%"_umi_extracted"};
       $DIRICORE_DIR/diricore/bin/extract_subsequences.py \
           -v \
           run \
           -o $hq_frame_file \
           -f 0 \
           ${INDEXDATAFN} \
           ${b},${bamfn}
    done
    echo "Done. Generated file: ${hq_frame_file}"

    echo "Extracting all subsequences"
    ls -1 ${INDIR}/*/accepted_hits.bam | sort -V | while read bamfn; do
        b=$(basename $(dirname "$bamfn"));
        b=${b%%.*};
        b=${b#"dem_"};
        b=${b%"_umi_extracted"};
        $DIRICORE_DIR/diricore/bin/extract_subsequences.py \
           -v \
           run \
           -o ${frame_file} \
           -f 0 \
           ${INDEXDATAFN} \
           ${b},${bamfn} \
       ;
    done
    echo "Done. Generated file: ${frame_file}"
fi

echo "Creating HQ plots"
# create subsequence shift plots
$DIRICORE_DIR/diricore/bin/plot_subsequence_shifts.py \
    -o ${PLOTDIR}/subsequence_shift_plots/${projectname}.hq.m${minreads}. \
    -m $minreads \
    --sample-names ${SAMPLENAMES} \
    --contrasts ${CONTRASTS} \
    ${hq_frame_file}
echo "Created plots in ${PLOTDIR}/subsequence_shift_plots"
###

echo "Creating ALL plots"
# create subsequence shift plots
$DIRICORE_DIR/diricore/bin/plot_subsequence_shifts.py \
    -o ${PLOTDIR}/subsequence_shift_plots/${projectname}.all.m${minreads}. \
    -m $minreads \
    --sample-names ${SAMPLENAMES} \
    --contrasts ${CONTRASTS} \
    ${frame_file}
echo "Created plots in ${PLOTDIR}/subsequence_shift_plots"

