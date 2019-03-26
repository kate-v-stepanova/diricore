#!/bin/bash

set -e;
set -u;

dataset_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"

OUTDIR="$PROJECT_DIR/analysis/output/subsequence_data";
INDIR="$PROJECT_DIR/analysis/output/tophat_out";
PLOTDIR="$PROJECT_DIR/analysis/output/figures";
SAMPLENAME_FILE="$PROJECT_DIR/analysis/input/metadata/subsequence_samplenames.tsv";
CONTRAST_FILE="$PROJECT_DIR/analysis/input/metadata/subsequence_contrasts.tsv";
RPF_CONTRASTS="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";
RPF_SAMPLENAME="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv"

species=$2;
projectname=$dataset_id

if [[ $# -eq 3 ]]; then
  minreads=$3;
else
  minreads=15
fi

DIRICORE_DIR="/home/e984a/diricore"
INDEXDATAFN="$DIRICORE_DIR/staticdata/${species}/subseq_index_data.pkl.gz";

if [ ! -f ${CONTRAST_FILE} ]; then
    cut -f1,2 ${RPF_CONTRASTS} > ${CONTRAST_FILE}
fi

if [ ! -f ${SAMPLENAME_FILE} ]; then
    cp ${RPF_SAMPLENAME} ${SAMPLENAME_FILE}
fi

###
rm -rf $PROJECT_DIR/analysis/output/subsequence_data/*
mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}/subsequence_shift_plots/"

echo "Extracting subsequences"
of="${OUTDIR}/${projectname}.subsequence_data.frame0.hdf5";
ls -1 ${INDIR}/*/accepted_hits.hqmapped_dedup.bam | sort -V | while read bamfn; do
    b=$(basename $(dirname "$bamfn"));
    b=${b%%.*};
    b=${b#"dem_"};
    b=${b%"_umi_extracted"};
   $DIRICORE_DIR/diricore/bin/extract_subsequences.py \
       -v \
       run \
       -o "${of}" \
       -f 0 \
       "${INDEXDATAFN}" \
       "${b},${bamfn}" \
   ;
done
echo "Done. Generated file: ${of}"

echo "Creating plots"
# create subsequence shift plots
$DIRICORE_DIR/diricore/bin/plot_subsequence_shifts.py \
    -o ${PLOTDIR}/subsequence_shift_plots/${projectname}.m${minreads}. \
    -m $minreads \
    --sample-names ${SAMPLENAME_FILE} \
    --contrasts ${CONTRAST_FILE} \
    "${of}"
echo "Created plots in ${PLOTDIR}/subsequence_shift_plots"
###
