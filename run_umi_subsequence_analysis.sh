#!/bin/bash

set -e;
set -u;


OUTDIR="./data/output/subsequence_data";
INDIR="./data/output/tophat_out";
PLOTDIR="./data/output/figures";
SAMPLENAME_FILE="./data/input/metadata/subsequence_samplenames.tsv";
CONTRAST_FILE="./data/input/metadata/subsequence_contrasts.tsv";
RPF_CONTRASTS="./data/input/metadata/rpf_density_contrasts.tsv";
RPF_SAMPLENAME="./data/input/metadata/rpf_density_samplenames.tsv"

species=$1;
projectname=$2;
minreads=$3;

INDEXDATAFN="./staticdata/${species}/subseq_index_data.pkl.gz";
of="${OUTDIR}/${projectname}.subsequence_data.frame0.hdf5";

if [ ! -f ${CONTRAST_FILE} ]; then
    cut -f1,2 ${RPF_CONTRASTS} > ${CONTRAST_FILE}
fi

if [ ! -f ${SAMPLENAME_FILE} ]; then
    cp ${RPF_SAMPLENAME} ${SAMPLENAME_FILE}
fi

###
rm -rf data/output/subsequence_data/*
mkdir -p ${OUTDIR} || true;
mkdir -p "${PLOTDIR}/subsequence_shift_plots/" || true;

ls -1 ${INDIR}/*/accepted_hits.hqmapped_dedup.bam | sort -V | while read bamfn; do
    b=$(basename $(dirname "$bamfn"));

    ./diricore/bin/extract_subsequences.py \
        -v \
        run \
        -o "${of}" \
        -f 0 \
        "${INDEXDATAFN}" \
        "${b},${bamfn}" \
    ;
done

# create subsequence shift plots
./diricore/bin/plot_subsequence_shifts.py \
    -o ${PLOTDIR}/subsequence_shift_plots/${projectname}.m${minreads}. \
    -m $minreads \
    --sample-names ${SAMPLENAME_FILE} \
    --contrasts ${CONTRAST_FILE} \
    "${of}"
###
