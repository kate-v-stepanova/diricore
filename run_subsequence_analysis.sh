#!/bin/bash


set -e;
set -u;

OUTDIR="./data/output/subsequence_data";
INDIR="./data/output/tophat_out";
PLOTDIR="./data/output/figures";
SAMPLENAME_FILE="./data/input/metadata/rpf_density_samplenames.tsv";
CONTRAST_FILE="./data/input/metadata/rpf_density_contrasts.tsv";

species=$1;
projectname=$2;
minreads=$3;

INDEXDATAFN="./staticdata/${species}/subseq_index_data.pkl.gz";
of="${OUTDIR}/${projectname}.subsequence_data.hdf5";
bin_extract="/home/e984a/diricore/diricore/bin/extract_subsequences.py"
bin_plot="/home/e984a/diricore/diricore/bin/plot_subsequence_shifts.py"

###
# rm -rf ${OUTDIR}/*
mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}/subsequence_shift_plots"

ls -1 ${INDIR}/*/accepted_hits.hqmapped.bam | sort -V | while read bamfn; do
    b=$(basename $(dirname "$bamfn"));
    ${bin_extract} -v run -o "${of}" -f 0 "${INDEXDATAFN}" "${b},${bamfn}";
done

# create subsequence shift plots
    $bin_plot \
    -o ${PLOTDIR}/subsequence_shift_plots/${projectname}.m${minreads}. \
    -m $minreads \
    --sample-names ${SAMPLENAME_FILE} \
    --contrasts ${CONTRAST_FILE} \
    "${of}"
###
