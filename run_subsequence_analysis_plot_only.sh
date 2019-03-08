#!/bin/bash

set -e;
set -u;


#This version of the script disables the mapping step, so it can be used with previously analyzed data without doing that step again.
#IMPORTANT: project name and sample name must be exactly as it was used to run the complete script, to match the names of the previously generated H5py file


OUTDIR="./data/output/subsequence_data/";
INDIR="./data/output/tophat_out/";
PLOTDIR="./data/output/figures/";
SAMPLENAME_FILE="./data/input/metadata/subsequence_samplenames.tsv";
CONTRAST_FILE="./data/input/metadata/subsequence_contrasts.tsv";

species=$1;
projectname=$2;
minreads=$3;

INDEXDATAFN="./staticdata/${species}/subseq_index_data.pkl.gz";
of="${OUTDIR}/${projectname}.subsequence_data.frame0.hdf5";


###
mkdir ${OUTDIR} || true;
mkdir -p "${PLOTDIR}/subsequence_shift_plots/" || true;


# create subsequence shift plots
./diricore/bin/plot_subsequence_shifts.py \
    -o ${PLOTDIR}/subsequence_shift_plots/${projectname}.m${minreads}. \
    -m $minreads \
    --y-limits " -1,1" \
    --sample-names ${SAMPLENAME_FILE} \
    --contrasts ${CONTRAST_FILE} \
    "${of}"
###
