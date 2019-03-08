#!/bin/bash

set -e;
set -u;


export OUTDIR="./data/output/subsequence_data/";
export INDIR="./data/output/tophat_out/";
export PLOTDIR="./data/output/figures/";
export SAMPLENAME_FILE="./data/input/metadata/subsequence_samplenames.tsv";
export CONTRAST_FILE="./data/input/metadata/subsequence_contrasts.tsv";

export species=$1;
export projectname=$2;
export minreads=$3;

export INDEXDATAFN="./staticdata/${species}/subseq_index_data.pkl.gz";
export of="${OUTDIR}/${projectname}.subsequence_data.frame0.hdf5";


###
mkdir ${OUTDIR} || true;
mkdir -p "${PLOTDIR}/subsequence_shift_plots/" || true;

extract_sub(){
	bamfn=$1;
	b=$(basename $(dirname "$bamfn"));
	./diricore/bin/extract_subsequences.py \
       	-v run \
        -o "${of}" \
        -f 0 \
        "${INDEXDATAFN}" \
        "${b},${bamfn}" \
	;
}

export -f extract_sub
ls -1 ${INDIR}/*/accepted_hits.hqmapped.bam | sort -V | parallel extract_sub


# create subsequence shift plots
./diricore/bin/plot_subsequence_shifts.py \
    -o ${PLOTDIR}/subsequence_shift_plots/${projectname}.m${minreads}. \
    -m $minreads \
    --sample-names ${SAMPLENAME_FILE} \
    --contrasts ${CONTRAST_FILE} \
    "${of}"
###
