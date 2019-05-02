#!/bin/bash

set -e;
set -u;

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
project_id=$1
project_dir="$BASE_DIR/$project_id"

INDIR="$project_dir/analysis/input/fastq";
OUTDIR="$project_dir/analysis/output/clean";

species=$2
# adapter=$3;

RRNA_REF="$BASE_DIR/static/${species}/rRNAs";
TRNA_REF="$BASE_DIR/static/${species}/tRNAs";
BOWTIE_BIN="/home/e984a/diricore/programs/bowtie2-2.0.6/bowtie2"

###
mkdir -p $OUTDIR || true;

ls -1 ${INDIR}/*.fastq.gz | while read fn; do
    echo "Starting preprocessing of file: ${fn}";

    bn=$(basename "$fn");
    b=${bn%%.*};

    of="${OUTDIR}/${b}.fastq.gz";
    ca_err="${OUTDIR}/${b}.clipadapt.err";
    ca2_err="${OUTDIR}/${b}.clipadapt_2.err";
    rrna_err="${OUTDIR}/${b}.rrna.err";
    trna_err="${OUTDIR}/${b}.trna.err";

    tmpfile="${INDIR}/${b}.rrna_cleaned.tmp.fastq.gz";
    rm -f $tmpfile;
    cp ${INDIR}/${b}.fastq.gz $tmpfile
    trap "{ rm -f ${tmpfile}; }" EXIT;
    cat "${fn}" | gzip -dc | ${BOWTIE_BIN} --seed 42 -p 1 --local --un-gz "${tmpfile}" "${RRNA_REF}" -  # > /dev/null
    cat "${tmpfile}" | gzip -dc  | ${BOWTIE_BIN} --seed 42 --local --un-gz "${of}" "${TRNA_REF}" - # > /dev/null # 2> "${trna_err}";
    rm -f ${tmpfile};
done

echo "Done with preprocessing";
###
