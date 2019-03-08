#!/bin/bash


set -e;
set -u;


export BOWTIE_BIN="./programs/bowtie2-2.0.6/bowtie2";
export INDIR="./data/input/fastq/";
export OUTDIR="./data/output/clean/";
echo "Settings:";
export species=$1; echo -e "\tspecies:$species";
export cores=${2:-1}; echo -e "\tparallel_processes: $cores";

export RRNA_REF="./staticdata/${species}/rRNAs";
export TRNA_REF="./staticdata/${species}/tRNAs";


###
mkdir -p $OUTDIR || true;

run_prep() {
    fn=$1
    bn=$(basename "$fn");
    b=${bn%%.*};

    of="${OUTDIR}/${b}.fastq.gz";
    rrna_err="${OUTDIR}/${b}.rrna.err";
    trna_err="${OUTDIR}/${b}.trna.err";

    tmpfile=$(tempfile -d "/tmp/" -s ".${b}.rrna_cleaned.tmp.fastq.gz");

    echo "Starting preprocessing of file: ${bn}";
    cat "${fn}" \
    | gzip -dc \
    | ${BOWTIE_BIN} --seed 42 -p 1 --local --un-gz "${tmpfile}" "${RRNA_REF}" - \
    > /dev/null 2> "${rrna_err}";

    cat "${tmpfile}" \
    | gzip -dc \
    | ${BOWTIE_BIN} --seed 42 --local --un-gz "${of}" "${TRNA_REF}" - \
    > /dev/null 2> "${trna_err}";

    rm ${tmpfile};
}
export -f run_prep
echo "Preprocessing files in parallel..."
ls ${INDIR}/*.fastq.gz | parallel run_prep

echo "Done with preprocessing";
###
