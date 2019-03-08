#!/bin/bash

set -e;
set -u;


BOWTIE_BIN="./programs/bowtie2-2.0.6/bowtie2";
INDIR="./data/input/fastq";
OUTDIR="./data/output/clean";

species=$1
adapter=$2;

RRNA_REF="./staticdata/${species}/rRNAs";
TRNA_REF="./staticdata/${species}/tRNAs";


###
mkdir -p $OUTDIR || true;

ls -1 ${INDIR}/*.fastq.gz | while read fn; do
    bn=$(basename "$fn");
    b=${bn%%.*};

    of="${OUTDIR}/${b}.fastq.gz";
    ca_err="${OUTDIR}/${b}.clipadapt.err";
    ca2_err="${OUTDIR}/${b}.clipadapt_2.err";
    rrna_err="${OUTDIR}/${b}.rrna.err";
    trna_err="${OUTDIR}/${b}.trna.err";

    tmpfile=$(mktemp "/tmp/"${b}".rrna_cleaned.tmp.fastq.gz");
    trap "{ rm -f /tmp/${b}.rrna_cleaned.tmp.fastq.gz;" } EXIT;

    echo "Starting preprocessing of file: ${bn}";
    cat "${fn}" \
    | gzip -dc \
    | cutadapt --quality-base=33 -a "${adapter}" -O 7 -e 0.15 -m 20 -q 5 --untrimmed-output=/dev/null - 2> "${ca_err}" \
    | cutadapt --quality-base=33 -a "GGCATTAACGCGAACTCGGCCTACAATAGT" -a "AAGCGTGTACTCCGAAGAGGATCCAAA" -O 7 -e 0.15 -m 20 - 2> "${ca2_err}" \
    | ${BOWTIE_BIN} --seed 42 -p 1 --local --un-gz "${tmpfile}" "${RRNA_REF}" - \
    > /dev/null 2> "${rrna_err}";

    cat "${tmpfile}" \
    | gzip -dc \
    | ${BOWTIE_BIN} --seed 42 --local --un-gz "${of}" "${TRNA_REF}" - \
    > /dev/null 2> "${trna_err}";

    rm ${tmpfile};
done

echo "Done with preprocessing";
###
