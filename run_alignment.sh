#! /bin/bash

#. ./diricore_virtualenv/bin/activate


set -e;
set -u;


BOWTIE_PATH="./programs/bowtie2-2.0.6";
TOPHAT_BIN="./programs/tophat-2.0.7/tophat2";

INDIR="./data/output/clean";
OUTDIR="./data/output/tophat_out";

species=$1;

REF="./staticdata/${species}/genome";
GTF="./staticdata/${species}/transcripts.gff";
TIDX="./staticdata/${species}/transcripts";


###
ls -1 ${INDIR}/*.fastq.gz | while read fn; do
    bn=$(basename "$fn");
    b=${bn%%.*};

    od="${OUTDIR}/${b}";
    mkdir -p "${od}" 2> /dev/null || true;
    thout="${OUTDIR}/${b}/tophat.out"
    therr="${OUTDIR}/${b}/tophat.err"

    echo "Starting alignment of file: $bn";

    PATH="${BOWTIE_PATH}:$PATH" \
        ${TOPHAT_BIN} --seed 42 -n 2 -m 1 \
        --no-novel-juncs --no-novel-indels --no-coverage-search \
        --segment-length 25 \
        --transcriptome-index "${TIDX}" -G "${GTF}" \
        -o "${od}" \
        -p 1 "${REF}" "${fn}" \
        > "${thout}" 2> "${therr}";
done

echo "Done with alignments; now filtering out non-primary/low-quality alignments";

# isolate "hqmapped" reads
ls -1 ${OUTDIR}/*/accepted_hits.bam | while read fn; do
    dn=$(dirname "$fn");
    of="${dn}/accepted_hits.hqmapped.bam";

    cat \
        <(samtools view -H "${fn}") \
        <(cat "${fn}" | samtools view -q10 -F260 -) \
    | samtools view -bS - \
    > "${of}";
done

echo "Done filtering";
###
