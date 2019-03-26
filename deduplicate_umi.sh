#!/bin/bash

set -e;
set -u;

dataset_id=$1

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
INDIR="$PROJECT_DIR/analysis/output/tophat_out"
OUTDIR="$PROJECT_DIR/analysis/output/tophat_out"

# index
for d in $(ls ${INDIR}); do
    f="${INDIR}/${d}/accepted_hits.hqmapped.bam";
    outfile="${OUTDIR}/${d}/accepted_hits.hqmapped.bam.bai";
    echo "Indexing ${f}";
    samtools index ${f};
done;

for d in $(ls ${INDIR}); do
    f="${INDIR}/${d}/accepted_hits.hqmapped.bam";
    echo "Deduplicating ${f}";
    fb=$(basename ${f});
    fb=${fb%.*};
    umi_tools dedup --log2stderr -L umidedup_log.txt -I ${f} -S ${OUTDIR}/${d}/${fb}_dedup.bam --output-stats=${OUTDIR}/${d}/${fb}_dedup.log
done;
