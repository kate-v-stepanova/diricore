#!/bin/bash

set -e;
set -u;

dataset_id=$1

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
INDIR="$PROJECT_DIR/analysis/output/tophat_out"
OUTDIR="$PROJECT_DIR/analysis/output/tophat_out"

## HQ mapped reads
echo "HQ mapped reads"
# index
for d in $(ls ${INDIR}); do
    f="${INDIR}/${d}/accepted_hits.hqmapped.bam";
    outfile="${OUTDIR}/${d}/accepted_hits.hqmapped.bam.bai";
    echo "Indexing ${f}";
    samtools index ${f}
done;

for d in $(ls ${INDIR}); do
    f="${INDIR}/${d}/accepted_hits.hqmapped.bam";
    echo "Deduplicating ${f}";
    fb=$(basename ${f});
    fb=${fb%.*};
    umi_tools dedup --log2stderr -L ${OUTDIR}/${d}/${fb}_umidedup_log.txt -I ${f} -S ${OUTDIR}/${d}/${fb}_dedup.bam --output-stats=${OUTDIR}/${d}/${fb}_dedup.log
done;


## ALL reads
echo "ALL reads"
for d1 in $(ls ${INDIR}); do
    f="${INDIR}/${d1}/accepted_hits.bam";
    outfile="${OUTDIR}/${d1}/accepted_hits.bam.bai";
    echo "Indexing ${f}";
    samtools index ${f}
done;

for d1 in $(ls ${INDIR}); do
   if [[ "$d1" != "logs" ]]; then 
        echo $d1;
	f1="${INDIR}/${d1}/accepted_hits.bam";
        echo "Deduplicating ${f1}";
        fb1=$(basename ${f1});
        fb1=${fb1%.*};
        umi_tools dedup --log2stderr -L ${OUTDIR}/${d1}/${fb1}_umidedup_log.txt -I ${f1} -S ${OUTDIR}/${d1}/${fb1}_dedup.bam --output-stats=${OUTDIR}/${d1}/${fb1}_dedup.log
    fi
done;
