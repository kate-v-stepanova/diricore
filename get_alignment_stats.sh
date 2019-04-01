#!/bin/bash

dataset_id=$1

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
INDIR="$PROJECT_DIR/analysis/output/tophat_out"
OUTDIR="$PROJECT_DIR/analysis/output/alignment_stats"

mkdir -p ${OUTDIR}

# Get useful alignment number (HQ_aligns) before dedup
rm -f ${OUTDIR}/alignment_hq_stats.txt
touch ${OUTDIR}/alignment_hq_stats.txt
for inf in $(ls ${INDIR}/*/accepted_hits.hqmapped.bam); do
 c=$(samtools view -c $inf)
 echo -e "${inf}\t${c}" >> ${OUTDIR}/alignment_hq_stats.txt
done;

# Get deduped alignment number
rm -f ${OUTDIR}/alignment_dedup_stats.txt
touch ${OUTDIR}/alignment_dedup_stats.txt
for inf in $(ls ${INDIR}/*/accepted_hits.hqmapped_dedup.bam); do
 c=$(samtools view -c $inf)
 echo -e "${inf}\t${c}" >> ${OUTDIR}/alignment_dedup_stats.txt
done;

# Get multimap stats
rm -f ${OUTDIR}/alignment_multimap_stats.txt
touch ${OUTDIR}/alignment_multimap_stats.txt
for inf in $(ls ${INDIR}/*/accepted_hits.bam); do
  bn=$(basename $(dirname $inf))
  echo -e "Doing file $inf"
  samtools view -F 0x100 $inf | awk '{print $5}' | sort | uniq -c  | sed "s,^,${bn} ," >> ${OUTDIR}/alignment_multimap_stats.txt
done
