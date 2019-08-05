#!/bin/bash

project_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
indir="$BASE_DIR/$project_id/analysis/output/demultiplexed"
outdir="$BASE_DIR/$project_id/analysis/output/trimmed"
adapter="TGGAATTCTCGGGTGCCA"
echo "Removing adapter"
for f in $(ls $indir/*.fastq.gz); do
    echo $f
    samplename=$(basename $f);
    samplename=${samplename%".fastq.gz"};
    echo $samplename
    gzip -dc $f | cutadapt -O 7 -m 20 -a $adapter --discard-untrimmed -o $outdir/${samplename}.fastq - > $BASE_DIR/$project_id/analysis/output/${samplename}_cutadapt_trimming_stats.txt
    echo "Removing adapter done. Outfiles: $outdir/*.fastq.";
done;
