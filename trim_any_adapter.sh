#!/bin/bash

project_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
indir="$BASE_DIR/$project_id/analysis/output/demultiplexed"
outdir="$BASE_DIR/$project_id/analysis/output/trimmed"


new_adapter="AGATCGGAAGAGCACACGTCTGAA"
old_adapter="TGGAATTCTCGGGTGCCA"
very_old_adapter="TCGTATGCCGTCTTCTGCTTGA"
yeast="CTGTAGGCACCATCAAT"
adapter=$new_adapter

# add to cutadapt parameter "-u 3" to remove 3 nts from 3' end - this is only for a NEW protocol
u=" -u 3 " 
if [[ $# -ge 2 ]]; then
    u="" # do not use parameter "-u 3" for any other protocol
    if [[ $2 == 'old' ]]; then
        adapter=$old_adapter
    elif [[ $2 == "old_adapter" ]]; then
        adapter=$old_adapter
    elif [[ $2 == "very_old" ]]; then
        adapter=$very_old_adapter
    elif [[ $2 == "very_old_adapter" ]]; then
        adapter=$very_old_adapter
    elif [[ $2 == "yeast" ]]; then
        adapter=$yeast
    fi
fi

read_length=30
if [[ $# -ge 3 ]]; then
    read_length=$3
fi

symlinks_dir="$BASE_DIR/$project_id/analysis/input/fastq"
script_dir="$BASE_DIR/tmp/$project_id/faster_trimming"

mkdir -p $script_dir
mkdir -p $outdir
mkdir -p $symlinks_dir

for f in $(ls $indir/*.fastq.gz); do
    samplename=$(basename $f);
    samplename=${samplename%".fastq.gz"};
    script_file="$script_dir/${samplename}.sh"
    echo "gzip -dc $f | cutadapt -O 7 -m $read_length $u -a $adapter --discard-untrimmed -o $outdir/${samplename}.fastq - > $BASE_DIR/$project_id/analysis/output/${samplename}_cutadapt_trimming_stats.txt" > $script_file
    echo "gzip $outdir/${samplename}.fastq" >> $script_file
    echo "ln -s $outdir/${samplename}.fastq.gz $symlinks_dir/${samplename}.fastq.gz" >> $script_file
    chmod +x $script_file
    echo "bsub -q medium $script_file"
done;

