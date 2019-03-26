#!/bin/bash


dataset_id=$1

if [[ $# -ne 2 ]]; then
  # I guess this is correct. 
  echo "WARNING: use $0 100000 to filter out the reads with less than 100000. For miseq run this number should be much less. If you don't enter this number, the default value (100 000) will be used. But then you might get an empty run_UMIs_top_rRNA_seqs.txt file"
  counts=100000
else
 counts=$2
fi

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
DIRICORE_DIR="/home/e984a/diricore"
RRNA_REF="$DIRICORE_DIR/staticdata/human/rRNAs"
BOWTIE="$DIRICORE_DIR/programs/bowtie2-2.0.6/bowtie2"
OUTDIR="$PROJECT_DIR/analysis/output/rrna_fragments"
INDIR="$PROJECT_DIR/analysis/input/fastq"

mkdir -p $OUTDIR

zcat $INDIR/*.fastq.gz | $BOWTIE --seed 42 -p 1 --local  $RRNA_REF - > $OUTDIR/bowtie_out.tmp.bam

samtools view -F 0x4 $OUTDIR/bowtie_out.tmp.bam | awk '{print $10}' | sort | uniq -c | awk '$1 >= '$counts' {print $1,$2}' | sort -k1nr > $OUTDIR/run_UMIs_top_rRNA_seqs.txt

rm -f $OUTDIR/bowtie_out.tmp.bam

echo "Created file: $OUTDIR/run_UMIs_top_rRNA_seqs.txt"
