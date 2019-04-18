#!/bin/bash


dataset_id=$1

if [[ $# -ge 2 ]]; then
  counts=$2
else
  # I guess this is correct. 
  echo "WARNING: use $0 100000 to filter out the reads with less than 100000. For miseq run this number should be much less. If you don't enter this number, the default value (100 000) will be used. But then you might get an empty run_UMIs_top_rRNA_seqs.txt file"
  counts=100
fi

prefix="rRNA"
trna=0
if [[ $# -ge 3 ]]; then
   trna=$3
   prefix="tRNA"
fi

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
DIRICORE_DIR="/home/e984a/diricore"
REF="$DIRICORE_DIR/staticdata/human/${prefix}s"
BOWTIE="$DIRICORE_DIR/programs/bowtie2-2.0.6/bowtie2"
INDIR="$PROJECT_DIR/analysis/input/fastq"

if [[ $trna == 0 ]]; then
    OUTDIR="$PROJECT_DIR/analysis/output/rrna_fragments"
else 
    OUTDIR="$PROJECT_DIR/analysis/output/trna_fragments"
fi

mkdir -p $OUTDIR

for f in $(ls ${INDIR}/*.fastq.gz); do
    echo "Processing $f"
    sample_name="$(basename $f)"
    sample_name="${sample_name%%.*}"
    tmp_file="$OUTDIR/${sample_name}_bowtie_out.tmp.bam"
    out_file="$OUTDIR/${sample_name}_top_${prefix}_seqs.txt"
    zcat $f | $BOWTIE --seed 42 -p 1 --local  $REF - > $tmp_file
    # echo "Counts	Seq" > out_file
    samtools view -F 0x4 $tmp_file | awk '{print $10}' | sort -T $OUTDIR | uniq -c | awk '$1 >= '$counts' {print $1,$2,$3,$4}' | sort -k1nr > $out_file
    rm -f $tmp_file
    echo "Created $out_file"
done
