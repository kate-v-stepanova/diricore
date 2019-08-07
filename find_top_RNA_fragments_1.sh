#!/bin/bash


dataset_id=$1
genome="hg19"
if [[ $# -ge 2 ]]; then
  genome=$2
fi


if [[ $# -ge 3 ]]; then
  counts=$3
else
  # I guess this is correct. 
  echo "WARNING: use $0 100000 to filter out the reads with less than 100000. For miseq run this number should be much less. If you don't enter this number, the default value (100 000) will be used. But then you might get an empty run_UMIs_top_rRNA_seqs.txt file"
  counts=1
fi

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
DIRICORE_DIR="/home/e984a/diricore"

BOWTIE="$DIRICORE_DIR/programs/bowtie2-2.0.6/bowtie2"
INDIR="$PROJECT_DIR/analysis/input/fastq"


prefix="rRNA"
trna=0
OUTDIR="$PROJECT_DIR/analysis/output/rrna_fragments"
if [[ $# -ge 3 ]]; then
   trna=$3
   if [[ $trna == 'trna' ]]; then
       prefix="tRNA"
       OUTDIR="$PROJECT_DIR/analysis/output/trna_fragments"
   else
       prefix="rRNA"
       OUTDIR="$PROJECT_DIR/analysis/output/rrna_fragments"
   fi
fi

REF="$BASE_DIR/static/$genome/${prefix}s"

mkdir -p $OUTDIR

for f in $(ls ${INDIR}/*.fastq.gz); do
    echo "Processing $f"
    sample_name="$(basename $f)"
    sample_name="${sample_name%%.*}"
    tmp_file="$OUTDIR/${sample_name}_bowtie_out.tmp.bam"
    out_file="$OUTDIR/${sample_name}_top_${prefix}_seqs.txt"
    rm $out_file
    touch $out_file
    zcat $f | $BOWTIE --seed 42 -p 1 --local  $REF - > $tmp_file
    # echo "Counts	Seq" > out_file
    samtools view -F 0x4 $tmp_file | awk '{print $10}' | sort -T $OUTDIR | uniq -c | awk '$1 >= '$counts' {print $1,$2,$3,$4}' | sort -k1nr >> $out_file
    rm -f $tmp_file
    echo "Created $out_file"
done
