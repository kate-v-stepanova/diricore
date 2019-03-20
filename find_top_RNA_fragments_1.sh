#!/bin/bash

if [[ $# -eq 1 ]]; then
  counts=$1
else
  # I guess this is correct. 
  echo "WARNING: use $0 100000 to filter out the reads with less than 100000. For miseq run this number should be much less. If you don't enter this number, the default value (100 000) will be used. But then you might get an empty run_UMIs_top_rRNA_seqs.txt file"
  counts=100000
fi

RRNA_REF="/home/e984a/diricore/staticdata/human/rRNAs"
OUTDIR="./data/output/rrna_fragments"
BOWTIE="/home/e984a/diricore/programs/bowtie2-2.0.6/bowtie2"
mkdir -p $OUTDIR

# zcat data/input/fastq/*.fastq.gz | ./programs/bowtie2-2.0.6/bowtie2 --seed 42 -p 1 --local "ParseError: KaTeX parse error: Expected 'EOF', got '\ ' at position 15: {RRNA_REF}" - | samtools vi…10}' | sort | uniq -c | awk 'ParseError: KaTeX parse error: Expected '}', got 'EOF' at end of input: …> 100000 {print1,$2}' | sort -k1nr ${OUTDIR}/run_UMIs_top_rRNA_seqs.txt

zcat data/input/fastq/*.fastq.gz | /home/e984a/diricore/programs/bowtie2-2.0.6/bowtie2 --seed 42 -p 1 --local  /home/e984a/diricore/staticdata/human/rRNAs - > $OUTDIR/bowtie_out.tmp.bam

samtools view -F 0x4 $OUTDIR/bowtie_out.tmp.bam | awk '{print $10}' | sort | uniq -c | awk '$1 >= '$counts' {print $1,$2}' | sort -k1nr > $OUTDIR/run_UMIs_top_rRNA_seqs.txt

rm -f $OUTDIR/bowtie_out.tmp.bam

echo "Created file: $OUTDIR/run_UMIs_top_rRNA_seqs.txt"
# samtools view -F 0x4 $OUTDIR/bowtie_out.tmp.bam | awk '{print $10}' | sort | uniq -c | awk '$1 > 100000 {print $1,$2}' | sort -k1nr > $OUTDIR/run_UMIs_top_rRNA_seqs.txt

# echo "zcat data/input/fastq/*.fastq.gz | ${BOWTIE} --seed 42 -p 1 --local  "${RRNA_REF}" - | samtools view -F 0x4 | awk '{print $10}' | sort | uniq -c | awk '$1 > 100000 {print $1,$2}' | sort -k1nr  > ${OUTDIR}/run_UMIs_top_rRNA_seqs.txt"
