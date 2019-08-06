#!/bin/bash

dataset_id=$1

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
INDIR="$PROJECT_DIR/analysis/output/alignments/toGenome"
#INDIR="$PROJECT_DIR/analysis/output/alignments/toTranscriptome"
OUTDIR="$PROJECT_DIR/analysis/output/alignment_stats"
clean_dir="$PROJECT_DIR/analysis/output/clean"

mkdir -p ${OUTDIR}

# Get HQ (with duplicates)
rm -f ${OUTDIR}/hq_stats.txt
for inf in $(ls ${INDIR}/*.hqmapped.bam); do
    fn=$(basename $inf)
    fn=${fn%"_toGenome.hqmapped.bam"}
    c=$(samtools view -c $inf)
    echo -e "${fn}\t${c}" >> ${OUTDIR}/hq_stats.txt
done;
echo "Created file: ${OUTDIR}/hq_stats.txt"

# Get HQ unique
rm -f ${OUTDIR}/hq_unique_stats.txt
for inf in $(ls ${INDIR}/*.hqmapped_dedup.bam); do
      fn=$(basename $inf)
      fn=${fn%"_toGenome.hqmapped_dedup.bam"}
      c=$(samtools view -c $inf)
      echo -e "${fn}\t${c}" >> ${OUTDIR}/hq_unique_stats.txt
done;
echo "Created file: ${OUTDIR}/hq_unique_stats.txt"

# Get ALL stats
rm -f ${OUTDIR}/all_stats.txt
for inf in $(ls $INDIR/*_toGenome.bam); do
    echo -e "Doing file $inf"
    bn=$(basename $inf)
    bn=${bn%"_toGenome.bam"}
    c=$(samtools view -c $inf)
    echo -e "${bn}\t$c" >> ${OUTDIR}/all_stats.txt
done
echo "Created file: ${OUTDIR}/all_stats.txt"

# Get ALL unique stats
rm -f ${OUTDIR}/all_unique_stats.txt
for inf in $(ls $INDIR/*_toGenome_dedup.bam); do
    bn=$(basename $inf)
    bn=${bn%"_toGenome_dedup.bam"}
    c=$(samtools view -c $inf)
    echo -e "${bn}\t$c" >> ${OUTDIR}/all_unique_stats.txt
done
echo "Created file: ${OUTDIR}/all_unique_stats.txt"


# Get rRNA stats
rm -f ${OUTDIR}/rrna_stats.txt
for f in $(ls $clean_dir/*.rrna.err); do
#for f in $(ls $clean_dir/*.fastq.gz); do
    fn=$(basename $f)
    fn=${fn%".rrna.err"}
#    fn=${fn%".fastq.gz"}
    one_time=$(grep "aligned exactly 1 time" < $f | cut -d " " -f5)
    more_times=$(grep "aligned >1 times" < $f | cut -d " " -f5)
    rrna=$(echo $one_time+$more_times | bc)
#    reads=$(zcat $f | wc -l)
#    by4=4
#    rrna=$((reads / by4))
    echo -e "$fn\t$rrna" >> ${OUTDIR}/rrna_stats.txt
done
echo "Created file: ${OUTDIR}/rrna_stats.txt"

# Get tRNA stats
rm -f ${OUTDIR}/trna_stats.txt
for f in $(ls $clean_dir/*.trna.err); do
    fn=$(basename $f)
    fn=${fn%".trna.err"}
    one_time=$(grep "aligned exactly 1 time" < $f | cut -d " " -f5)
    more_times=$(grep "aligned >1 times" < $f | cut -d " " -f5)
    rrna=$(echo $one_time+$more_times | bc)
    echo -e "$fn\t$rrna" >> ${OUTDIR}/trna_stats.txt
done
echo "Created file: ${OUTDIR}/trna_stats.txt"
