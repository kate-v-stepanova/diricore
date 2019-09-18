#!/bin/bash

set -e;
set -u;

dataset_id=$1
genome_dir="toGenome"
if [[ $# -ge 3 ]]; then
   if [[ $3 == "--transcriptome" ]]; then
       genome_dir="toTranscriptome"
   elif [[ $3 == "--toTranscriptome" ]]; then
       genome_dir="toTranscriptome"
   elif [[ $3 == "--totranscriptome" ]]; then
       genome_dir="toTranscriptome"
   elif [[ $3 == "--trans" ]]; then
       genome_dir="toTranscriptome"
   elif [[ $3 == "trans" ]]; then
       genome_dir="toTranscriptome"
   fi
fi
echo $genome_dir    

bam_type="hq"
# can be: hq, all
if [[ $# -ge 2 ]]; then
   if [[ $2 == "all" ]]; then
       bam_type="all"
   fi
fi 
echo $bam_type

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
INDIR="$PROJECT_DIR/analysis/output/alignments/$genome_dir"
OUTDIR="$PROJECT_DIR/analysis/output/alignments/$genome_dir"

mkdir -p $OUTDIR/logs


if [[ $bam_type == "hq" ]]; then
  ## HQ mapped reads
  echo "HQ mapped reads"
  # index
  for f in $(ls ${INDIR}/*.hqmapped.bam); do
    samplename=$(basename $f);
    samplename=${samplename%".hqmapped.bam"}
    outfile="${OUTDIR}/${samplename}.hqmapped.bam.bai";
    echo "Indexing ${f}";
    samtools index ${f}
  done;

  for f in $(ls ${INDIR}/*.hqmapped.bam); do
    samplename=$(basename $f)
    samplename=${samplename%".bam"}
    echo "Deduplicating ${f}";
    umi_tools dedup --log2stderr -L ${OUTDIR}/logs/${samplename}_umidedup_log.txt -I ${f} -S ${OUTDIR}/${samplename}_dedup.bam --output-stats=${OUTDIR}/logs/${samplename}_dedup.log
  done;
fi

if [[ $bam_type == "all" ]]; then
  ## All reads
  echo "All reads"
  # index
  for f in $(ls ${INDIR}/*${genome_dir}.bam); do
    samplename=$(basename $f);
    samplename=${samplename%"${genome_dir}.bam"}
    outfile="${OUTDIR}/${samplename}_${genome_dir}.bam.bai";
    echo "Indexing ${f}";
    samtools index ${f}
  done;

  for f in $(ls ${INDIR}/*${genome_dir}.bam); do
    samplename=$(basename $f)
    samplename=${samplename%".bam"}
    echo "Deduplicating ${f}";
    umi_tools dedup --log2stderr -L ${OUTDIR}/logs/${samplename}_umidedup_log.txt -I ${f} -S ${OUTDIR}/${samplename}_dedup.bam --output-stats=${OUTDIR}/logs/${samplename}_dedup.log
  done;
fi
