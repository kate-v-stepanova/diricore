#!/bin/bash

set -e
set -u

dataset_id=$1
genome=$2

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
INDIR="$PROJECT_DIR/analysis/input/fastq"
STAR_GENOME_DIR="/icgc/dkfzlsdf/analysis/OE0532/static/$genome"
OUTDIR=$PROJECT_DIR/analysis/output/alignments

echo "Unzip the fastq files in $INDIR to proceed"

module load STAR
mkdir -p $OUTDIR

for f in $(ls ${INDIR}/*.fastq.gz); do
    b=$(basename ${f});
    b=${b%%.*};
    STAR --genomeDir $STAR_GENOME_DIR \
	--runThreadN 100 \
	--readFilesIn ${f} \
	--outFileNamePrefix $OUTDIR/star_${b}_ \
	--outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM \
        --readFilesCommand zcat
done;

echo "Sorting bam files"
module load samtools/1.6
for f in $(ls ${OUTDIR}/*.bam); do
    b=$(basename ${f});
    b=${b%%.bam}
    sorted=$OUTDIR/$b.sorted.bam
    samtools sort $f -o $sorted
    mv $sorted $f
done;

mkdir -p $OUTDIR/toTranscriptome
mkdir -p $OUTDIR/toGenome
for f in $(ls ${OUTDIR}/*Aligned.toTranscriptome.out.bam); do
    b=$(basename ${f});
    b=${b%%_umi_extract*}; # will get a sample name: GFP_IP_RNAse_3X
    b=${b#*star_dem_};
    b=${b%%_Aligned.toTranscriptome*}; # will get star_dem_GFP_IP_RNAse_3X_umi_extracted
    echo "mv $f ${OUTDIR}/toTranscriptome/${b}_toTranscriptome.bam"
    mv $f ${OUTDIR}/toTranscriptome/${b}_toTranscriptome.bam 
done;

for f in $(ls ${OUTDIR}/*Aligned.out.bam); do
    b=$(basename $f);
    b=${b%%_umi_extract*}
    b=${b#*star_dem_}
    b=${b%%_Aligned.out.bam}
    echo "mv $f ${OUTDIR}/toGenome/${b}_toGenome.bam"
    mv $f ${OUTDIR}/toGenome/${b}_toGenome.bam
done

echo "Done with alignments; now filtering out non-primary/low-quality alignments (for Transciptome files)";
# isolate "hqmapped" reads
ls -1 ${OUTDIR}/toTranscriptome/*toTranscriptome.bam | while read fn; do
    echo "Processing $fn "
    b=$(basename ${fn})
    b=${b%%".bam"}
    of="${OUTDIR}/toTranscriptome/${b}.hqmapped.bam"
#    dn=$(dirname "$fn");
#    of="${dn}/accepted_hits.hqmapped.bam";
    cat \
        <(samtools view -H "${fn}") \
        <(cat "${fn}" | samtools view -q10 -F260 -) \
    | samtools view -bS - \
    > "${of}";

  #  samtools sort $fn > "${OUTDIR}/${b}.sorted.bam"
  #  samtools sort $of > "${OUTDIR}/${b}.hqmapped.sorted.bam"

done
echo "Done filtering";

echo "Done with alignments; now filtering out non-primary/low-quality alignments (for Genome files)";
# isolate "hqmapped" reads
ls -1 ${OUTDIR}/toGenome/*toGenome.bam | while read fn; do
    echo "Processing $fn "
    b=$(basename ${fn})
    b=${b%%".bam"}
    of="${OUTDIR}/toGenome/${b}.hqmapped.bam"
#    dn=$(dirname "$fn");
#    of="${dn}/accepted_hits.hqmapped.bam";
    cat \
        <(samtools view -H "${fn}") \
        <(cat "${fn}" | samtools view -q10 -F260 -) \
    | samtools view -bS - \
    > "${of}";

  #  samtools sort $fn > "${OUTDIR}/${b}.sorted.bam"
  #  samtools sort $of > "${OUTDIR}/${b}.hqmapped.sorted.bam"

done
echo "Done filtering";
