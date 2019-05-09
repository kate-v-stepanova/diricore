#!/bin/bash

set -e
set -u

dataset_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
INDIR="$PROJECT_DIR/analysis/input/fastq"
STAR_GENOME_DIR="/icgc/dkfzlsdf/analysis/OE0532/star/star_index"

if [[ "$#" -eq 2 ]]; then
    OUTDIR=$2
else
    OUTDIR=$PROJECT_DIR/analysis/output/alignments
fi

echo "Unzip the fastq files in $INDIR to proceed"

module load STAR
mkdir -p $OUTDIR

for f in $(ls ${INDIR}); do
    b=$(basename ${f});
    b=${b%%.*};
    STAR --genomeDir $STAR_GENOME_DIR \
	--runThreadN 100 \
	--readFilesIn ${INDIR}/${f} \
	--outFileNamePrefix $OUTDIR/star_${b}_ \
	--outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM
done;


for f in $(ls ${OUTDIR}/*Aligned.toTranscriptome.out.bam); do
    b=$(basename ${f});
    b=${b%%_umi_extract*}; # will get a sample name: GFP_IP_RNAse_3X
    b=${b#*star_dem_};
    b=${b%%_Aligned.toTranscriptome*}; # will get star_dem_GFP_IP_RNAse_3X_umi_extracted
    mv $f ${OUTDIR}/${b}_toTranscriptome.bam 
done;

echo "Done with alignments; now filtering out non-primary/low-quality alignments";
# isolate "hqmapped" reads
ls -1 ${OUTDIR}/*toTranscriptome.bam | while read fn; do
    echo "Processing $fn "
    b=$(basename ${fn})
    b=${b%%".bam"}
    of="${OUTDIR}/${b}.hqmapped.bam"
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
