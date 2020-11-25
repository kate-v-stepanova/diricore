#!/bin/bash

set -e
set -u

dataset_id=$1
if [ "$#" -ge 2 ]; then
   genome=$2
else
   genome="hg19"
fi

module load bedtools

bam_type="hq_unique"
bam_pattern="_toGenome.hqmapped_dedup.bam"
# can be: hq, hq_unique, all, all_unique
if [ $# -ge 3 ]; then
   bam_type=$3
fi

if [[ $bam_type == "hq" ]]; then
   bam_pattern="_toGenome.hqmapped.bam"
elif [[ $bam_type == "all" ]]; then
   bam_pattern="_toGenome.bam"
elif [[ $bam_type == "all_unique" ]]; then
   bam_pattern="_toGenome_dedup.bam"
fi

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
DIRICORE_DIR="/icgc/dkfzlsdf/analysis/OE0532/software/diricore"

chrom_sizes="$BASE_DIR/static/$genome/$genome.chrom.sizes"
#INDIR="$PROJECT_DIR/analysis/output/tophat_out"
INDIR="$PROJECT_DIR/analysis/output/alignments/toGenome"
OUTDIR="$PROJECT_DIR/analysis/output/gen_tracks/$bam_type"

BDG2BW="$DIRICORE_DIR/utils/bdg2bw.sh"


mkdir -p $OUTDIR
# Getting bdg files
for f in $(ls $INDIR/*"$bam_pattern"); do
  echo "Calculating normalization factor for $f"
  samplename=$(basename $f)
  samplename=${samplename%"$bam_pattern"}
  echo $samplename
  fact=$(samtools view -c $f) &&
  fact=$(echo "scale=6; 1000000.0 / $fact" | bc)
  echo "Done: ${fact}. Writing bdg files: $OUTDIR/${samplename}_minus.bdg $OUTDIR/${samplename}_plus.bdg"
  samtools view -b $f | tee >(bedtools genomecov -ibam /dev/stdin -g $chrom_sizes -scale $fact -bg -split -strand + | sort -k1,1 -k2,2n > $OUTDIR/${samplename}_plus.bdg) | bedtools genomecov -ibam /dev/stdin -g $chrom_sizes -scale $fact -bg -split -strand - | sort -k1,1 -k2,2n > $OUTDIR/${samplename}_minus.bdg
done;

# Convert to bigwig
for f in $(ls ${OUTDIR}/*.bdg); do
  bn=$(basename ${f})
  echo "Converting to bigwig "
  $BDG2BW $f $chrom_sizes
  echo "Done: $(ls ${OUTDIR}/${bn%.*}.bw)"
done
