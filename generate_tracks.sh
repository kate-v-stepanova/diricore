#!/bin/bash

set -e
set -u

dataset_id=$1
if [ "$#" -ge 2 ]; then
   genome=$2
else
   genome="hg19"
fi

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
DIRICORE_DIR="/home/e984a/diricore"

chrom_sizes="$BASE_DIR/static/$genome/$genome.chrom.sizes"
INDIR="$PROJECT_DIR/analysis/output/tophat_out"
OUTDIR="$PROJECT_DIR/analysis/output/gen_tracks"

BDG2BW="$DIRICORE_DIR/utils/bdg2bw.sh"

mkdir -p $OUTDIR
# Getting bdg files
for f in $(ls ${INDIR}); do
  infile="${INDIR}/$f/accepted_hits.hqmapped_dedup.bam"
  echo "Calculating normalization factor for $f"
  samplename=${f#dem_}
  samplename=${samplename%_umi_extracted}
  echo $samplename
  fact=$(samtools view -c ${infile}) &&
  fact=$(echo "scale=6; 1000000.0 / $fact" | bc)
  echo "Done: ${fact}. Writing bdg files: $OUTDIR/${samplename}_minus.bdg $OUTDIR/${samplename}_plus.bdg"
  samtools view -b $infile | tee >(bedtools genomecov -ibam /dev/stdin -g $chrom_sizes -scale $fact -bg -split -strand + | sort -k1,1 -k2,2n > $OUTDIR/${samplename}_plus.bdg) | bedtools genomecov -ibam /dev/stdin -g $chrom_sizes -scale $fact -bg -split -strand - | sort -k1,1 -k2,2n > $OUTDIR/${samplename}_minus.bdg
done;

# Convert to bigwig
for f in $(ls ${OUTDIR}/*.bdg); do
  bn=$(basename ${f})
  echo "Converting to bigwig "
  $BDG2BW $f $chrom_sizes
  echo "Done: $(ls ${OUTDIR}/${bn%.*}.bw)"
done
