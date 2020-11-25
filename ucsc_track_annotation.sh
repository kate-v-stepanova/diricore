#!/bin/bash

project_id=$1
species=$2
bam_type="hq_unique"
if [[ $# -ge 3 ]]; then
    bam_type=$3
fi

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$project_id"
INDIR="$PROJECT_DIR/analysis/output/gen_tracks/$bam_type"
OUTFILE="$INDIR/ucsc_track_annotation.txt"

echo "browser" > $OUTFILE
for f in $(ls $INDIR/*.bw); do
    sample=$(basename $f)
    sample=${sample%".bw"}
    echo "track name=$sample db=$species type=bigWig bigDataUrl=ftp://ftp.dkfz-heidelberg.de/outgoing/B250/ucsc/$project_id/$bam_type/${sample}.bw visibility=full alwaysZero=on yLineOnOff=on yLineMark=0 maxHeightPixels=40" >> $OUTFILE
done
echo "Created file: $OUTFILE"
