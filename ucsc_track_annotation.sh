#!/bin/bash

project_id=$1
species=$2

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$project_id"
INDIR="$PROJECT_DIR/analysis/output/gen_tracks"
OUTFILE="$INDIR/ucsc_track_annotation.txt"

echo "browser" > $OUTFILE
for f in $(ls $INDIR/*.bw); do
    sample=$(basename $f)
    sample=${sample%".bw"}
    echo "track name=$sample db=$species type=bigWig bigDataUrl=ftp://ftp.dkfz-heidelberg.de/outgoing/B250/ucsc/$project_id/${sample}.bw" >> $OUTFILE
done
echo "Created file: $OUTFILE"
