#!/bin/bash

set -e
set -u

dataset_id=$1

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"

OUTDIR="$PROJECT_DIR/analysis/output/gen_tracks"
bc_file="$PROJECT_DIR/analysis/input/metadata/bc_file.txt"

rm -f ${OUTDIR}/UCSC_TrackDb.txt
touch ${OUTDIR}/UCSC_TrackDb.txt

echo -e "
  track run_UMIs
  superTrack on show
  shortLabel run_UMIs
  longLabel run_UMIs
  html ../track_help.html
  autoScale On
  yLineOnOff on
  alwaysZero on
  " > $OUTDIR/UCSC_TrackDb.txt

samples=$(awk '{print $1}' ${bc_file})

for i in $samples; do
  echo -e "
  \ttrack $i
  \tparent run_UMIs
  \tcontainer multiWig
  \tshortLabel $i
  \tlongLabel $i
  \ttype bigWig
  \taggregate solidOverlay
  \tshowSubtrackColorOnUi on
  \tvisibility full
  \tmaxHeightPixels 100:40:10

  \t\ttrack ${i}_plus.bw
  \t\tbigDataUrl run_UMIs/${i}_plus.bw
  \t\tparent $i
  \t\ttype bigWig
  \t\tshortLabel ${i}_+
  \t\tLongLabel ${i}_+
  \t\tcolor 19,80,181

  \t\ttrack ${i}_minus.bw
  \t\tbigDataUrl run_UMIs/${i}_minus.bw
  \t\tparent $i
  \t\ttype bigWig
  \t\tshortLabel ${i}_-
  \t\tLongLabel ${i}_-
  \t\tnegateValues on
  \t\taltColor 206,58,0

  \t\t####################################
  "
done >> ${OUTDIR}/UCSC_TrackDb.txt
