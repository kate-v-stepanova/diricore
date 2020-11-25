#!/bin/bash

set -e
set -u

dataset_id=$1

outdir="/icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/output"
stats_file="$outdir/bc_split_stats.txt"

total=$(cat $outdir/*_bc_split_stats.txt | grep total | cut -f 2 | paste -s -d+ - | bc)
for f in $(ls $outdir/*_bc_split_stats.txt); do
     sample=$(basename $f);
     sample=${sample%"_bc_split_stats.txt"}
     unmatched_str=$(cat $f | grep unmatched)
     to_replace=$(echo ${unmatched_str//unmatched/"$sample"_unmatched})
     if [[ ! -f $stats_file ]]; then
         cat $f | grep -v total | grep -v unmatched > $stats_file
     else
         cat $f | grep -v total | grep -v unmatched | tail -n +2 >> $stats_file
     fi
     echo $to_replace | tr ' ' '\t' >> $stats_file
done

echo -e "total\t$total" >> $stats_file

echo "Created file:  $stats_file"
