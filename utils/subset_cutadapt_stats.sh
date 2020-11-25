#!/bin/bash

set -e
set -u

dataset_id=$1
subset=$2

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
project_dir="$BASE_DIR/$dataset_id"
#cutadapt_file="$project_dir/analysis/output/cutadapt_trimming_stats.txt"
stats_file="$project_dir/analysis/output/${subset}_cutadapt_plot_stats.txt"
tmp_file1="$project_dir/analysis/output/${subset}_tmp_file1.txt"
tmp_file2="$project_dir/analysis/output/${subset}_tmp_file2.txt"

rm -f $tmp_file1
rm -f $tmp_file2
for f in $(ls $project_dir/analysis/output/${subset}_cutadapt_trimming_stats*.txt); do
     if [[ ! -f $tmp_file1 ]]; then
         awk '/Total reads processed/,/Reads written/' $f | cut -d ":" -f2 | cut -d "(" -f1 | tr -d , | tr -d " " | tr -d "\t"> $tmp_file1
     else
         awk '/Total reads processed/,/Reads written/' $f | cut -d ":" -f2 | cut -d "(" -f1 | tr -d , | tr -d " " | tr -d "\t" > $tmp_file2
         paste -d+ $tmp_file1 $tmp_file2 | bc > $stats_file
         cp $stats_file $tmp_file1
     fi
done
echo -e "Total reads:\nWith adapters:\nToo short:\nPassed filters:" > ${tmp_file2}
paste ${tmp_file2} $tmp_file1 > $stats_file

rm -f $tmp_file1
rm -f $tmp_file2

echo "Created file:  $stats_file"
