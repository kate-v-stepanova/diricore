#!/bin/bash

set -e
set -u

dataset_id=$1

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
project_dir="$BASE_DIR/$dataset_id"
cutadapt_file="$project_dir/analysis/output/cutadapt_trimming_stats.txt"
stats_file="$project_dir/analysis/output/cutadapt_plot_stats.txt"

echo "Total reads:
With adapters:
Too short:
Passed filters:" > /tmp/header.txt
awk '/Total reads processed/,/Reads written/' $cutadapt_file | cut -d ":" -f2 | cut -d "(" -f1 | tr -d , | tr -d " " > /tmp/reads.txt # $stats_file
paste /tmp/header.txt /tmp/reads.txt > $stats_file
rm /tmp/header.txt /tmp/reads.txt
echo "Created file: $stats_file"
