#!/bin/bash

# exit on error:
set -e

project_id=$1
flowcell_id=$2

BASE_PATH="/icgc/dkfzlsdf/analysis/OE0532"
project_path="$BASE_PATH/$project_id"

merged_path="$project_path/analysis/input/merged"
merged_file="${merged_path}/${flowcell_id}.fastq.gz"

# Creating directory structure (if not exist)
mkdir -p $merged_path

echo "Merging $flowcell_path/*/*/*.fastq.gz into one file"
cat $flowcell_path/*/*/*.fastq.gz > $merged_file
echo "Merging done. Created file: ${merged_file}"
