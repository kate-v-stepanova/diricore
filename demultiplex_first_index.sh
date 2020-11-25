#!/bin/bash

# exit on error:
set -e

project_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
project_path="$BASE_DIR/$project_id"
# Check if merged file is present
merged_path="$project_path/analysis/input/merged"
echo "Checking if $merged_path is not empty:"
if [[ $(ls -A $merged_path) ]]; then
    echo "Fine"
else
    echo "Dir is empty. Exiting"
    exit
fi

#demultiplexed_path="$project_path/analysis/output/demultiplexed"
diricore_path="/icgc/dkfzlsdf/analysis/OE0532/software/diricore"
demultiplexed_path=$merged_path/dmx_1st_index


bc_path="$project_path/analysis/input/metadata/"

stats_dir="$demultiplexed_path"
bc_split_stats="$stats_dir/bc_split_stats.txt"

bc_pattern="NNNNNNNNNN" # random sequence 5nt long + barcode
barcode_splitter="$diricore_path/programs/fastx_toolkit/fastx_barcode_splitter.pl"

# Creating directory structure (if not exist)
mkdir -p $demultiplexed_path

# remove adapter
min_len=30
if [[ $# -ge 2 ]]; then
   min_len=$2
fi
for f in $(ls $merged_path/*.fastq.gz); do
    echo "Processing $f"
    fn=$(basename $f);
    fn=${fn%.fastq.gz}
    bc_split_stats="$stats_dir/${fn}_bc_split_stats.txt"
    bc_file="$bc_path/bc_file_${fn}.txt"
    echo "Demultiplexing: cat $f | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_ --suffix .fastq --eol --mismatches 1 > $bc_split_stats"
    zcat $f | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_ --suffix .fastq --eol --mismatches 1 > $bc_split_stats
    echo "Demultiplexing done"
    mv $demultiplexed_path/dem_unmatched.fastq $demultiplexed_path/dem_unmatched_${fn}.fastq
done
# cat $trimmed_file | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_ --suffix .fastq --eol > $bc_split_stats

