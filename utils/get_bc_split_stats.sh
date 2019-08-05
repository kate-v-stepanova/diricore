#!/bin/bash

project_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
project_dir="$BASE_DIR/$project_id"
#dem_path="$project_dir/analysis/input/fastq"
#dem_path="$project_dir/analysis/output/demultiplexed"
bc_split_file="$project_dir/analysis/output/bc_split_stats.txt"

dem_path="$project_dir/analysis/output/demultiplexed"
# $2 can be: demultiplexed, trimmed
if [[ $# -ge 2 ]]; then
    if [[ $2 == "trimmed" ]]; then
          dem_path="$project_dir/analysis/output/trimmed"
    fi
fi

only_stdout=0
if [[ $# -ge 3 ]]; then
   if [[ "$3" == "only_stdout" ]]; then
       only_stdout=1
   fi
fi
if [[ $only_stdout == 0 ]]; then
    echo -e "Barcode\tCount\tLocation" > $bc_split_file
fi

for f in $(ls $dem_path/*fastq.gz); do
    lines=$(zcat $f | wc -l)
    divide_by=4
    reads=$((lines / divide_by))
    samplename=$(basename $f)
    samplename=${samplename#"dem_"}
    samplename=${samplename%".fastq.gz"}
    echo "$samplename	$reads	$f"
    if [[ $only_stdout == 0 ]]; then
        echo -e "$samplename\t$reads\t$f" >> $bc_split_file
        echo "Created file: $bc_split_file"
    fi
done
