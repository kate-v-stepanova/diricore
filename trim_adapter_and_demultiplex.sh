#!/bin/bash

# exit on error:
set -e


flowcell_path=$1
flowcell_id=$(basename ${flowcell_path})
project_path=$(dirname ${flowcell_path})

diricore_path="/home/e984a/diricore"
merged_path="$project_path/analysis/input/merged"
merged_file="${merged_path}/${flowcell_id}.fastq.gz"
demultiplexed_path="$project_path/analysis/output/demultiplexed"
umi_extract_path="$project_path/analysis/output/umi_extract"
umi_extract_logs="${umi_extract_path}/logs"
fastq_path="$project_path/analysis/input/fastq"

trimmed_file="$merged_path/${flowcell_id}_trimmed.fastq"

bc_file="$project_path/analysis/input/metadata/bc_file.txt"
cutadapt_trimming_stats="$project_path/analysis/output/cutadapt_trimming_stats.txt"
bc_split_stats="$project_path/analysis/output/bc_split_stats.txt"

adapter_sequence="AGATCGGAAGAGCACACGTCTGAAX"
bc_pattern="NNNNNNNNNN" # random sequence 5nt long
# bc_pattern="NNNNN"
barcode_splitter="$diricore_path/programs/fastx_toolkit/fastx_barcode_splitter.pl"

# Creating directory structure (if not exist)
mkdir -p $merged_path
mkdir -p $demultiplexed_path
mkdir -p $umi_extract_logs
mkdir -p $fastq_path

# merge fastq files into one
echo "Merging $flowcell_path/*/*/*.fastq.gz into one file"
cat $flowcell_path/*/*/*.fastq.gz > $merged_file
echo "Merging done. Created file: ${merged_file}"

# remove adapter
echo "Removing adapter"
echo "gzip -dc $merged_file | cutadapt -u 3 -O 7 -m 30 -a $adapter_sequence --discard-untrimmed -o $trimmed_file - 2> $cutadapt_trimming_stats"

echo "cat $trimmed_file | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_ --suffix .fastq --eol > $bc_split_stats "
# both commands together
# gzip -dc $merged_file | cutadapt -u 3 -O 7 -m 30 -a $adapter_sequence --discard-untrimmed - 2> $cutadapt_trimming_stats | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_ --suffix .fastq --eol > $bc_split_stats

echo "Removing adapter done. Outfiles: $demultiplexed_path/*.fastq"  
# exit
# remove umi (any random sequence of 5-nt length)
echo "Removing UMIs"
for f in $demultiplexed_path/*.fastq; do
    fb=$(basename ${f});
  b=${fb%%.*};
  if [[ "$b" != "dem_unmatched" ]]; then
      echo "Processing ${f}";
      umi_tools extract --extract-method=string --3prime --bc-pattern=$bc_pattern --log=$umi_extract_logs/${b}.log --log2stderr --stdin=$demultiplexed_path/${fb} --stdout=$umi_extract_path/${b}_umi_extracted.fastq.gz
      echo "Done ${f}"

  fi
  done;

# the result of umi_extract will become input for the next step (preprocessing)
echo "Moving files to $fastq_path"
# mv or cp?
mv $umi_extract_path/*.fastq.gz $fastq_path
echo "Done"
