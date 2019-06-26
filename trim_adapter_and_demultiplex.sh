#!/bin/bash

# exit on error:
set -e

project_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
project_path="$BASE_DIR/$project_id"
# Check if merged file is present
merged_path="$project_path/analysis/input/merged"
merged_file="${merged_path}/*.fastq.gz" # should be just 1 file
echo "Checking if $merged_file is present:"
ls $merged_file 
echo "Fine"

demultiplexed_path="$project_path/analysis/output/demultiplexed"
umi_extract_path="$project_path/analysis/output/umi_extract"
umi_extract_logs="${umi_extract_path}/logs"
fastq_path="$project_path/analysis/input/fastq"
diricore_path="/home/e984a/diricore"

sample_id=$(basename $merged_file)
sample_id=${sample_id%".fastq.gz"}
trimmed_file="$merged_path/${sample_id}_trimmed.fastq"

bc_file="$project_path/analysis/input/metadata/bc_file.txt"
cutadapt_trimming_stats="$project_path/analysis/output/cutadapt_trimming_stats.txt"
bc_split_stats="$project_path/analysis/output/bc_split_stats.txt"

adapter_sequence="AGATCGGAAGAGCACACGTCTGAA"
bc_pattern="NNNNNNNNNN" # random sequence 5nt long + barcode
barcode_splitter="$diricore_path/programs/fastx_toolkit/fastx_barcode_splitter.pl"

# Creating directory structure (if not exist)
mkdir -p $demultiplexed_path
mkdir -p $umi_extract_logs
mkdir -p $fastq_path

# remove adapter
echo "Removing adapter"

gzip -dc $merged_file | cutadapt -u 3 -O 7 -m 30 -a $adapter_sequence --discard-untrimmed -o $trimmed_file - > $cutadapt_trimming_stats
cat $trimmed_file | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_ --suffix .fastq --eol > $bc_split_stats

echo "Removing adapter done. Outfiles: $demultiplexed_path/*.fastq. Stats file: $cutadapt_trimming_stats"
  
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
echo "Creating symlinks to $fastq_path"
# mv or cp?: ln -s!!
for f in $(ls $umi_extract_path/*fastq.gz); do
  fn=$(basename $f)
  fn=${fn%"_umi_extracted.fastq.gz"}
  fn=${fn#"dem_"}
  ln -s $f $fastq_path/${fn}.fastq.gz
done
echo "Done"
