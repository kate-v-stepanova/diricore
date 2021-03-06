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

subset=""
if [[ $# -ge 2 ]]; then
    subset=$2
fi

demultiplexed_path="$project_path/analysis/output/demultiplexed"
umi_extract_path="$project_path/analysis/output/umi_extract"
umi_extract_logs="${umi_extract_path}/logs"
fastq_path="$project_path/analysis/input/fastq"
diricore_path="/icgc/dkfzlsdf/analysis/OE0532/software/diricore"


bc_path="$project_path/analysis/input/metadata/"

cutadapt_trimming_stats="$project_path/analysis/output/cutadapt_trimming_stats.txt"
bc_split_stats="$project_path/analysis/output/bc_split_stats.txt"
stats_dir="$project_path/analysis/output"

adapter_sequence="AGATCGGAAGAGCACACGTCTGAA"
# adapter_sequence="TTCAGACGTGTGCTCTTCCGATCT" # reverse

bc_pattern="NNNNNNNNNN" # random sequence 5nt long + barcode
barcode_splitter="$diricore_path/programs/fastx_toolkit/fastx_barcode_splitter.pl"

# Creating directory structure (if not exist)
mkdir -p $demultiplexed_path
mkdir -p $umi_extract_logs
mkdir -p $fastq_path

# remove adapter
min_len=30
if [[ $# -ge 2 ]]; then
   min_len=$2
fi

for f in $(ls $merged_path/*$subset*.fastq.gz); do
    echo "Processing $f"
    fn=$(basename $f);
    fn=${fn%.fastq.gz}
    trimmed_file="$merged_path/${fn}_trimmed.fastq.gz"
    if [[ -f $trimmed_file ]]; then
        echo "File exists! $trimmed_file"
    else
        cutadapt_trimming_stats="$stats_dir/${fn}_cutadapt_trimming_stats.txt"
        bc_split_stats="$stats_dir/${fn}_bc_split_stats.txt"
        bc_file="$bc_path/bc_file_${fn}.txt"
    #    echo "gzip -dc $merged_file | cutadapt -u 3 -O 7 -m $min_len -j 10 -a $adapter_sequence --discard-untrimmed -o $trimmed_file - > $cutadapt_trimming_stats"
    #    gzip -dc $merged_file | cutadapt -u 3 -O 7 -m $min_len -j 10 -a $adapter_sequence --discard-untrimmed -o $trimmed_file - > $cutadapt_trimming_stats
        echo "Removing adapter: gzip -dc $f | cutadapt -u 3 -O 7 -m $min_len -j 10 -a $adapter_sequence --discard-untrimmed -o $trimmed_file - > $cutadapt_trimming_stats"

        gzip -dc $f | cutadapt -u 3 -O 5 -m $min_len -j 10 -a $adapter_sequence --discard-untrimmed -o $trimmed_file - > $cutadapt_trimming_stats
        echo "Removing adapter done. Outfiles: $trimmed_file. Stats file: $cutadapt_trimming_stats"
    fi
    echo "Demultiplexing: cat $trimmed_file | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_ --suffix .fastq --eol > $bc_split_stats"
    zcat $trimmed_file | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_  --suffix .fastq --eol > $bc_split_stats
    echo "Demultiplexing done"
    mv $demultiplexed_path/dem_unmatched.fastq $demultiplexed_path/dem_unmatched_${fn}.fastq
done
# cat $trimmed_file | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_ --suffix .fastq --eol > $bc_split_stats

# remove umi (any random sequence of 5-nt length)
echo "Removing UMIs"
for f in $demultiplexed_path/*.fastq; do
    fb=$(basename ${f});
  b=${fb%%.*};
  if [[ "$b" != "dem_unmatched*" ]]; then
      echo "Processing ${f}";
      umi_tools extract --extract-method=string --3prime --bc-pattern=$bc_pattern --log=$umi_extract_logs/${b}.log --log2stderr --stdin=$demultiplexed_path/${fb} --stdout=$umi_extract_path/${b}_umi_extracted.fastq.gz
      echo "Done ${f}"
  fi
done;

# the result of umi_extract will become input for the next step (preprocessing)
echo "Creating symlinks to $fastq_path"
for f in $(ls $umi_extract_path/*fastq.gz); do
  fn=$(basename $f)
  fn=${fn%"_umi_extracted.fastq.gz"}
  fn=${fn#"dem_"}
  ln -s $f $fastq_path/${fn}.fastq.gz
done
echo "Done"
