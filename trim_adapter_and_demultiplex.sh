#!/bin/bash

# exit on error:
set -e
# exit on Ctrl-C
trap 'exit' INT

flowcell_path=$1

# merge fastq files into one
echo "Merging $flowcell_path/*/*/*.fastq.gz into one file"
cat $flowcell_path/*/*/*.fastq.gz > ./data/input/merged/merged.fastq.gz
echo "Merging done"

# merge fastq files into one
cat ./data/input/fastq/*.fastq.gz > ./data/input/merged/merged.fastq.gz
rm -f ./data/input/fastq/*.fastq.gz

# remove adapter
gzip -dc data/input/merged/*merged.fastq.gz | cutadapt -u 3 -O 7 -m 30 -a AGATCGGAAGAGCACACGTCTGAA --discard-untrimmed - 2> data/output/cutadapt.log | ../fastx_toolkit/fastx_barcode_splitter.pl --bcfile data/input/metadata/bc_file.txt --prefix data/output/demultiplexed/dem_ --suffix .fastq --eol 2> data/output/bc_split.log

# remove umi (any random sequence of 5-nt length)
mkdir -p ./data/output/umi_extract
mkdir -p ./data/output/umi_extract/logs
for f in ./data/output/demultiplexed/*.fastq; do
  echo "Processing ${f}";
  fb=$(basename ${f});
  b=${fb%%.*};
  umi_tools extract --extract-method=string --3prime --bc-pattern=NNNNN --log=data/output/umi_extract/logs/${b}.log --log2stderr --stdin=data/output/demultiplexed/${fb} --stdout=data/output/umi_extract/${b}_umi_extracted.fastq.gz
  echo "Done ${f}"
done;

# the result of umi_extract will become input for the next step (preprocessing)
echo "Moving files to data/input/fastq"
# mv or cp?
mv data/output/umi_extract/*.fastq.gz data/input/fastq
echo "Done"

