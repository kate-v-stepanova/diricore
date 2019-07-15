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
trimmed_path="$project_path/analysis/output/trimmed"
umi_extract_path="$project_path/analysis/output/umi_extract"
umi_extract_logs="${umi_extract_path}/logs"
fastq_path="$project_path/analysis/input/fastq"
diricore_path="/home/e984a/diricore"


bc_file="$project_path/analysis/input/metadata/bc_file.txt"
cutadapt_trimming_stats="$project_path/analysis/output/cutadapt_trimming_stats.txt"
bc_split_stats="$project_path/analysis/output/bc_split_stats.txt"

adapter_sequence="AGATCGGAAGAGCACACGTCTGAA"
bc_pattern="NNNNN" # random sequence 5nt long + barcode
barcode_splitter="$diricore_path/programs/fastx_toolkit/fastx_barcode_splitter.pl"

# Creating directory structure (if not exist)
mkdir -p $demultiplexed_path
mkdir -p $umi_extract_logs
mkdir -p $fastq_path

echo "Removing adapter"
# remove adapter based on this:
# first 5 nt is a barcode, the next ones are the first nts of the adapter
#cutadapt -a P2_1_ATCGT_SF050_LARS=ATCGTAGATCGGAA -a P2_2_AGCTA_SF051_LARS=AGCTAAGATCGGAA -a P2_3_CGTAA_SF054_HC=CGTAAAGATCGGAA -a P2_4_CTAGA_SF054_HC_noco=CTAGAAGATCGGAA -a P2_5_GATCA_SF055_HC=GATCAAGATCGG -m 20 --discard-untrimmed --times=2 -e 0 -o $demultiplexed_path/{name}_trimmed.fastq.gz $merged_file > $cutadapt_trimming_stats
part_of_adapter=$(echo $adapter_sequence | cut -c1-5 )
echo $part_of_adapter
command="cutadapt "
while read -r line; do
    sample=$(echo $line | cut -d ' ' -f1)
    barcode=$(echo $line | cut -d ' ' -f2)
    to_append="-a $sample=${barcode}${part_of_adapter}"
    command="${command} $to_append"
done < $bc_file

command="$command -m 20 --untrimmed-output $demultiplexed_path/untrimmed.fastq.gz --times=2 -e 0 -o $demultiplexed_path/{name}_trimmed.fastq.gz $merged_file"
echo $command
a=$($command)

## gzip -dc $merged_file | cutadapt -u 3 -O 7 -m 20 -a $adapter_sequence --discard-untrimmed -o $trimmed_file - > $cutadapt_trimming_stats
## cat $trimmed_file | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_ --suffix .fastq --eol > $bc_split_stats

echo "Removing adapter done. Outfiles: $demultiplexed_path/*.fastq.gz. Stats file: $cutadapt_trimming_stats"
  
# remove umi (any random sequence of 5-nt length)
echo "Removing UMIs"
for f in $(ls $demultiplexed_path/*.fastq.gz); do
   
    fb=$(basename ${f});
    b=${fb%%.*};
    #echo $b
#  if [[ "$b" != "dem_unmatched" ]]; then
     #echo "Processing ${f}";
     umi_tools extract --extract-method=string --3prime --bc-pattern=$bc_pattern --log=$umi_extract_logs/${b}.log --log2stderr --stdin=$demultiplexed_path/${fb} --stdout=$umi_extract_path/${b}_umi_extracted.fastq.gz
      # echo "Done ${f}"
 # fi
done;

