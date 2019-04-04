#!/bin/bash

set -e
set -u


dataset_id=$1
flowcell_id=$2
species=$3
minreads=$4

echo "Checking required files"
ls /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata/rpf_density_contrasts.tsv
ls /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata/rpf_density_samplenames.tsv
# ls /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata/rpf_transcript_distribution_sampleinfo.tsv
ls /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata/bc_file.txt
#ls /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata/subsequence_contrasts.tsv
# ls /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata/subsequence_samplenames.tsv

ls  /home/e984a/diricore/staticdata/human/gencode.v29.basic.annotation.gtf.gz

echo "./trim_adapter_and_demultiplex.sh $dataset_id $flowcell_id"
./trim_adapter_and_demultiplex.sh $dataset_id $flowcell_id

echo "./run_umi_preprocessing.sh $species $dataset_id"
./run_umi_preprocessing.sh $species $dataset_id

echo "./run_alignment.sh $species $dataset_id"
./run_alignment.sh $species $dataset_id

echo "./deduplicate_umi.sh $dataset_id"
./deduplicate_umi.sh $dataset_id

echo "Alignment stats"
echo "./get_alignment_stats.sh $dataset_id"
./get_alignment_stats.sh $dataset_id
echo "./get_RNAcleup_info.r $dataset_id"
./get_RNAcleup_info.r $dataset_id
echo "./plot_alignment_stats.r $dataset_id"
./plot_alignment_stats.r $dataset_id


echo "./get_counts.sh $dataset_id"
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/output/counts
./get_counts.sh $dataset_id

echo "./analyze_counts.r $dataset_id"
module load R
./analyze_counts.r $dataset_id

echo "./fragment_size_1.sh $dataset_id"
./fragment_size_1.sh $dataset_id

echo "./fragment_size_2.r $dataset_id"
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/output/figures/fragment_size
./fragment_size_2.r $dataset_id

echo "./generate_tracks.sh $dataset_id"
./generate_tracks.sh $dataset_id
echo "./generate_trackDB.sh $dataset_id"
./generate_trackDB.sh $dataset_id


echo "./find_top_RNA_fragments_1.sh $dataset_id"
./find_top_RNA_fragments_1.sh $dataset_id
echo "./find_top_RNA_fragments_2.r $dataset_id"
./find_top_RNA_fragments_2.r $dataset_id

echo "./run_umi_rpf_density_analysis.sh $dataset_id $species $minreads"
./run_umi_rpf_density_analysis.sh $dataset_id $species $minreads
echo "./run_umi_subsequence_analysis.sh $dataset_id $species $minreads"
./run_umi_subsequence_analysis.sh $dataset_id $species $minreads

echo "./plot_rpf_transcript_distribution.sh $dataset_id $minreads"
./plot_rpf_transcript_distribution.sh $dataset_id $minreads


