#!/bin/bash

set -e
set -u

dataset_id=$1

# trim and demultiplex
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/output/demultiplexed


mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/fastq
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/merged
touch /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata/bc_file.txt
touch /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata/rpf_density_contrasts.tsv
touch /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata/rpf_density_samplenames.tsv

mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/output/demultiplexed
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/output/trimmed
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/output/umi_extract/logs
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/output/figures/
