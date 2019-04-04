#!/bin/bash

set -e
set -u

if [[ "$#" -eq 1 ]]; then
    dataset_id=$1
else
    echo "Usage example: $0 14548
	where 14548 is dataset_id" 
    exit 
fi;


analysis_path="/icgc/dkfzlsdf/analysis/OE0532"
mkdir -p $analysis_path/$dataset_id

# cp -R /midterm/$dataset_id/data/ $analysis_path/$dataset_id
sftp -r ad+e984a@ftp4midterm.dkfz.de:0$dataset_id/data/* $analysis_path/$dataset_id

# Create directory structure
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/fastq
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/metadata
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/input/merged

mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/output/demultiplexed
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/$dataset_id/analysis/output/umi_extract/logs


