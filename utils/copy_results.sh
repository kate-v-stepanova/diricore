#!/bin/bash

set -e
set -u

dataset_id=$1

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
project_dir="$BASE_DIR/$dataset_id"
cluster="odcf-lsf02.dkfz.de"
username="e984a"

out_dir="/Users/b250-admin/analysis/$dataset_id"

scp -r $username@$cluster:$project_dir/analysis/output/figures $out_dir
