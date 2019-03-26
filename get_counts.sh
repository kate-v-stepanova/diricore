#!/bin/bash

set -e;
set -u;

dataset_id=$1

DIRICORE_PATH="/home/e984a/diricore"
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
counts_dir="$PROJECT_DIR/analysis/output/counts"

mkdir -p $counts_dir

gtf="$DIRICORE_PATH/staticdata/human/gencode.v29.basic.annotation.gtf.gz"

ls $gtf


files=$(ls -1 $PROJECT_DIR/analysis/output/tophat_out/*/accepted_hits.hqmapped_dedup.bam | tr '\n' ' ')
htseq-count -f bam -r pos -s yes -t exon -i gene_id --additional-attr=gene_name --additional-attr=gene_type -m union $files $gtf > $counts_dir/htseq_counts_dedup.txt

# Add the header (NOTE the double quotes to pass var to sed)
# Note, the first step removes a trailing space
files=$(echo $files | sed 's/[[:blank:]]+$//')
header="Gene_id\tGene_name\tGene_type\t${files// /\\t}"
sed -i "1i$header" $counts_dir/htseq_counts_dedup.txt

# repeat the same steps using all the HQ algs (before deduping)
files=$(ls -1 $PROJECT_DIR/analysis/output/tophat_out/*/accepted_hits.hqmapped.bam | tr '\n' ' ')
htseq-count -f bam -r pos -s yes -t exon -i gene_id --additional-attr=gene_name --additional-attr=gene_type -m union $files $gtf > $counts_dir/htseq_counts.txt

# Add the header (NOTE the double quotes to pass var to sed)
# Note, the first step removes a trailing space
files=$(echo $files | sed 's/[[:blank:]]+$//')
header="Gene_id\tGene_name\tGene_type\t${files// /\\t}"
sed -i "1i$header" $counts_dir/htseq_counts.txt


