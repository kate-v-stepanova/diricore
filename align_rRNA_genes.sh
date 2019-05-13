#!/bin/bash


project_id=$1
DIRICORE_PATH="/home/e984a/diricore"
blat_path="$DIRICORE_PATH/programs/blat_for_linux"
BASE_PATH="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_PATH="$BASE_PATH/$project_id"
PSL_PARSER="$DIRICORE_PATH/utils/parse_psl.py"

genome="hg19"
if [[ $# -ge 3 ]]; then
   genome=$3
fi

GENOME_PATH="$BASE_PATH/static/$genome/rRNA_genes.2bit"


trna=0
prefix="rRNA"
if [[ $# -ge 2 ]]; then
   trna=$2
   if [[ $trna == 'trna' ]]; then
       prefix="tRNA"
       INDIR="$PROJECT_PATH/analysis/output/trna_fragments"
   else
      prefix="rRNA"
      INDIR="$PROJECT_PATH/analysis/output/rrna_fragments"
   fi
fi

for fasta_file in $(ls $INDIR/*.fasta); do
    sample_name=$(basename $fasta_file)
    sample_name=${sample_name%top_${prefix}_seqs.fasta}
    psl_file="$INDIR/${sample_name}.${prefix}_aligned.with_header.psl"
    psl_no_header="$INDIR/${sample_name}.${prefix}_aligned.psl"
    echo "Alignment of $fasta_file";
    $blat_path -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 $GENOME_PATH $fasta_file $psl_file
    # remove first 5 lines (header)
    tail -n +6 $psl_file > $psl_no_header
    rm $psl_file
    echo "Alignment done. Processing $psl_no_header"
    output_file="$INDIR/${sample_name}.parsed_psl.txt"
    $PSL_PARSER $psl_no_header > $output_file
    echo "Created $output_file"
done;
