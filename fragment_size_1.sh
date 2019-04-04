#!/bin/bash

set -e
set -u

dataset_id=$1

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
OUTDIR="$PROJECT_DIR/analysis/output/fragment_size"
SEQTK_DIR="/home/e984a/seqtk"

mkdir -p $OUTDIR

for f in $(ls ${PROJECT_DIR}/analysis/input/fastq/*.fastq.gz); do
  echo "Processing ${f}";
  fb=$(basename ${f});
  b="${fb%%.*}";
  tmp_file="${b}_tmp.fastq"
  gzip -dc $f | cutadapt -u 3 -O 5 -a AGATCGGAAGAGCACACGTCTGAA --discard-untrimmed -  > ${OUTDIR}/${tmp_file}
  echo "Extracting fragment size from ${b}"
  ${SEQTK_DIR}/seqtk seq -a ${OUTDIR}/${tmp_file} | grep -v '^>' | awk '{print length}' | sort | uniq -c | sed -e 's/^ \+//g' | tr ' ' '\t' > ${OUTDIR}/${b}_UMIs_fragment_size.txt
  echo "Written file: ${OUTDIR}/${b}_UMIs_fragment_size.txt"
  echo "Removing tmp file ${tmp_file}"
 #  rm ${OUTDIR}/${tmp_file}
  # trap "rm -f ${OUTDIR}/${tmp_file}" EXIT
done;

echo "Done"
