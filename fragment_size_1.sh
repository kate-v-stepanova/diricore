#!/bin/bash

set -e
set -u

OUTDIR="./data/output/fragment_size"
SEQTK_DIR="/home/e984a/seqtk"
for f in $(ls ./data/input/fastq/run_umi*.fastq.gz); do
 fb=$(basename ${f})
 b=${fb%%.*};
 tmp_file=${b}_tmp.fastq
gzip -dc $f | cutadapt -u 3 -O 5 -a AGATCGGAAGAGCACACGTCTGAAX --discard-untrimmed - 2> /dev/null > ${OUTDIR}/${tmp_file}

${SEQTK_DIR}/seqtk seq -a ${OUTDIR}/${tmp_file} | grep -v '^>' | awk '{print length}' | sort | uniq -c | sed -e 's/^ \+//g' | tr ' ' '\t' > ${OUTDIR}/${b}_UMIs_fragment_size.txt
rm ${OUTDIR}/${tmp_file}
trap "rm -f ${OUTDIR}/${tmp_file}" EXIT

#  gzip -dc $f | cutadapt -u 3 -O 5 -a AGATCGGAAGAGCACACGTCTGAAX --discard-untrimmed - 2> /dev/null | fastq_to_fasta -Q33 | grep -v '^>' | awk '{print length}' | sort | uniq -c | sed -e 's/^ \+//g' | tr ' ' '\t' > ./data/output/fragment_size/${fb}_UMIs_fragment_size.txt
done;
#gzip -dc Data/runUmis.fastq.gz | cutadapt -u 3 -O 5 -a AGATCGGAAGAGCACACGTCTGAAX --discard-untrimmed - 2> /dev/null | fastq_to_fasta -Q33 | grep -v '^>' | awk '{print length}' | sort | uniq -c | sed -e 's/^ \+//g' | tr ' ' '\t' > run_UMIs_fragment_sizes.txt
