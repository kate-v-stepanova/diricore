#!/bin/bash


. ./diricore_virtualenv/bin/activate

set -e;
set -u;


export BOWTIE_PATH="./programs/bowtie2-2.0.6/";
export TOPHAT_BIN="./programs/tophat-2.0.7/tophat2";

export INDIR="./data/output/clean/";
export OUTDIR="./data/output/tophat_out/";

echo "settings:";
export species=${1:-human}; echo -e "\tspecies: $species";
export cores=${2:-1}; echo -e "\tparallel processes: $cores";
export threads=${3:-2}; echo -e "\ttophat threads: $threads";

export REF="./staticdata/${species}/genome";
export GTF="./staticdata/${species}/transcripts.gff";
export TIDX="./staticdata/${species}/transcripts";

echo "Doing alignments in parallel..."
##
run_align(){
	fn=$1
	bn=$(basename "$fn");
    	b=${bn%%.*};
	od="${OUTDIR}/${b}";
    	mkdir -p "${od}" 2> /dev/null || true;
   	thout="${OUTDIR}/${b}/tophat.out"
    	therr="${OUTDIR}/${b}/tophat.err"
	echo "Starting alignment of file: $bn";
  		 PATH="${BOWTIE_PATH}:$PATH" \
       		 ${TOPHAT_BIN} --seed 42 -n 2 -m 1 -p "${threads}"\
       		 --no-novel-juncs --no-novel-indels --no-coverage-search \
       		 --segment-length 25 \
       		 --transcriptome-index "${TIDX}" -G "${GTF}" \
       		 -o "${od}" \
       		 "${REF}" "${fn}" \
       		 > "${thout}" 2> "${therr}";
	echo "Done with file: $bn";
}

export -f run_align
ls -1 ${INDIR}/*.fastq.gz | parallel -j $cores run_align

echo "Done with alignments; now filtering out non-primary/low-quality alignments...";

# isolate "hqmapped" reads
run_filter(){
	samtools_bin="./programs/tophat-2.0.7/samtools"
	fn=$1
	dn=$(dirname "$fn");
    	of="${dn}/accepted_hits.hqmapped.bam";
	cat \
        <(${samtools_bin} view -H "${fn}") \
        <(cat "${fn}" | ${samtools_bin} view -q10 -F260 -) \
	| ${samtools_bin} view -bS - \
   	> "${of}";
}
export -f run_filter
ls -1 ${OUTDIR}/*/accepted_hits.bam | parallel -j $cores run_filter

echo "Done filtering";
###

