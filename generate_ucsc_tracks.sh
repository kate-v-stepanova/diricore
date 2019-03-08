#!/bin/bash

set -e;
set -u;

export project=$1;
export cores=${2:-1};
export OUTDIR="./data/output/ucsc_tracks/";
export INDIR="./data/output/tophat_out/";
export METAFILE="./data/input/metadata/ucsc_track_samplenames.tsv";
export bedtools_bin="./programs/bedtools"


mkdir ${OUTDIR} || true;



bamtobed(){
	sampname=$(echo $1 | cut -d ' ' -f1)
	samplabel=$(echo $1 | cut -d ' ' -f2)
	sampcol=$(echo $1 | cut -d ' ' -f3)
	samppath="${INDIR}/${sampname}/accepted_hits.hqmapped.bam"
	for strand in "+" "-"; do
	trackopt="name=${samplabel}_${strand} visibility=2 color=${sampcol} maxHeightPixels=50"
	outf="${OUTDIR}/${project}_${sampname}_${samplabel}_${strand}.bdg.gz"
	${bedtools_bin} genomecov \
	-ibam ${samppath} -bg -split \
	-strand "${strand}" \
	-trackline -trackopts "${trackopt}" \
	| gzip > ${outf}
	echo "Done with file: ${outf}"
	done
	}
export -f bamtobed

echo "Preparing UCSC tracks..."
cat $METAFILE | parallel bamtobed
echo "Done"
