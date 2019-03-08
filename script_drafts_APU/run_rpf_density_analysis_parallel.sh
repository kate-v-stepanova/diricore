#!/bin/bash

set -e;
set -u;


export OUTDIR="./data/output/rpf_5p_density/";
export INDIR="./data/output/tophat_out/";
export PLOTDIR="./data/output/figures/";
export SAMPLENAME_FILE="./data/input/metadata/rpf_density_samplenames.tsv";
export CONTRAST_FILE="./data/input/metadata/rpf_density_contrasts.tsv";

export species=$1;
export projectname=$2;
export minreads=$3;

export INDEXDATA_FILE="./staticdata/${species}/transcript_data.hdf5";
export MAPS_FILE="./staticdata/${species}/codon_regions.width_61.hdf5";
export MAPSSTART_FILE="./staticdata/${species}/codon_regions.START_Other_ATG.width_61.hdf5";


###
# { setup
export of="${OUTDIR}/${projectname}.txcoord_counts.hdf5";

mkdir ${OUTDIR} || true;
mkdir -p "${PLOTDIR}/rpf_5p_density_plots/" || true;

# map RPFs to transcriptome coordinates
run_map(){
	inf=$1
	echo "Processing file: ;${inf}..."
	bname=$(basename $(dirname $inf))
	bname="${bname%%.*}"
	echo ./diricore/bin/map_rpfs_to_transcriptome_positions_APU.py \
	    -t "${INDEXDATA_FILE}" \
	    -o "${of}" \
	    -b "${inf}" \
	    -n "${bname}"
	echo "Done with file: ${inf}..."
}
export -f run_map
ls -1 ./data/output/tophat_out/*/accepted_hits.hqmapped.bam | sort -V | parallel run_map

# }

# { generate RPF density shift plots
export datafile="${of}";

plot_rpf(){
	inf=$1
	aa=$(echo -e $inf | cut -f1)
	codons=$(echo -e $inf | cut -f2 | sed 's/,/ /g')
	outf="${PLOTDIR}/rpf_5p_density_plots/${projectname}.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf"
	echo python ./diricore/bin/plot_rpf_5p_density.py \
		-c "${CONTRAST_FILE}" \
		-n "${SAMPLENAME_FILE}" \
		-o "${outf}" \
		-m "${minreads}" \
		"${datafile}" \
		"${MAPS_FILE}" \
		${codons}
}
export -f plot_rpf

echo -n "\
Ala\tGCA,GCC,GCG,GCT
Arg\tCGA,CGC,CGG,CGT,AGA,AGG
Asn\tAAC,AAT
Asp\tGAC,GAT
Cys\tTGC,TGT
Gln\tCAA,CAG
Glu\tGAA,GAG
Gly\tGGA,GGC,GGG,GGT
His\tCAC,CAT
Ile\tATA,ATC,ATT
Leu\tCTA,CTC,CTG,CTT,TTA,TTG
Lys\tAAA,AAG
Met\tATG
Phe\tTTC,TTT
Pro\tCCA,CCC,CCG,CCT
Ser\tTCA,TCC,TCG,TCT,AGC,AGT
Thr\tACA,ACC,ACG,ACT
Trp\tTGG
Tyr\tTAC,TAT
Val\tGTA,GTC,GTG,GTT
" | parallel plot_rpf 


# RPF density at special codons
# (START/other ATG)
echo -n "\
ATG_split\tSTART_ATG,Other_ATG\
" | parallel plot_rpf


