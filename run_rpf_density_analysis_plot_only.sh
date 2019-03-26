#!/bin/bash

set -e;
set -u;

# This modification disable the mapping to transcriptome coordinates so the script can be used to make plots from 
# previously analyzed data. 
# IMPORTANT: project name and sample names must be exactly as used when the complete script was run

dataset_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"

OUTDIR="$PROJECT_DIR/analysis/output/rpf_5p_density";
INDIR="$PROJECT_DIR/analysis/output/tophat_out";
PLOTDIR="$PROJECT_DIR/analysis/output/figures";
SAMPLENAME_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv";
CONTRAST_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";

species=$2;
projectname=$dataset_id;
minreads=$3;

DIRICORE_DIR="/home/e984a/diricore"
INDEXDATA_FILE="$DIRICORE_DIR/staticdata/${species}/transcript_data.hdf5";
MAPS_FILE="$DIRICORE_DIR/staticdata/${species}/codon_regions.width_61.hdf5";
MAPSSTART_FILE="$DIRICORE_DIR/staticdata/${species}/codon_regions.START_Other_ATG.width_61.hdf5";


###
# { setup
of="${OUTDIR}/${projectname}.txcoord_counts.hdf5";

mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}/rpf_5p_density_plots/"

# map RPFs to transcriptome coordinates
#./diricore/bin/map_rpfs_to_transcriptome_positions.py \
#    -t "${INDEXDATA_FILE}" \
#    -o "${of}" \
#    -b <(\
#        ls -1 ${INDIR}/*/accepted_hits.hqmapped.bam | sort -V | while read fn; do
#                b=$(basename $(dirname "$fn"));
#                b=${b%%.*};
#
#                echo -e "${b}\t${fn}";
#            done \
#        )
# }

# { generate RPF density shift plots
datafile="${of}";
echo -ne "\
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
" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/${projectname}.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
	--y-limits " -1,1" \
        -m ${minreads} \
        "${datafile}" \
        ${MAPS_FILE} \
        ${codons}
done

# RPF density at special codons
# (START/other ATG)
echo -ne "\
ATG_split\tSTART_ATG,Other_ATG
" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/${projectname}.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        --y-limits " -1,1" \
        -m ${minreads} \
        "${datafile}" \
        ${MAPSSTART_FILE} \
        ${codons}
done
