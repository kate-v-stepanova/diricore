#!/bin/bash

set -e;
set -u;

project_id=$1
species=$2
minreads=$3

plots_only=0
if [[ $# -ge 4 ]]; then
  plots_only=$4
  if [[ "$plots_only" == "plots_only" ]]; then
      plots_only=1
  fi
fi

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$project_id"
OUTDIR="$PROJECT_DIR/analysis/output/rpf_5p_density";
INDIR="$PROJECT_DIR/analysis/output/tophat_out";
PLOTDIR="$PROJECT_DIR/analysis/output/figures/rpf_density";
SAMPLENAME_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv";
CONTRAST_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";

DIRICORE_PATH="/home/e984a/diricore"
INDEXDATA_FILE="$DIRICORE_PATH/staticdata/${species}/transcript_data.hdf5";
MAPS_FILE="$DIRICORE_PATH/staticdata/${species}/codon_regions.width_61.hdf5";
MAPSSTART_FILE="$DIRICORE_PATH/staticdata/${species}/codon_regions.START_Other_ATG.width_61.hdf5";


###
# { setup
of="${OUTDIR}/${project_id}.txcoord_counts.hdf5";

mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}"

if [[ $plots_only == 0 ]]; then
  echo "Mapping RPFs to transcriptome coordinates"
  # map RPFs to transcriptome coordinates
  $DIRICORE_PATH/diricore/bin/map_rpfs_to_transcriptome_positions.py \
    -t "${INDEXDATA_FILE}" \
    -o "${of}" \
    -b <(\
        ls -1 ${INDIR}/*.bam | sort -V | while read fn; do
                b=$(basename "$fn");
                b=${b%%.*};

                echo -e "${b}\t${fn}";
            done \
        )
  # }
  echo "Done. Created file: $of"
fi

# { generate RPF density shift plots
echo "Generating RPF density shift plots"
output="\
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
"

datafile="${of}";

echo -ne "$output" | while read aa codongroupstr; do
    of="${PLOTDIR}/${project_id}.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');
    python $DIRICORE_PATH/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
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
    of="${PLOTDIR}/${project_id}.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_PATH/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${datafile}" \
        ${MAPSSTART_FILE} \
        ${codons}
done
