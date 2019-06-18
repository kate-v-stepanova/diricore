#!/bin/bash

set -e;
set -u;

dataset_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"

OUTDIR="$PROJECT_DIR/analysis/output/rpf_5p_density";
INDIR="$PROJECT_DIR/analysis/output/alignments/toGenome";
PLOTDIR="$PROJECT_DIR/analysis/output/figures";
SAMPLENAME_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv";
CONTRAST_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";

species=$2;
projectname=$dataset_id;

if [[ $# -ge 3 ]]; then
  minreads=$3
else
  minreads=100
fi

plots_only=0
if [[ $# -ge 4 ]]; then
  plots_only=$4
  if [[ "$plots_only" == "plots_only" ]]; then
      plots_only=1
  fi
fi

DIRICORE_DIR="/home/e984a/diricore"
INDEXDATA_FILE="$BASE_DIR/static/${species}/transcript_data.hdf5";
MAPS_FILE="$BASE_DIR/static/${species}/codon_regions.width_61.hdf5";
MAPSSTART_FILE="$BASE_DIR/static/${species}/codon_regions.START_Other_ATG.width_61.hdf5";

# Checking that files are present
ls $INDEXDATA_FILE
ls $MAPS_FILE
ls $MAPSSTART_FILE

###
# { setup
of_hq_unique="${OUTDIR}/${projectname}.txcoord_counts.hq.dedup.${minreads}.hdf5";
of_hq="${OUTDIR}/${projectname}.txcoord_counts.hq.${minreads}.hdf5";
of_all_unique="${OUTDIR}/${projectname}.txcoord_counts.all.dedup.${minreads}.hdf5";
of_all="${OUTDIR}/${projectname}.txcoord_counts.all.${minreads}.hdf5";


mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}/rpf_5p_density_plots/"

if [[ $plots_only == 0 ]]; then
  if [[ ! -f $of_hq ]]; then
    echo "Mapping RPFs to transcripttome coordinates (hq)"
    tmp_file="$OUTDIR/tmp_file.${minreads}.${species}.txt"
    rm -f $tmp_file
    touch $tmp_file
    for fn in $(ls ${INDIR}/*.hqmapped.bam); do
        b=$(basename $fn);
        b=${b#"star_"};
        b=${b%%"_toGenome.hqmapped.bam"};
        echo -e "${b}\t${fn}" >> $tmp_file
    done;

    $DIRICORE_DIR/diricore/bin/map_rpfs_to_transcriptome_positions.py \
       -t "${INDEXDATA_FILE}" \
       -o "${of_hq}" \
       -b $tmp_file

    echo "Mapping done. Created file: $of_hq"
  fi
  
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

# RPF density at special codons
# (START/other ATG)
if [[ -f $of_hq ]]; then
  echo "Generating RPF density shift plots (hq)"
  echo "Data file: $of_hq"
  mkdir -p $PLOTDIR/rpf_5p_density_plots/hq_with_dup
  echo -ne "$output" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/hq_with_dup/${projectname}.hq.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${of_hq}" \
        ${MAPS_FILE} \
        ${codons}
    echo "Shift plots done (hq). Created file: $of"
  done

  # RPF density at special codons
  # (START/other ATG)
  echo "RPF density at special codons (hq)"
  echo -ne "\
  ATG_split\tSTART_ATG,Other_ATG
  " | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/hq_with_dup/${projectname}.hq.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${of_hq}" \
        ${MAPSSTART_FILE} \
        ${codons}
    echo "Special codons done (hq). Created file: $of"
  done
else
  echo "File is missing! Skipping: $of_hq"
fi
