#!/bin/bash

set -e;
set -u;

dataset_id=$1
BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"

OUTDIR="$PROJECT_DIR/analysis/output/rpf_5p_density";
#INDIR="$PROJECT_DIR/analysis/output/tophat_out";
INDIR="$PROJECT_DIR/analysis/output/alignments/toGenome"
PLOTDIR="$PROJECT_DIR/analysis/output/figures";
SAMPLENAME_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv";
CONTRAST_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";

species=$2;
projectname=$dataset_id;
project_id=$dataset_id

if [[ $# -ge 3 ]]; then
  minreads=$3
else
  minreads=100
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
#  setup
of_hq_unique="${OUTDIR}/${projectname}.txcoord_counts.hq.dedup.${minreads}.hdf5";
of_hq="${OUTDIR}/${projectname}.txcoord_counts.hq.${minreads}.hdf5";
of_all_unique="${OUTDIR}/${projectname}.txcoord_counts.all.dedup.${minreads}.hdf5";
of_all="${OUTDIR}/${projectname}.txcoord_counts.all.${minreads}.hdf5";

bam_type="hq_unique"
# can be: hq, hq_unique, all, all_unique
if [[ $# -ge 4 ]]; then 
    bam_type=$4;
    echo "bam_type: $bam_type"
fi
echo "$bam_type"

plots_only=0
if [[ $# -ge 5 ]]; then
  plots_only=$5
  if [[ "$plots_only" == "plots_only" ]]; then
      plots_only=1
  fi
fi


mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}/rpf_5p_density_plots/"

script_dir="/home/e984a/faster_diricore/rpf_density/$project_id"
mkdir -p $script_dir

if [ "$bam_type" == "hq" ]; then
    outfile="$of_hq"
    bam_pattern="_toGenome.hqmapped.bam"
elif [ "$bam_type" == "all_unique" ]; then
    outfile=$of_all_unique
    bam_pattern="_toGenome_dedup.bam"
elif [ $bam_type == "all" ]; then
    outfile=$of_all
    bam_pattern="_toGenome.bam"
else
    outfile=$of_hq_unique
    bam_pattern="_toGenome.hqmapped_dedup.bam"
fi


if [ "$plots_only" -eq 0 ]; then
 # map RPFs to transcriptome coordinates
 echo "Mapping RPFs to transcripttome coordinates ($bam_type)"
 bam_files=$(ls ${INDIR}/*$bam_pattern);
 echo $bam_files;
 ii=$(for f in $bam_files; do b=$(basename $f); b=${b%_$bam_pattern}; echo -e "${b}\t$f"; done);
 $DIRICORE_DIR/diricore/bin/map_rpfs_to_transcriptome_positions.py  -t ${INDEXDATA_FILE}  -o ${outfile} -b $ii
 #echo "$DIRICORE_DIR/diricore/bin/map_rpfs_to_transcriptome_positions.py  -t ${INDEXDATA_FILE}  -o ${outfile} -b <(ls -1 ${INDIR}/*$bam_pattern | sort -V | while read fn; do b=$(basename $fn); b=${b%_$bam_pattern}; echo -e '${b}\t${fn}'; done )"
 echo "Mapping done. Created file: $outfile"
fi
exit

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

echo "Generating RPF density shift plots ($bam_type)"
mkdir -p ${PLOTDIR}/rpf_5p_density_plots/$bam_type
echo -ne "$output" | while read aa codongroupstr; do
    plot_file="${PLOTDIR}/rpf_5p_density_plots/$bam_type/${projectname}.${bam_type}.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    if [[ ! -f $plot_file ]]; then
      codons=$(echo $codongroupstr | sed 's/,/ /g');

      python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${plot_file}" \
        -m ${minreads} \
        "${outfile}" \
        ${MAPS_FILE} \
        ${codons}
      echo "Shift plots done ($bam_type). Created file: $plot_file"
    fi
done

# RPF density at special codons
# (START/other ATG)
echo "RPF density at special codons ($bam_type)"
echo -ne "\
ATG_split\tSTART_ATG,Other_ATG
" | while read aa codongroupstr; do
    plot_file="${PLOTDIR}/rpf_5p_density_plots/$bam_type/${projectname}.$bam_type.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    if [[ ! -f $plot_file ]]; then
      codons=$(echo $codongroupstr | sed 's/,/ /g');
      python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${plot_file}" \
        -m ${minreads} \
        "${outfile}" \
        ${MAPSSTART_FILE} \
        ${codons}
      echo "Special codons done ($bam_type). Created file: $plot_file"
    fi
done
