#!/bin/bash

set -e;
set -u;

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
INDEXDATA_FILE="$DIRICORE_DIR/staticdata/${species}/transcript_data.hdf5";
MAPS_FILE="$DIRICORE_DIR/staticdata/${species}/codon_regions.width_61.hdf5";
MAPSSTART_FILE="$DIRICORE_DIR/staticdata/${species}/codon_regions.START_Other_ATG.width_61.hdf5";

# Checking that files are present
ls $INDEXDATA_FILE
ls $MAPS_FILE
ls $MAPSSTART_FILE

###
# { setup
of_hq_unique="${OUTDIR}/${projectname}.txcoord_counts.hq.dedup.hdf5";
of_hq="${OUTDIR}/${projectname}.txcoord_counts.hq.hdf5";
of_all_unique="${OUTDIR}/${projectname}.txcoord_counts.all.dedup.hdf5";
of_all="${OUTDIR}/${projectname}.txcoord_counts.all.hdf5";


mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}/rpf_5p_density_plots/"

if [[ $plots_only == 0 ]]; then
# map RPFs to transcriptome coordinates
echo "Mapping RPFs to transcripttome coordinates (hq unique)"
$DIRICORE_DIR/diricore/bin/map_rpfs_to_transcriptome_positions.py \
   -t "${INDEXDATA_FILE}" \
   -o "${of_hq_unique}" \
   -b <(\
        ls -1 ${INDIR}/*/accepted_hits.hqmapped_dedup.bam | sort -V | while read fn; do
               b=$(basename $(dirname "$fn"));
               b=${b%%.*};
               b=${b#"dem_"};
               b=${b%"_umi_extracted"};
               echo -e "${b}\t${fn}";
           done \
       )
echo "Mapping done. Created file: $of_hq_unique"

echo "Mapping RPFs to transcripttome coordinates (hq)"
$DIRICORE_DIR/diricore/bin/map_rpfs_to_transcriptome_positions.py \
   -t "${INDEXDATA_FILE}" \
   -o "${of_hq}" \
   -b <(\
        ls -1 ${INDIR}/*/accepted_hits.hqmapped.bam | sort -V | while read fn; do
               b=$(basename $(dirname "$fn"));
               b=${b%%.*};
               b=${b#"dem_"};
               b=${b%"_umi_extracted"};
               echo -e "${b}\t${fn}";
           done \
       )
echo "Mapping done. Created file: $of_hq"

echo "Mapping RPFs to transcripttome coordinates (all)"
$DIRICORE_DIR/diricore/bin/map_rpfs_to_transcriptome_positions.py \
   -t "${INDEXDATA_FILE}" \
   -o "${of_all}" \
   -b <(\
        ls -1 ${INDIR}/*/accepted_hits.bam | sort -V | while read fn; do
               b=$(basename $(dirname "$fn"));
               b=${b%%.*};
               b=${b#"dem_"};
               b=${b%"_umi_extracted"};
               echo -e "${b}\t${fn}";
           done \
       )
echo "Mapping done. Created file: $of_all"

echo "Mapping RPFs to transcripttome coordinates (all, unique)"
$DIRICORE_DIR/diricore/bin/map_rpfs_to_transcriptome_positions.py \
   -t "${INDEXDATA_FILE}" \
   -o "${of_all_unique}" \
   -b <(\
        ls -1 ${INDIR}/*/accepted_hits_dedup.bam | sort -V | while read fn; do
               b=$(basename $(dirname "$fn"));
               b=${b%%.*};
               b=${b#"dem_"};
               b=${b%"_umi_extracted"};
               echo -e "${b}\t${fn}";
           done \
       )
echo "Mapping done. Created file: $of_all_unique"
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
echo "Generating RPF density shift plots (hq, unique)"
mkdir -p ${PLOTDIR}/rpf_5p_density_plots/hq_unique
echo -ne "$output" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/hq_unique/${projectname}.hq.unique.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${of_hq_unique}" \
        ${MAPS_FILE} \
        ${codons}
    echo "Shift plots done (hq, unique). Created file: $of"
done


# RPF density at special codons
# (START/other ATG)
echo "RPF density at special codons (hq, unique)"
echo -ne "\
ATG_split\tSTART_ATG,Other_ATG
" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/hq_unique/${projectname}.hq.unique.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${of_hq_unique}" \
        ${MAPSSTART_FILE} \
        ${codons}
    echo "Special codons done (hq, unique). Created file: $of"
done


echo "Generating RPF density shift plots (hq)"
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


echo "Generating RPF density shift plots (all)"
mkdir -p $PLOTDIR/rpf_5p_density_plots/all_with_dup
echo -ne "$output" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/all_with_dup/${projectname}.all.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${of_all}" \
        ${MAPS_FILE} \
        ${codons}
echo "Shift plots done (all). Created file: $of"
done


# RPF density at special codons
# (START/other ATG)
echo "RPF density at special codons (all)"
echo -ne "\
ATG_split\tSTART_ATG,Other_ATG
" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/all_with_dup/${projectname}.all.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${of_all}" \
        ${MAPSSTART_FILE} \
        ${codons}
    echo "Special codons done (all). Created file: $of"
done


echo "Generating RPF density shift plots (all, unique)"
mkdir -p $PLOTDIR/rpf_5p_density_plots/all_unique
echo -ne "$output" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/all_unique/${projectname}.all.unique.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${of_all_unique}" \
        ${MAPS_FILE} \
        ${codons}
    echo "Shift plots done (all, unique). Created file: $of"
done

# RPF density at special codons
# (START/other ATG)
echo "RPF density at special codons (all, unique)"
echo -ne "\
ATG_split\tSTART_ATG,Other_ATG
" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/all_unique/${projectname}.all.unique.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python $DIRICORE_DIR/diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${of_all_unique}" \
        ${MAPSSTART_FILE} \
        ${codons}
    echo "Special codons done (all, unique). Created file: $of"
done
