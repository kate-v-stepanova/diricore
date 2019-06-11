#!/bin/bash


set -e;
set -u;

dataset_id=$1
species=$2;
projectname=$dataset_id

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

y_limits=1
if [[ $# -ge 5 ]]; then
    y_limits=$5
fi

echo "$y_limits"

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"

OUTDIR="$PROJECT_DIR/analysis/output/subsequence_data";
INDIR="$PROJECT_DIR/analysis/output/alignments/toGenome"
PLOTDIR="$PROJECT_DIR/analysis/output/figures";
SAMPLENAME_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv";
CONTRAST_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";

DIRICORE_DIR="/home/e984a/diricore"
INDEXDATAFN="$BASE_DIR/static/${species}/subseq_index_data.pkl.gz";
of="${OUTDIR}/${projectname}.subsequence_data.hq.${minreads}.hdf5";
of_all="${OUTDIR}/${projectname}.subsequence_data.all.${minreads}.hdf5"
bin_extract="$DIRICORE_DIR/diricore/bin/extract_subsequences.py"
bin_plot="$DIRICORE_DIR/diricore/bin/plot_subsequence_shifts.py"


###
mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}/subsequence_shift_plots"

if [[ $plots_only == 0 ]]; then
    echo "Extracting subsequences (hqmapped)"
    ls -1 ${INDIR}/*.hqmapped.bam | sort -V | while read bamfn; do
        b=$(basename "$bamfn");
        b=${b%%_toGenome.hqmapped.bam};
       $DIRICORE_DIR/diricore/bin/extract_subsequences.py \
           -v \
           run \
           -o $of \
           -f 0 \
           ${INDEXDATAFN} \
           ${b},${bamfn}
    done
    echo "Done"

#    echo "Extracting subsequences (all)"
#    ls -1 ${INDIR}/*_toGenome.bam | sort -V | while read bamfn; do
#        b=$(basename "$bamfn");
#        b=${b%%_toGenome.bam};
#       $DIRICORE_DIR/diricore/bin/extract_subsequences.py \
#           -v \
#           run \
#           -o $of_all \
#           -f 0 \
#           ${INDEXDATAFN} \
#           ${b},${bamfn}
#    done
#    echo "Done"
fi

# create subsequence shift plots
if [[ -f $of ]]; then
  echo "Plotting subsequences (hq)"
  $bin_plot \
    -o ${PLOTDIR}/subsequence_shift_plots/${projectname}.m${minreads}.hq. \
    -m $minreads \
    --sample-names ${SAMPLENAME_FILE} \
    --contrasts ${CONTRAST_FILE} \
    --y-limits $y_limits \
    ${of}
  echo "Done. Generated file: ${PLOTDIR}/subsequence_shift_plots/${projectname}.m${minreads}.hq."
else
  echo "File $of is missing! Skipping"
fi

#if [[ -f $of ]]; then
#  echo "Plotting subsequences (all)"
#  $bin_plot \
#    -o ${PLOTDIR}/subsequence_shift_plots/${projectname}.m${minreads}.all. \
#    -m $minreads \
#    --sample-names ${SAMPLENAME_FILE} \
#    --contrasts ${CONTRAST_FILE} \
#    --y-limits $y_limits \
#    ${of_all}
#  echo "Done. Generated file: ${PLOTDIR}/subsequence_shift_plots/${projectname}.m${minreads}.all."
#else
#  echo "File $of_all is missing! Skipping"
#fi
###
