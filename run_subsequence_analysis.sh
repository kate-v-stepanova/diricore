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

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"

OUTDIR="$PROJECT_DIR/analysis/output/subsequence_data";
INDIR="$PROJECT_DIR/analysis/output/tophat_out";
PLOTDIR="$PROJECT_DIR/analysis/output/figures";
SAMPLENAME_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv";
CONTRAST_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";

DIRICORE_DIR="/home/e984a/diricore"
INDEXDATAFN="$DIRICORE_DIR/staticdata/${species}/subseq_index_data.pkl.gz";
of="${OUTDIR}/${projectname}.subsequence_data.hdf5";
bin_extract="$DIRICORE_DIR/diricore/bin/extract_subsequences.py"
bin_plot="$DIRICORE_DIR/diricore/bin/plot_subsequence_shifts.py"


###
mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}/subsequence_shift_plots"

if [[ $plots_only == 0 ]]; then
    echo "Exctracting subsequences"
    ls -1 ${INDIR}/*.bam | sort -V | while read bamfn; do
        b=$(basename "$bamfn");
        b=${b%%.*};
       $DIRICORE_DIR/diricore/bin/extract_subsequences.py \
           -v \
           run \
           -o $of \
           -f 0 \
           ${INDEXDATAFN} \
           ${b},${bamfn}
    done

#  ls -1 ${INDIR}/*.bam | sort -V | while read bamfn; do
#    b=$(basename $(dirname "$bamfn"));
#      b=$(basename $bamfn)
#      b=${b%%.*}
#      ${bin_extract} -v run -o ${of} -f 0 ${INDEXDATAFN} ${b},${bamfn}
#  done
  echo "Done"
fi

# create subsequence shift plots
echo "Plotting subsequences"
$bin_plot \
    -o ${PLOTDIR}/subsequence_shift_plots/${projectname}.m${minreads}. \
    -m $minreads \
    --sample-names ${SAMPLENAME_FILE} \
    --contrasts ${CONTRAST_FILE} \
    ${of}
echo "Done. Generated file: ${of}"
###
