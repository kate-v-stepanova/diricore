#!/bin/bash


set -e;
set -u;

dataset_id=$1
species=$2

BASE_PATH="/icgc/dkfzlsdf/analysis/OE0532"
DIRICORE_PATH="/home/e984a/diricore"

BOWTIE_BIN="$DIRICORE_PATH/programs/bowtie2-2.0.6/bowtie2";

PROJECT_PATH="$BASE_PATH/$dataset_id"
INDIR="$PROJECT_PATH/analysis/input/fastq";
OUTDIR="$PROJECT_PATH/analysis/output/clean";

RRNA_REF="$BASE_PATH/static/${species}/rRNAs";
TRNA_REF="$BASE_PATH/static/${species}/tRNAs";

TMP_DIR="$BASE_PATH/tmp"

project_id=$dataset_id
###
mkdir -p $OUTDIR

run_prep() {
    fn=$1
    bn=$(basename "$fn");
    b=${bn%%.*};

    of="${OUTDIR}/${b}.fastq.gz";
    rrna_err="${OUTDIR}/${b}.rrna.err";
    trna_err="${OUTDIR}/${b}.trna.err";

    tmpfile="$TMP_DIR/${project_id}_${b}.rrna_cleaned.tmp.fastq.gz";
    rm -f $tmpfile;
    $(cp ${INDIR}/${b}.fastq.gz $tmpfile);
    trap "{ rm -f ${tmpfile}; }" EXIT;

    echo "Starting preprocessing of file: ${bn}";
    # added -p 50 (number of processors)
    cat "${fn}" | gzip -dc | ${BOWTIE_BIN} --seed 42 -p 50 --local --un-gz "${tmpfile}" "${RRNA_REF}" -  > /dev/null 2> "${rrna_err}";

    cat "${tmpfile}" | gzip -dc  | ${BOWTIE_BIN} --seed 42 -p 50 --local --un-gz "${of}" "${TRNA_REF}" -  > /dev/null 2> "${trna_err}";

    rm ${tmpfile};
}
export -f run_prep
for f in `ls ${INDIR}/*.fastq.gz`; do
  fn=$(basename $f)
  of=${OUTDIR}/$fn
  if [[ ! -f $of ]]; then
    echo "Preprocessing ${f}";
    run_prep ${f};
    echo "Done ${f}";
  else
    echo "File exists! Skipping: $of"
  fi
done;

echo "Done with preprocessing";
###
