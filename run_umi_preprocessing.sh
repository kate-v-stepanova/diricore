#!/bin/bash


set -e;
set -u;

species=$1
dataset_id=$2

BASE_PATH="/icgc/dkfzlsdf/analysis/OE0532"
DIRICORE_PATH="/home/e984a/diricore"

BOWTIE_BIN="$DIRICORE_PATH/programs/bowtie2-2.0.6/bowtie2";

PROJECT_PATH="$BASE_PATH/$dataset_id"
INDIR="$PROJECT_PATH/analysis/input/fastq";
OUTDIR="$PROJECT_PATH/analysis/output/clean";

RRNA_REF="$DIRICORE_PATH/staticdata/${species}/rRNAs";
TRNA_REF="$DIRICORE_PATH/staticdata/${species}/tRNAs";


###
mkdir -p $OUTDIR

run_prep() {
    fn=$1
    bn=$(basename "$fn");
    b=${bn%%.*};

    of="${OUTDIR}/${b}.fastq.gz";
    rrna_err="${OUTDIR}/${b}.rrna.err";
    trna_err="${OUTDIR}/${b}.trna.err";

    tmpfile="/tmp/${b}.rrna_cleaned.tmp.fastq.gz";
    rm -f $tmpfile;
    $(cp ${INDIR}/${b}.fastq.gz $tmpfile);
    trap "{ rm -f ${tmpfile}; }" EXIT;

    echo "Starting preprocessing of file: ${bn}";
    cat "${fn}" | gzip -dc | ${BOWTIE_BIN} --seed 42 -p 1 --local --un-gz "${tmpfile}" "${RRNA_REF}" -  > /dev/null 2> "${rrna_err}";

    cat "${tmpfile}" | gzip -dc  | ${BOWTIE_BIN} --seed 42 --local --un-gz "${of}" "${TRNA_REF}" -  > /dev/null 2> "${trna_err}";

    rm ${tmpfile};
}
export -f run_prep
for f in `ls ${INDIR}/*.fastq.gz`; do
  echo "Preprocessing ${f}";
  run_prep ${f};
  echo "Done ${f}";
done;

echo "Done with preprocessing";
###
