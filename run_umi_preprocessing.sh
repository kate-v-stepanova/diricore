#!/bin/bash


set -e;
set -u;


export BOWTIE_BIN="./programs/bowtie2-2.0.6/bowtie2";
export INDIR="./data/input/fastq";
export OUTDIR="./data/output/clean";
echo "Settings:";
export species=$1; echo -e "\tspecies:$species";

export RRNA_REF="./staticdata/${species}/rRNAs";
export TRNA_REF="./staticdata/${species}/tRNAs";


###
mkdir -p $OUTDIR || true;

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
    cat "${fn}" \
    | gzip -dc \
    | ${BOWTIE_BIN} --seed 42 -p 1 --local --un-gz "${tmpfile}" "${RRNA_REF}" - \
    > /dev/null 2> "${rrna_err}";

    cat "${tmpfile}" \
    | gzip -dc \
    | ${BOWTIE_BIN} --seed 42 --local --un-gz "${of}" "${TRNA_REF}" - \
    > /dev/null 2> "${trna_err}";

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
