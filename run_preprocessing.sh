#!/bin/bash

set -e;
set -u;

dataset_id=$1
species=$2

BASE_PATH="/icgc/dkfzlsdf/analysis/OE0532"
DIRICORE_PATH="/icgc/dkfzlsdf/analysis/OE0532/software/diricore"

BOWTIE_BIN="$DIRICORE_PATH/programs/bowtie2-2.0.6/bowtie2";

PROJECT_PATH="$BASE_PATH/$dataset_id"
INDIR="$PROJECT_PATH/analysis/input/fastq";
OUTDIR="$PROJECT_PATH/analysis/output/clean";

RRNA_REF="$BASE_PATH/static/${species}/rRNAs";
TRNA_REF="$BASE_PATH/static/${species}/tRNAs";

project_id=$dataset_id
p_id=${project_id/"/"/"_"}

TMP_DIR="$BASE_PATH/tmp/$p_id"
script_dir="$TMP_DIR/faster_preprocessing/$dataset_id"


###
mkdir -p $OUTDIR
mkdir -p $script_dir
mkdir -p $TMP_DIR

run_prep() {
    fn=$1
    bn=$(basename "$fn");
    b=${bn%%.*};

    of="${OUTDIR}/${b}.fastq.gz";
    rrna_err="${OUTDIR}/${b}.rrna.err";
    trna_err="${OUTDIR}/${b}.trna.err";

    tmpfile="$TMP_DIR/${p_id}_${b}.rrna_cleaned.tmp.fastq.gz";
    
    script_file=$script_dir/${b}.sh

    rm -f $tmpfile;
    echo "#!/bin/bash" > $script_file
    echo "set -e;" >> $script_file
    echo "cp ${INDIR}/${b}.fastq.gz $tmpfile" >> $script_file
    echo "cat ${fn} | gzip -dc | ${BOWTIE_BIN} --seed 42 -p 10 --local --un-gz ${tmpfile} ${RRNA_REF} -  > /dev/null 2> ${rrna_err}" >> $script_file
    echo "cat ${tmpfile} | gzip -dc  | ${BOWTIE_BIN} --seed 42 -p 10 --local --un-gz ${of} ${TRNA_REF} -  > /dev/null 2> ${trna_err}" >> $script_file
    echo "rm ${tmpfile}" >> $script_file
    chmod +x $script_file
    echo "bsub -q long -R \"rusage[mem=10G]\" $script_file"
}
export -f run_prep
for f in `ls ${INDIR}/*.fastq.gz`; do
  run_prep ${f};
done;

