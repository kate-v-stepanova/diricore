#!/bin/bash

set -e
set -u

dataset_id=$1
genome=$2

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR="$BASE_DIR/$dataset_id"
#INDIR="$PROJECT_DIR/analysis/input/fastq"
INDIR="$PROJECT_DIR/analysis/output/clean"
STAR_GENOME_DIR="/icgc/dkfzlsdf/analysis/OE0532/static/$genome"
OUTDIR=$PROJECT_DIR/analysis/output/alignments # after preprocessing
script_dir="/icgc/dkfzlsdf/analysis/OE0532/tmp/faster_alignment/$dataset_id"

mkdir -p $script_dir
mkdir -p $OUTDIR/toGenome
mkdir -p $OUTDIR/toTranscriptome
mkdir -p $OUTDIR/reads_per_gene
mkdir -p $OUTDIR/logs

for f in $(ls ${INDIR}/*.fastq.gz); do
    b=$(basename ${f});
    b=${b%%.*};
    mkdir -p $script_dir
    outfile="$script_dir/$b.sh"
    echo "#!/bin/bash" > $outfile
    echo "set -e" >> $outfile # exit on error, don't continue
    # align
    echo "module load STAR" >> $outfile
    echo "STAR --genomeDir $STAR_GENOME_DIR --runThreadN 100 --readFilesIn ${f} --outFileNamePrefix $OUTDIR/${b}_ --outSAMtype BAM Unsorted --readFilesCommand zcat  --quantMode TranscriptomeSAM GeneCounts" >> $outfile
	# sort transcriptome
	echo "module load samtools/1.6" >> $outfile
	sorted="$OUTDIR/$b.sorted.trans.bam"
    # samtools sort: option -T is important!! this is where the tmp will be written. Otherwise there will be problems with:
	#  1) running 2 datasets with equal sample names, the bam files will be mixed up!
    #  2) hard drive quota will be exceeded.
	echo "samtools sort $OUTDIR/${b}_Aligned.toTranscriptome.out.bam -o $sorted -T $OUTDIR/$b" >> $outfile
	echo "mv $sorted $OUTDIR/${b}_Aligned.toTranscriptome.out.bam" >> $outfile
	echo "mv $OUTDIR/${b}_Aligned.toTranscriptome.out.bam ${OUTDIR}/toTranscriptome/${b}_toTranscriptome.bam" >> $outfile
	# sort genome
	sorted="$OUTDIR/$b.sorted.bam"
	echo "samtools sort $OUTDIR/${b}_Aligned.out.bam -o $sorted -T $OUTDIR/$b" >> $outfile
	echo "mv $sorted ${b}_Aligned.out.bam" >> $outfile
	echo "mv ${OUTDIR}/${b}_Aligned.out.bam ${OUTDIR}/toGenome/${b}_toGenome.bam" >> $outfile
	# mv log files
	echo "mv $OUTDIR/${b}*ReadsPerGene.out.tab $OUTDIR/reads_per_gene" >> $outfile
	echo "mv $OUTDIR/${b}*.out $OUTDIR/logs/" >> $outfile
	echo "mv $OUTDIR/${b}*.tab $OUTDIR/logs/" >> $outfile
	## extract hq reads
	gen_bam="$OUTDIR/toGenome/${b}_toGenome.bam"
	gen_hq_bam="$OUTDIR/toGenome/${b}_toGenome.hqmapped.bam"
	trans_bam="${OUTDIR}/toTranscriptome/${b}_toTranscriptome.bam"
	trans_hq_bam="${OUTDIR}/toTranscriptome/${b}_toTranscriptome.hqmapped.bam"
	echo "cat <(samtools view -H $gen_bam)  <(cat $gen_bam | samtools view -q10 -F260 -) | samtools view -bS - > $gen_hq_bam " >> $outfile
	echo "cat <(samtools view -H $trans_bam)  <(cat $trans_bam | samtools view -q10 -F260 -) | samtools view -bS - > $trans_hq_bam " >> $outfile
	
	chmod +x $outfile
	echo "bsub -q long $outfile"
done

