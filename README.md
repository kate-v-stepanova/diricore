#  Diricore pipeline installation
This is a diricore version from E.Stepanova. For the other versions there is other documentation (in docs directory).

```
git clone https://github.com/kate-v-stepanova/diricore.git
cd diricore
conda create -n diricore python=2.7
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
pip install -r requirements.txt
```

In `diricore/programs` there are installed packages required by diricore.  In ideal case everything will work smoothly just like that. If not, install the packages manually:

1. Samtools: version 0.1.19

`https://sourceforge.net/projects/samtools/files/samtools/0.1.19/`

2. cutadapt 

```
pip install cutadapt
```

3. Bowtie: version 2.0.6

```
https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.0.6/
unzip bowtie2-2.0.6*.zip
cp -R bowtie2-2.0.6 diricore_pipeline/programs
```

4. One of the scripts requires [seqtk](https://github.com/lh3/seqtk)

```
git clone https://github.com/lh3/seqtk
cd seqtk
make
```
# Steps of the diricore
## Preprocessing
This includes demultiplexing, trimming adapter, and extracting UMIs

The command is the following: `trim_adapter_and_demultiplex.sh <project_id>`

This script requires the following input files:
* fastq file located at `<BASE_DIR>/<project_id>/analysis/input/merged/<filename>.fastq.gz`
* Barcodes information stored at `<BASE_DIR>/<project_id>/analysis/input/metadata/bc_file.txt`
The file should contain a sample name and a barcode sequence divided by **tab**. For example: 
```
Control ATCGT
CDK1_inh        AGCTA
Roca    CTAGA
```

The script includes the following steps:
1. Trimming an adapter, which is the same for all datasets, so it is hardcoded in the script. This is done with `cutadapt`

The output file be in `<BASE_DIR>/<project_id>/analysis/input/merged/<filename>.trimmed.fastq`
Additionally, we get the stats file `<BASE_DIR>/<project_id>/analysis/output/cutadapt_trimming_stats.txt`, which contains the number of reads of the input file, number of reads with adapter, number of reads too short, and number of reads in the output file.

What is important to know, it also removes 3 extra bases from the beginning of each read (5' end). 

2. Demultiplexing.

This is done with [`fastx_toolkit/fastx_barcode_splitter.pl`](http://hannonlab.cshl.edu/fastx_toolkit/)

Input file: `<BASE_DIR>/<project_id>/analysis/input/merged/<filename>.trimmed.fastq` 

The output files will be in `<BASE_DIR>/<project_id>/analysis/output/demultiplexed/dem_*.fastq`, where `dem_` is a prefix indicating that the sample has been demultiplexed and `*` is a samplename (from `bc_file.txt`). 

Additionally we get a stats file `<BASE_DIR>/<project_id>/analysis/output/bc_split_stats.txt`, which contains number of reads per sample, number of discarded reads, and total number of reads.

3. Extracting UMIs

This done with `umi_tools`. 

Input files: `<BASE_DIR>/<project_id>/analysis/output/demultiplexed/dem_*.fastq`

Output files: `<BASE_DIR>/<project_id>/analysis/output/umi_extracted/dem_*_umi_extracted.fastq`

Here the files can be finally gzipped: `gzip <BASE_DIR>/<project_id>/analysis/output/umi_extracted/*.fastq`

The other ones too, for the sake of storage space, but not for diricore.

4. Linking final fastq files to `input` directory. 

It is not (yet) in the script, has to be done manually:
```
for f in $(ls <BASE_DIR>/<project_id>/analysis/output/umi_extracted/*.fastq.gz); do samplename=$(basename $f); samplename=${samplename%"_umi_extracted.fastq.gz"}; samplename=${samplename#"dem_"}; echo "ln -s $f <BASE_DIR>/<project_id>/analysis/input/fastq/${samplename}.fastq.gz"; done
```

Input: `<BASE_DIR>/<project_id>/analysis/output/umi_extracted/*.fastq.gz`

Output: `<BASE_DIR>/<project_id>/analysis/input/fastq/*.fastq.gz`

## Removing rRNAs and tRNAs

Here we remove the reads aligning to rRNA and tRNA genes. 

Input: `<BASE_DIR>/<project_id>/analysis/input/fastq/*.fastq.gz`

Output: `<BASE_DIR>/<project_id>/analysis/output/clean/*.fastq.gz`

## Alignment

### Tophat

Previously alignment was done with `tophat`. In this case run the following command: `run_alignment.sh <project_id> <genome>`
This script requires a bunch of bowtie index files which were together with diricore, the origins of these files are unknown, the only thing I was able to figure out, they belong to `hg19`. 

For the current version of diricore these files have to be located at `<BASE_DIR>/static/<genome>`, where genome is `hg19`. 
This script will create a subdirectory for each of the samples, and will have 2 bam files per sample: `<sample_dir>/accepted_hits.bam`, `<sample_dir>/accepted_hits.hqmapped.bam`, where hqmapped files are the ones that have a quality scores above a certain threshold. 

Input: `<BASE_DIR>/<project_id>/analysis/output/clean/*.fastq.gz`

Output: `<BASE_DIR>/<project_id>/analysis/output/tophat_out/<sample_dir>/accepted_hits.bam` and `<BASE_DIR>/<project_id>/analysis/output/tophat_out/<sample_dir>/accepted_hits.hqmapped.bam`

### STAR

This method is the one I personally prefer, so all the next scripts were recently adapted for the output of STAR alignment.
Before running the alignment, the index has to be created:

#### Preparing the index

```
hg19
https://www.gencodegenes.org/human/release_19.html

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

STAR --runThreadN 100 --runMode genomeGenerate --genomeDir <BASE_DIR>/static/hg19 --genomeFastaFiles /icgc/dkfzlsdf/analysis/OE0532/static/hg19/GRCh37.p13.genome.fa --sjdbGTFfile /icgc/dkfzlsdf/analysis/OE0532/static/hg19/gencode.v19.annotation.gtf --sjdbOverhang 100
```
The same for hg38. For the consistency (so that both versions have the same filenames), I created softlinks: 
```
ln -s gencode.v19.annotation.gtf gencode.annotation.gtf
ln -s GRCh37.p13.genome.fa genome.fa
```
I don't rename/remove the original files, so that I (or someone else) will be able to figure out where these files are coming from. 

### Running the alignment

Should be run with the following command: `align_to_transcriptome.sh <project_id> <genome>`. 

The script has a possibility to align both to the genome and to the transcriptome, for now I switch between those options by just commenting out a line with `--quantMode TranscriptomeSAM GeneCounts`. Let's say at some point I added a parameter to the script, e.g. `to_genome/to_transcriptome`.

To align to genome only: `align_to_transcriptome.sh <project_id> <genome> to_genome`

To align to both genome and transcriptome: `align_to_transcriptome.sh <project_id> <genome> to_transcriptome`

Input: `<BASE_DIR>/<project_id>/analysis/output/clean/*.fastq.gz`

Output (to genome): `<BASE_DIR>/<project_id>/analysis/output/alignments/toGenome/<samplename>.bam` and `<BASE_DIR>/<project_id>/analysis/output/alignments/toGenome/<samplename>.hqmapped.bam`

Output (to transcriptome): `<BASE_DIR>/<project_id>/analysis/output/alignments/toTranscriptome/<samplename>.bam` and `<BASE_DIR>/<project_id>/analysis/output/alignments/toTranscriptome/<samplename>.ReadsPerGene.out.tab`

## Deduplication

Final step of the pre-processing.

Deduplication should be run with the command: `bsub -q long ./deduplicate_umi.sh <project_id>`

Input: `<BASE_DIR>/<project_id>/analysis/output/alignments/toGenome/<samplename>.bam`

Output: `<BASE_DIR>/<project_id>/analysis/output/alignments/toGenome/<samplename>_dedup.bam`

