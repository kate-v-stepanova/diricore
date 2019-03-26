# Analysis of 13728 dataset - HiSeq

The diricore version used for the analysis:

```
commit 6b1d20fed0c93b053ba6fc3c221a15a2cb98e905 (HEAD -> master, origin/master, origin/HEAD)
Author: b250-admin <b250-admin@dkfz-vpn94.inet.dkfz-heidelberg.de>
Date:   Tue Mar 26 11:31:26 2019 +0100

    Second working version. Tested on Hiseq dataset
```


## Download data

Go to: `https://ilseform.dkfz.de/ilse-wui/secure/mySubmissions.jsf` -> Download -> Outside DKFZ
Copy & paste link on the cluster

````bash
mkdir /icgc/dkfzlsdf/analysis/OE0532/13728
sftp -r ad+e984a@ftp4midterm.dkfz.de:013728/data/* /icgc/dkfzlsdf/analysis/OE0532/13728
````

## Create the directory structure

```
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/input
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/input/fastq
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/input/metadata
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/input/merged

mkdir -p /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/output/demultiplexed
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/output/umi_extract/logs

```

## Generate bc_file

```
cd ~/diricore
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/input/metadata/
Rscript utils/make_bc_file.r 13728
```

The output file will be stored in `data/input/metadata/bc_file.txt`


## Redemultiplex files
The facility attempted to demultiplex using the internal barcodes, but it failed. So, we need to merge all the files into one, remove the adapter and demultiplex them again. For all these steps execute a command:

```
bsub -q long ./trim_adapter_and_demultiplex.sh /icgc/dkfzlsdf/analysis/OE0532/13728/050106_D00636_0075_BCCDYVANXX
```

This script includes 3 steps:
1. Merging files from the flowcell_path into one: `./data/input/merged/050106_D00636_0075_BCCDYVANXX.fastq.gz`
2. Removing adapter `AGATCGGAAGAGCACACGTCTGAA` and splitting the merged file according to the given barcodes. The barcodes must be defined in `./data/input/metadata/bc_file.txt`. And the input file must be in  `./data/input/merged/050106_D00636_0075_BCCDYVANXX.fastq.gz`. The output will be: `./data/output/demultipexed/*.fastq`. This will also include `unmached.fastq`, which maybe a good idea to remove
3. Removing UMIs
The given pattern is `NNNNN`, which is a random sequence of the length 5.
After the removing UMIs, the output files will be written into `analysis/output/umi_extract_path`, but to save the space, it will be moved (not copied) to `analysis/input/fastq`
Maybe need to optimize this step and write directly to `analysis/input/fastq`

Output:

```
ls -lh analysis/output/demultiplexed/
total 13G
-rw-r--r--. 1 e984a B250 756M Mar 22 14:16 dem_CDK1_inh.fastq
-rw-r--r--. 1 e984a B250 1.6G Mar 22 14:16 dem_C.fastq
-rw-r--r--. 1 e984a B250 1.8G Mar 22 14:16 dem_Control.fastq
-rw-r--r--. 1 e984a B250 2.0G Mar 22 14:16 dem_Ha.fastq
-rw-r--r--. 1 e984a B250 1.1G Mar 22 14:16 dem_Ha_Roca.fastq
-rw-r--r--. 1 e984a B250 2.2G Mar 22 14:16 dem_Roca.fastq
-rw-r--r--. 1 e984a B250 711M Mar 22 14:16 dem_unmatched.fastq
```
dem_Control.fastq can be ignored for contrasts

## rRNA clean up - run_umi_preprocessing

````bash
./run_umi_preprocessing.sh human 13728
````

---
## Alignment
The alignment can be run normally using the diricore script
````bash
bsub -q long ./run_alignment.sh human 13728
````

## Deduplicate by UMI

```
bsub ./deduplicate_umi.sh 13728
```

NEXT TIME: do the RPF densities first, then subsequence analysis, then counts and so on.

## Analyse counts

```
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/output/counts/
bsub ./get_counts.sh 13728
module load R
bsub ./analyze_counts.r 13728
```

Expected output from `get_counts.sh`:  

```
ls /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/output/counts/
htseq_counts_dedup.txt  htseq_counts.txt
``` 

Expected output from `analyze_counts.r`: 

```
ls /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/output/counts/
run_UMIs_htseq_cpms_prot_cod.tsv  run_UMIs_htseq_dedup_cpms_prot_cod.tsv
```

## Analyse fragments

```
bsub ./fragment_size_1.sh 13728
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/13728/analysis/output/figures/fragment_size
bsub ./fragment_size_2.r 13728

```

## Analyse alignments

DOESN'T work. Need to see why the `bc_stats.txt` files from test dataset is different from the current dataset. Skipping for now. 

```
bsub ./get_alignment_stats.sh 13728
bsub ./get_RNAcleup_info.r 13728
bsub ./plot_alignment_stats.r 13728
```

## Generate tracks

```
bsub ./generate_tracks.sh 13728
bsub ./generate_trackDB.sh 13728
```

## Find the top rRNA fragments for substraction

The first one takes some time.

```
bsub -q long ./find_top_RNA_fragments_1.sh 13728
module load R
bsub ./find_top_RNA_fragments_2.r 13728
```

## RPF Density
Requires `analysis/input/metadata/rpf_density_contrasts.tsv`:

```
Ha      C   #ff99bb
Roca    C   #80e5ff
Ha_Roca C   #99ffdd
```

and `analysis/input/metadata/rpf_density_samplenames.tsv` :

```
C       Ctrl
Ha      Harr
Roca    Harr
Ha_Roca Harr
```

Run the analysis:

```
bsub -q medium ./run_umi_rpf_density_analysis.sh 13728 human 25
```


## Subsequence analysis

```
bsub ./run_umi_subsequence_analysis.sh 13728 human 25
```

## RPF transcript distribution

Requires `/icgc/dkfzlsdf/analysis/OE0532/13728/analysis/input/metadata/rpf_transcript_distribution_sampleinfo.tsv`:

```
Ha      C       #ff99bb
Roca    C       #80e5ff
Ha_Roca C       #99ffdd
```

Running script:

```
./plot_rpf_transcript_distribution.sh 13728 25
```
