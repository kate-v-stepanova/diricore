# Analysis of 14548 dataset - DMS
Version of diricore from `~/diricore_dev`

## Download the data

```
./utils/get_data.sh 14548
```
This script will copy data from the `midterm` to our group directory on the cluster `/icgc/dkfzlsdf/analysis/OE0532`

## Generate bc_file

First, make sure that the directories exist. 

```
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/14548/analysis/input/metadata
```

```
./utils/make_bc_file.r 14548
```

## Trim and demultiplex

```
bsub -q long ./trim_adapter_and_demultiplex.sh /icgc/dkfzlsdf/analysis/OE0532/14548/190308_M01688_0156_000000000-C9C58
```
--
Results are the following:

```
ls -lh analysis/output/demultiplexed/
total 3.0G
-rw-r--r--. 1 e984a B250 1014M Mar 22 13:31 dem_Control_DMS.fastq
-rw-r--r--. 1 e984a B250  1.2G Mar 22 13:31 dem_Control.fastq
-rw-r--r--. 1 e984a B250  101M Mar 22 13:31 dem_Noco_DMS.fastq
-rw-r--r--. 1 e984a B250  134M Mar 22 13:31 dem_Roca_DMS.fastq
-rw-r--r--. 1 e984a B250   29M Mar 22 13:31 dem_unmatched.fastq
```

```
ls -lh analysis/input/fastq/
total 463M
-rw-r--r--. 1 e984a B250 139M Mar 22 13:35 dem_Control_DMS_umi_extracted.fastq.gz
-rw-r--r--. 1 e984a B250 190M Mar 22 13:41 dem_Control_umi_extracted.fastq.gz
-rw-r--r--. 1 e984a B250  15M Mar 22 13:42 dem_Noco_DMS_umi_extracted.fastq.gz
-rw-r--r--. 1 e984a B250  20M Mar 22 13:43 dem_Roca_DMS_umi_extracted.fastq.gz
-rw-r--r--. 1 e984a B250 5.1M Mar 22 13:44 dem_unmatched_umi_extracted.fastq.gz
```

## rRNA and tRNA clean up - UMI preprocessing

```
bsub -q long ./run_umi_preprocessing.sh human 14548

```
Input from `/icgc/dkfzlsdf/analysis/OE0532/14548/analysis/input/fastq`
Output in `/icgc/dkfzlsdf/analysis/OE0532/14548/analysis/output/clean`


## Alignment

```
bsub -q long ./run_alignment.sh human 14548
```


Done. Proceed with DMS pipeline

