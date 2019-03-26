#  Analysis of dataset 14464 (miseq_data)

## Fetching the data

```
./utils/get_data.sh 14464
```

## Generating bc_file

```
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/14548/analysis/input/metadata
./utils/make_bc_file.r 14548
```


## Trim and demultiplex

```
bsub ./trim_adapter_and_demultiplex.sh /icgc/dkfzlsdf/analysis/OE0532/14464/190306_M01688_0155_000000000-D5J58
```

## Cleanup rRNAs

```
bsub ./run_umi_preprocessing.sh human 14464
```

## Alignment

```
bsub ./run_alignment.sh human 14464
```

## Deduplicate

This script usually runs very fast, no need to submit it as a job, but anyway:

```
bsub ./deduplicate_umi.sh 14464
```

## Get counts

Required `staticdata/human/gencode.v29.basic.annotation.gtf` file. Make sure that the file is present and then run script:

```
bsub -q medium ./get_counts.sh 14464
```

## Analyze counts

First, make sure that the output dir exists:

```
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/14464/analysis/output/counts/
```

Then, run the analysis:

```
bsub -q medium ./analyze_counts.r 14464
```
It will finish in a few seconds, so maybe no need to `bsub`. 


## Get and plot fragments 
```
bsub ./fragment_size_1.sh 14464
mkdir -p /icgc/dkfzlsdf/analysis/OE0532/14464/analysis/output/figures/fragment_size
bsub ./fragment_size_2.r 14464
```

## Alignment stats
```
./get_alignment_stats.sh 14464
./get_RNAcleup_info.r 14464
./get_bc_split_stats.sh 14464

```
