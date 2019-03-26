#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
   stop(paste("Usage example:", args[0], "14548"))
} else if (length(args)==1) {
   dataset_id=args[1]
}

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR=paste(BASE_DIR, dataset_id, sep="/")
counts_file=paste(PROJECT_DIR, "analysis/output/counts/htseq_counts.txt", sep="/")
counts_dedup_file=paste(PROJECT_DIR, "analysis/output/counts/htseq_counts_dedup.txt", sep="/")

# Reading input files
print("Reading input files")
require(data.table)
dt <- fread(counts_file)
dt_dedup <- fread(counts_dedup_file)


# Fix colnames
print("Fixing colnames")
last_col=ncol(dt)
last_col_dedup=ncol(dt_dedup)
colnames(dt)[4:last_col] <- gsub('.*run_umi_(.*)_umied.*',"\\1",colnames(dt)[4:last_col])
colnames(dt_dedup)[4:last_col_dedup] <- gsub('.*run_umi(.*)_umied.*', "\\1",colnames(dt_dedup)[4:last_col_dedup])

# Filter stats
print("Filtering stats")
dt <- dt[!grepl('^__', Gene_id)]
dt_dedup <- dt_dedup[!grepl('^__', Gene_id)]

# Filter low counts
print("Filtering low counts")
dt <- dt[rowSums(dt[,4:last_col]) > 20]
 dt_dedup <- dt_dedup[rowSums(dt_dedup[,4:last_col_dedup]) > 20]

# Filter coding genes
print("Filtering coding genes")
dt <- dt[Gene_type == 'protein_coding']
dt_dedup <- dt_dedup[Gene_type == 'protein_coding']

# Normalize
print("Normalization")
cols <- colnames(dt)[4:last_col]
cols_dedup = colnames(dt_dedup)[4:last_col_dedup]
dt[,(cols):=lapply(.SD, function(x) x*100000/sum(x)),.SDcols=cols]
dt_dedup[,(cols_dedup):=lapply(.SD, function(x) x*100000/sum(x)), .SDcols=cols_dedup]

# Write results
print("Writing results")
write.table(dt, paste(PROJECT_DIR, 'analysis/output/counts/run_UMIs_htseq_cpms_prot_cod.tsv', sep="/"), sep='\t',quote=F, row.names=F)
write.table(dt_dedup, paste(PROJECT_DIR, 'analysis/output/counts/run_UMIs_htseq_dedup_cpms_prot_cod.tsv', sep="/"), sep='\t', quote=F, row.names=F)
print("Done")
