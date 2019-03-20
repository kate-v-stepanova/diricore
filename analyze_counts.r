#!/usr/bin/env Rscript

counts_file <- "./data/output/counts/htseq_counts.txt"
counts_dedup_file <- "./data/output/counts/htseq_counts_dedup.txt"

require(data.table)
dt <- fread(counts_file)
dt_dedup <- fread(counts_dedup_file)

# Fix colnames
#:ncols(dt)-1
last_col=ncol(dt)
last_col_dedup=ncol(dt_dedup)
print(last_col)
print(last_col_dedup)
colnames(dt)[4:last_col] <- gsub('.*run_umi_(.*)_umied.*',"\\1",colnames(dt)[4:last_col])
colnames(dt_dedup)[4:last_col_dedup] <- gsub('.*run_umi(.*)_umied.*', "\\1",colnames(dt_dedup)[4:last_col_dedup])
# Filter stats
dt <- dt[!grepl('^__', Gene_id)]
dt_dedup <- dt_dedup[!grepl('^__', Gene_id)]

# Filter low counts
dt <- dt[rowSums(dt[,4:last_col]) > 20]
 dt_dedup <- dt_dedup[rowSums(dt_dedup[,4:last_col_dedup]) > 20]

# Filter coding genes
dt <- dt[Gene_type == 'protein_coding']
dt_dedup <- dt_dedup[Gene_type == 'protein_coding']

# Normalize
cols <- colnames(dt)[4:last_col]
cols_dedup = colnames(dt_dedup)[4:last_col_dedup]
dt[,(cols):=lapply(.SD, function(x) x*100000/sum(x)),.SDcols=cols]
dt_dedup[,(cols_dedup):=lapply(.SD, function(x) x*100000/sum(x)), .SDcols=cols_dedup]
write.table(dt, './data/output/counts/run_UMIs_htseq_cpms_prot_cod.tsv', sep='\t',quote=F, row.names=F)
write.table(dt_dedup, './data/output/counts/run_UMIs_htseq_dedup_cpms_prot_cod.tsv', sep='\t', quote=F, row.names=F)

