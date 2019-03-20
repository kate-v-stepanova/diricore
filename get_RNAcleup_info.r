#!/usr/bin/env Rscript

require(data.table)

INDIR="./data/output/clean"
OUTDIR="./data/output/alignment_stats"

list_of_files <- list.files(INDIR, full.names=T)
names(list_of_files) <- basename(list_of_files)

# Left after rrna and used for alignment
sel <- list_of_files[grepl('*trna.err', list_of_files)]
dt <- do.call('rbind', lapply(sel, function(x){
      x <- readLines(x)
      data.table('rrnaleft' = as.numeric(gsub('([0-9]+).*',"\\1",x[1])),
                  'trnaleft' = as.numeric(gsub('([0-9]+).*',"\\1",x[3])))
    }))
dt$file = gsub('\\.trna\\.err',"",basename(sel))
write.table(dt, paste(OUTDIR, '/RNA_cleup_stats.txt', sep=""), quote=F, row.names=F, sep='\t')
