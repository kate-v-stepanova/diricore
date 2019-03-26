#!/usr/bin/env Rscript
require(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
   stop(paste("Usage example:", args[0], "14548"))
} else if (length(args)==1) {
   dataset_id=args[1]
}

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR=paste(BASE_DIR, dataset_id, sep="/")
INDIR=paste(PROJECT_DIR, "analysis/output/clean", sep="/")
OUTDIR=paste(PROJECT_DIR, "analysis/output/alignment_stats", sep="/")
outfile=paste(OUTDIR, "RNA_clenup_stats.txt", sep="/")

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
print( paste("Writing", outfile))
write.table(dt, outfile, quote=F, row.names=F, sep='\t')
print( "Done")
