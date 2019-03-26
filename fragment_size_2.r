#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop(paste("Usage example:", args[0], "14548"))
} else if (length(args)==1) {
    dataset_id=args[1]
}

require(data.table)
require(ggplot2)
require(tools)


BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR=paste(BASE_DIR, dataset_id, sep="/")

INDIR=paste(PROJECT_DIR, "analysis/output/fragment_size", sep="/")
OUTDIR=paste(PROJECT_DIR, "analysis/output/figures/fragment_size", sep="/")
if (!file.exists(OUTDIR)) {
 dir.create(OUTDIR)
}
fragment_files=list.files(INDIR, pattern="*_UMIs_fragment_size.txt", full.names=T)
print(fragment_files)

for (f in fragment_files) {
 # plot everything without filtering
 dt <- fread(f, col.names=c('Counts', 'Size'))
 plot <- ggplot(dt, aes(x=Size, y=Counts))
 plot <- plot + geom_bar(stat='identity', width=1, fill='#4a7dc4', col='gray40')+theme_bw()+ylab('Milion Reads')+xlab('Size of protected fragment(nt)')
ggsave(paste(OUTDIR, "/", file_path_sans_ext(basename(f)), ".pdf", sep=""), plot)
}
