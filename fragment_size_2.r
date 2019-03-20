#!/usr/bin/env Rscript

require(data.table)
require(ggplot2)
require(tools)

INDIR="./data/output/fragment_size"
OUTDIR="./data/output/figures/fragment_size"
if (!file.exists(OUTDIR)) {
 dir.create(OUTDIR)
}
fragment_files=list.files(INDIR, pattern="*_UMIs_fragment_size.txt", full.names=T)
print(fragment_files)

for (f in fragment_files) {
 dt <- fread(f, col.names=c('Counts', 'Size'))
 dt[Size >= 10, Size:=Size - 10L]
 dt <- dt[, .(Counts=sum(Counts)), by=Size]
 plot <- ggplot(dt, aes(x=Size, y=Counts))
 plot <- plot + geom_bar(stat='identity', width=1, fill='#4a7dc4', col='gray40')+theme_bw()+ylab('Milion Reads')+xlab('Size of protected fragment(nt)')
ggsave(paste(OUTDIR, "/", file_path_sans_ext(basename(f)), ".pdf", sep=""), plot)
}
