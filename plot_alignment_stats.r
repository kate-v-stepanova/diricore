#!/usr/bin/env Rscript

require(ggplot2)
require(data.table)
require(RColorBrewer)

INDIR="./data/output/alignment_stats"
OUTDIR="./data/output/figures"
BC_SPLIT_FILE="./data/output/bc_split_stats.txt"

# Manually input data from Logs/cutadapt_trimming_stats.txt
print("Manually input data from Logs/cutadapt_trimming_stats.txt")

mydt <- data.table(
      'Reads' = c('No_adapt','Too_short','Passed'),
      'Counts' = c(263519465-233487367, 140206555,93280812))

myplot <- ggplot(mydt,aes(x='run_umis',y=Counts/1000000, fill=Reads))
myplot <- myplot + geom_bar(stat='identity',position='stack')+
      theme_bw()+xlab(NULL)+ ylab('Million reads')
ggsave(paste(OUTDIR, '/Cutadapt_stats.pdf', sep=""), myplot, width = 2.5, height = 4)

# BC stats
print("BC stats")

mydt <- fread(BC_SPLIT_FILE, fill=T)
# mydt <- fread('Logs/bc_split_stats.txt',fill=T)
mydt <- mydt[1:nrow(mydt) -1]
myplot <- ggplot(mydt, aes(x='run_umis', y=Count/1000000, fill=Barcode))
myplot <- myplot + geom_bar(stat='identity',position='stack')+
                scale_fill_manual(values=c(brewer.pal(6,'Set3'),'gray30', 'blue'))+
                xlab(NULL) + ylab('Million reads') + theme_bw() 
print('ggsave doesnt work')
ggsave(paste(OUTDIR, '/BCsplit_stats.pdf', sep=""), myplot, width = 2.5, height = 4)

# Pipeline stats

# HQ Aligns
print("HQ Aligns")
mydt <- fread(paste(INDIR, '/alignment_hq_stats.txt', sep=""), col.names=c('file', 'HQ_algs'))
mydt[,file:=dirname(file)]

# Dedup Aligns
print("Dedup Aligns")
dt <- fread(paste(INDIR, '/alignment_dedup_stats.txt', sep=""), col.names=c('file', 'Dedup_algs'))
dt[,file:=dirname(file)]
mydt <- merge(mydt,dt, by='file')

# Discarded alignments
print("Discarded alignments")
dt <- fread(paste(INDIR, '/alignment_multimap_stats.txt', sep=""), col.names=c('file', 'N','MAPQ'))
dt <-  dt[,.(Multimappers=sum(N[MAPQ < 50])), by=.(file)]
mydt <- merge(mydt,dt, by='file')

# rRNA cleanup
print("rRNA cleanup")
dt <- fread(paste(INDIR, '/RNA_cleup_stats.txt', sep=""))
mydt <- merge(mydt,dt, by='file')

# Fix file name
print("Fix file name")
# print(file)
write.table(mydt, 'debug.txt',sep='\t', quote=F, row.names=F)

# print(mydt[,file:=gsub('run_umi_(.*)_umied', "\\1", file)])
# mydt[,file:= gsub('run_umi_(.*)_umied',"\\1",file)]

# Add initial reads (after bc split)
print("Add initial reads (after bc split)")
dt <- fread(BC_SPLIT_FILE,fill=T,select = c('Barcode','Count'))
dt <- dt[1:nrow(mydt) -1]
colnames(dt) <- c('file','Initial_reads')
mydt <- merge(mydt,dt, by='file')

# Create stats
print("Create stats")
mysum <- mydt[,.(
      'file'=file,
      'rRNA_reads'= Initial_reads - rrnaleft,
      'tRNA_reads' = rrnaleft - trnaleft,
      'Multimapper_reads' = Multimappers,
      'Other_causes'= trnaleft - HQ_algs - Multimappers,
      'Dup_HQalgs' = HQ_algs - Dedup_algs,
      'Uniq_HQals' = Dedup_algs)]

# Check
rowSums(mysum[,-1]) == mydt$Initial_reads # True

# Arrange and save
print("Arrange and save")
mysum <- melt(mysum,id.vars='file')
write.table(mysum, paste(INDIR, '/diricore_stats.txt',sep=""), sep='\t', quote=F, row.names=F)

# Plot
print("Plot")
myplot <- ggplot(mysum, aes(x=file, y=value/1000000,fill=variable))
myplot <- myplot + geom_bar(stat='identity')+
      theme_bw(8) + coord_flip() + ylab('Million reads') + xlab(NULL)
myplot
ggsave(paste(OUTDIR, '/Diricore_stats.pdf', sep=""), myplot,width = 5, height = 3)
