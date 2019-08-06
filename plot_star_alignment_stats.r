#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(RColorBrewer)


args = commandArgs(trailingOnly=TRUE)
 if (length(args)==0) {
    stop(paste("Usage example:", args[0], "14548"))
 } else if (length(args)==1) {
     dataset_id=args[1]
 }

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR=paste(BASE_DIR, dataset_id, sep="/")

INDIR=paste(PROJECT_DIR, "analysis/output/alignment_stats", sep="/")
OUTDIR=paste(PROJECT_DIR, "analysis/output/figures", sep="/")
BC_SPLIT_FILE=paste(PROJECT_DIR, "analysis/output/bc_split_stats.txt", sep="/")
cutadapt_file=paste(PROJECT_DIR, "analysis/output/cutadapt_plot_stats.txt", sep="/")
hq_file = paste(INDIR, "hq_stats.txt", sep="/")
hq_unique = paste(INDIR, "hq_unique_stats.txt", sep="/")
all_stats = paste(INDIR, "all_stats.txt", sep="/")
all_unique = paste(INDIR, "all_unique_stats.txt", sep="/")
rrna_file = paste(INDIR, "rrna_stats.txt", sep="/")
trna_file = paste(INDIR, "trna_stats.txt", sep="/")


print("Parsing cutadapt_plot_stats.txt")
lines <- read.csv(file=cutadapt_file, sep="\t", colClasses=c("NULL", NA), header=F)
no_adapt = lines[1,] - lines[2,]
too_short = strtoi(lines[3,])
passed = strtoi(lines[4,])
mydt <- data.table(
      'Reads' = c('No_adapt','Too_short','Passed'),
      'Counts' = c(no_adapt, too_short, passed))

myplot <- ggplot(mydt,aes(x=dataset_id,y=Counts/1000000, fill=Reads))
myplot <- myplot + geom_bar(stat='identity',position='stack')+
     theme_bw()+xlab(NULL)+ ylab('Milion Reads')
ggsave(paste(OUTDIR, '/Cutadapt_stats.pdf', sep=""), myplot, width = 2.5, height = 4)

# BC stats
print("BC stats")

mydt <- fread(BC_SPLIT_FILE, sep="\t", fill=T)
# mydt <- mydt[1:nrow(mydt) -1]
mydt <-mydt[!(mydt$Barcode=="total"),]
myplot <- ggplot(mydt, aes(x=dataset_id, y=Count/1000000, fill=Barcode)) + 
                geom_bar(stat='identity',position='stack')+
                scale_fill_manual(values=c(brewer.pal(10,'Set3'),'gray30', 'blue'))+
                xlab(NULL) + ylab('Million reads') + theme_bw() 
print("Saving BCsplit_stats.pdf")
ggsave(paste(OUTDIR, '/BCsplit_stats.pdf', sep=""), myplot, width = 2.5, height = 4)

# Alignment stats

# rRNA stats
print("rRNA stats")
mydt <- fread(rrna_file, col.names = c('sample', 'rrna'))

# tRNA
print("tRNA stats")
dt <- fread(trna_file, col.names = c('sample', 'trna'))
mydt <- merge(mydt,dt, by='sample')


# HQ unique
print("HQ unique")
if(!file.exists(hq_unique)) {
  mydt$hq_unique <- 0 
} else {
  dt <- fread(hq_unique, col.names=c('sample', 'hq_unique'))
  mydt <- merge(mydt,dt, by='sample')
}

# HQ with dup
print("HQ stats")
dt <- fread(hq_file, col.names=c('sample', 'hq_with_dup'))
mydt <- merge(mydt, dt, by="sample")
mydt$hq_with_dup <- mydt$hq_with_dup - mydt$hq_unique # hq_with_dup

# All unique
print("All unique")
if(!file.exists(all_unique)) {
  mydt$all_unique <- 0 
} else {
  dt <- fread(all_unique, col.names = c('sample', 'all_unique'))
  mydt <- merge(mydt,dt, by='sample')
}
mydt$lq_unique <- mydt$all_unique - mydt$hq_unique
mydt$all_unique <- NULL

# LQ with dup
print("LQ stats")
dt <- fread(all_stats, col.names=c('sample', 'all'))
mydt <- merge(mydt,dt, by='sample')
mydt$lq_with_dup <- mydt$all - mydt$hq_with_dup - mydt$hq_unique - mydt$lq_unique
mydt$all <- NULL


# total reads
print("total reads")
dt <- fread(BC_SPLIT_FILE, sep="\t", fill=T, col.names = c('sample', 'bc_split', 'file'))
dt$file <- NULL
mydt <- merge(mydt, dt, by='sample')

col_order = c("sample", "rrna", "trna", "lq_with_dup", "lq_unique", "hq_with_dup", "hq_unique", "bc_split")
mydt <- as.data.frame(mydt)[, col_order]

# Create stats
print("Create diricore stats")
write.table(mydt, paste(INDIR, '/diricore_stats.txt',sep=""), sep='\t', quote=F, row.names=F)

# Check
rowSums(mydt[,-1]) == mydt$bc_split # True
mydt$bc_split <- NULL
# Arrange and save
print("Arrange and save")
mysum <- melt(mydt,id.vars='sample')

# Plot
print("Plot")
myplot <- ggplot(mysum, aes(x=sample, y=value/1000000,fill=variable))
myplot <- myplot + geom_bar(stat='identity')+ ggtitle("Unique HQ mapped reads") +
      theme_bw(8) + coord_flip() + ylab('Million reads') + xlab(NULL)
ggsave(paste(OUTDIR, '/Diricore_stats.pdf', sep=""), myplot,width = 5, height = 3)
