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

print("Please create analysis/output/cutadapt_plot_stats.txt file with the command ./utils/extract_cutadapt_stats.sh")
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

# Pipeline stats

# HQ Aligns
print("HQ Aligns")
mydt <- fread(paste(INDIR, '/alignment_hq_stats.txt', sep=""), col.names=c('file', 'HQ_algs'))
mydt[,file:=basename(file)]
mydt[,file:= gsub('_toGenome.hqmapped.bam',"\\1",file)]

# Dedup Aligns
print("Dedup Aligns")
dedup_file = paste(INDIR, 'alignment_dedup_stats.txt', sep="/")
if(!file.exists(dedup_file) | file.info(dedup_file)$size == 0) {
  mydt$Dedup_algs <- 0 
} else {
  dt <- fread(dedup_file, col.names=c('file', 'Dedup_algs'))
  dt[, file:=basename(file)]
  dt[,file:= gsub('_toGenome.hqmapped_dedup.bam',"\\1",file)]
  mydt <- merge(mydt,dt, by='file')
}

# Discarded alignments
print("Discarded alignments")
# dt <- fread(paste(INDIR, '/alignment_multimap_stats.txt', sep=""), col.names=c('file', 'N','MAPQ'))
# dt <- fread(paste(INDIR, '/alignment_multimap_stats.txt', sep=""), col.names=c('file', 'Multimap'))
# dt[, file:=basename(file)]
# dt[,file:= gsub('_toGenome.bam',"\\1",file)]
# dt <- dt[, .(seqnames, file, N)]
### WHAT IS < 50 ??
#dt <-  dt[,.(Multimappers=sum(N[MAPQ < 50])), by=.(file)]

dt <- fread(paste(INDIR, '/alignment_all_stats.txt', sep=""), col.names=c('file', 'No_HQ'))
dt[, file:=basename(file)]
dt[,file:= gsub('_toGenome.bam',"\\1",file)]

mydt <- merge(mydt,dt, by='file')

mydt$No_HQ <- mydt$No_HQ - mydt$HQ_algs

# rRNA cleanup
print("rRNA cleanup")
dt <- fread(paste(INDIR, '/RNA_clenup_stats.txt', sep=""))
dt[,file:= gsub('dem_(.*)_umi_extracted',"\\1",file)]
mydt <- merge(mydt,dt, by='file')

# Fix file name
print("Fix file name")
mydt[,file:= gsub('dem_(.*)_umi_extracted',"\\1",file)]

# Add initial reads (after bc split)
print("Add initial reads (after bc split)")
dt[,file:= gsub('(.*)_umi_extracted',"\\1",file)]
dt <- fread(BC_SPLIT_FILE,fill=T,select = c('Barcode','Count'))
# dt <- dt[1:nrow(mydt) -1]
colnames(dt) <- c('file','Initial_reads')
dt[,file:= gsub('(.*)_trimmed',"\\1",file)]
# mydt <- merge(mydt,dt, by.x = 0)
#mydt[["Initial_reads"]] = dt[["Initial_reads"]]
mydt <- merge(x = mydt, y = dt, by = "file", all.x = TRUE)

# Create stats
print("Create diricore stats")
# print(mydt)
#' mysum <- mydt[,.(
#'       'file'=file,
#'       'rRNA_reads'= Initial_reads - rrnaleft,
#'       'tRNA_reads' = rrnaleft - trnaleft,
#'       'Multimapper_reads' = Multimappers,
#'       #'Other_causes'= trnaleft - HQ_algs - Multimappers,
#'       'Dup_HQalgs' = HQ_algs - Dedup_algs,
#'       'Uniq_HQals' = Dedup_algs)]
mysum <- mydt[,.(
  'sample'=file,
  'rRNA_reads'= Initial_reads - rrnaleft,
  'tRNA_reads' = rrnaleft - trnaleft,
  'No_HQ' = No_HQ,
  #'Other_causes'= trnaleft - HQ_algs - Multimap,
  'HQ_with_duplicates' = HQ_algs - Dedup_algs,
  'Unique_HQ' = Dedup_algs)]

mysum$No_HQ <- mydt$Initial_reads - mysum$rRNA_reads - mysum$tRNA_reads - mysum$HQ_with_duplicates - mysum$Unique_HQ
mysum$Initial_reads <- mydt$Initial_reads
write.table(mysum, paste(INDIR, '/diricore_stats.txt',sep=""), sep='\t', quote=F, row.names=F)

mysum$Initial_reads <- NULL
# Check
rowSums(mysum[,-1]) == mydt$Initial_reads # True
# mysum <- mydt[,.('file'=file, 'value' = HQ_algs)]
# Arrange and save
print("Arrange and save")
mysum <- melt(mysum,id.vars='sample')

# Plot
print("Plot")
myplot <- ggplot(mysum, aes(x=sample, y=value/1000000,fill=variable))
myplot <- myplot + geom_bar(stat='identity')+ ggtitle("Unique HQ mapped reads") +
      theme_bw(8) + coord_flip() + ylab('Million reads') + xlab(NULL)
myplot
ggsave(paste(OUTDIR, '/Diricore_stats.pdf', sep=""), myplot,width = 5, height = 3)
