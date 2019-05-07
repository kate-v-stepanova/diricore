#!/usr/bin/env Rscript

require(ggplot2)
require(data.table)
require(RColorBrewer)


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

myplot <- ggplot(mydt,aes(x='run_umis',y=Counts/1000000, fill=Reads))
myplot <- myplot + geom_bar(stat='identity',position='stack')+
     theme_bw()+xlab(NULL)+ ylab('Milion Reads')
ggsave(paste(OUTDIR, '/Cutadapt_stats.pdf', sep=""), myplot, width = 2.5, height = 4)

# BC stats
print("BC stats")

mydt <- fread(BC_SPLIT_FILE, fill=T)
mydt <- mydt[1:nrow(mydt) -1]
myplot <- ggplot(mydt, aes(x='run_umis', y=Count/1000000, fill=Barcode))
myplot <- myplot + geom_bar(stat='identity',position='stack')+
                scale_fill_manual(values=c(brewer.pal(6,'Set3'),'gray30', 'blue'))+
                xlab(NULL) + ylab('Million reads') + theme_bw() 
print("Saving BCsplit_stats.pdf")
ggsave(paste(OUTDIR, '/BCsplit_stats.pdf', sep=""), myplot, width = 2.5, height = 4)

# Pipeline stats

# HQ Aligns
print("HQ Aligns")
mydt <- fread(paste(INDIR, '/alignment_hq_stats.txt', sep=""), col.names=c('file', 'HQ_algs'))
mydt[,file:=basename(dirname(file))]
# Dedup Aligns
print("Dedup Aligns")
dt <- fread(paste(INDIR, '/alignment_dedup_stats.txt', sep=""), col.names=c('file', 'Dedup_algs'))
dt[, file:=basename(dirname(file))]
mydt <- merge(mydt,dt, by='file')

# Discarded alignments
print("Discarded alignments")
dt <- fread(paste(INDIR, '/alignment_multimap_stats.txt', sep=""), col.names=c('file', 'N','MAPQ'))
### WHAT IS < 50 ??
dt <-  dt[,.(Multimappers=sum(N[MAPQ < 50])), by=.(file)]
mydt <- merge(mydt,dt, by='file')
# rRNA cleanup
print("rRNA cleanup")
dt <- fread(paste(INDIR, '/RNA_clenup_stats.txt', sep=""))
mydt <- merge(mydt,dt, by='file')

# Fix file name
print("Fix file name")
mydt[,file:= gsub('dem_(.*)_umi_extracted',"\\1",file)]

# Add initial reads (after bc split)
print("Add initial reads (after bc split)")
dt <- fread(BC_SPLIT_FILE,fill=T,select = c('Barcode','Count'))
# dt <- dt[1:nrow(mydt) -1]
colnames(dt) <- c('file','Initial_reads')
mydt <- merge(mydt,dt, by='file')

# Create stats
print("Create diricore stats")
# print(mydt)
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
