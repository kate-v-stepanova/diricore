library(data.table)
library(DESeq2)
library(stats)
library(pheatmap)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
genome_id=args[2]

# cpm - counts per million
BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = paste(BASE_DIR, project_id, sep="/")
READS_DIR = paste(project_dir, "analysis/output/alignments/reads_per_gene", sep="/") # has to be NOT normalized!
contrasts_file = paste(project_dir, "analysis/input/metadata/rpf_density_contrasts.tsv", sep="/")
gene_names = paste(BASE_DIR, "static", genome_id, "gene_names.txt", sep="/")

data_dir = paste(project_dir, "analysis/output/cpm_heatmap", sep="/")
plot_dir = paste(project_dir, "analysis/output/figures/cpm_heatmap", sep="/")

# we can pass a file which contains list of samples to be plotted (one line per sample)
if (length(args) >= 3) {
  infile_path <- args[3]
  selected_samples <- scan(infile_path, what=character())
  filename_preffix = args[4] # e.g. 1 or 2, so that the output files will be called 1_heatmap.csv
  filename_preffix = paste(filename_preffix, "_", sep="") # so that we separate preffix from samplename by "_"
} else {
  filename_preffix = ""
  selected_samples = c()
}


print("Creating output directories")
dir.create(plot_dir)
dir.create(data_dir)

# Getting table with gene_id and gene_name
gene_names_dt <- as.data.frame(fread(gene_names, sep="|", header=F))
colnames(gene_names_dt) <- c("transcript_id", "gene_id", "gene_name")
gene_names_dt <- gene_names_dt[!duplicated(gene_names_dt$gene_id), ]
gene_names_dt$transcript_id <- NULL # drop transcript id
rownames(gene_names_dt) <- gene_names_dt$gene_id


print("Correlation matrix")
cpm = NULL
samples = c()
coding_dt = NULL
noncoding_dt = NULL
for (counts_file in list.files(READS_DIR, full.names=T)) {
  samplename = basename(counts_file)
  samplename = gsub("_ReadsPerGene.out.tab", "", samplename)
  samples = c(samples, samplename)
  counts <- fread(counts_file, header = F, sep = "\t", skip=4)
  counts$V3 <- NULL
  counts$V4 <- NULL
  colnames(counts) <- c("gene_id", "reads")
  counts$reads <- as.numeric(counts$reads)
  # merge with gene_names
  counts = merge(counts, gene_names_dt, by="gene_id", all.x = T)
  
  # remove rows that don't have gene_name -> assuming they are non-coding
  coding <- counts[!is.na(counts$gene_name),]
  coding <- coding[, list(reads=sum(reads)),by=gene_name]
  rownames(coding) <- coding[["gene_name"]]
  names(coding) <- gsub("reads", samplename, names(coding))
  if(is.null(coding_dt)) {
    coding_dt <- coding
  } else {
    coding_dt <- merge(coding_dt, coding, by="gene_name")
  }
  
  # non-coding + coding
  counts$gene_name <- ifelse(is.na(counts$gene_name), counts$gene_id, counts$gene_name)
  counts <- counts[, list(reads=sum(reads)),by=gene_name]
  rownames(counts) <- counts$gene_name
  names(counts) <- gsub("reads", samplename, names(counts))
  if(is.null(noncoding_dt)) {
    noncoding_dt <- counts
  } else {
    noncoding_dt <- merge(noncoding_dt, counts, by="gene_name")
  }
}

if (length(selected_samples) == 0) {
  selected_samples = samples
}

# CODING
coding_cpm <- as.data.frame(coding_dt)
rownames(coding_cpm) <- coding_cpm$gene_name
coding_cpm$gene_name <- NULL
# drop rows where sum of counts in the row is less than 10 per sample
coding_cpm <- coding_cpm[as.logical(rowSums(coding_cpm >= length(samples) * 10)), ] 
# normalization to cpm (counts per million)
coding_cpm <- apply(coding_cpm, 2, function(x) x/sum(as.numeric(x)) * 10^6)

# Sort by variance (coding)
# 1 means rows (2 means columns, c(1, 2) means both), var - is a function computing variance across the data frame
V <- apply(coding_cpm, 1, var)
row_order = names(V[order(V, decreasing = T)])
coding_cpm <- coding_cpm[row_order,]
coding_cpm <- as.data.frame(coding_cpm)

# save dataframe (coding)
# add gene_name column for .tsv file
coding_cpm$gene_name <- rownames(coding_cpm)
data_file = paste(data_dir, "cpm_coding_genes.tsv", sep="/")
print(paste("Writing CPM for coding genes:", data_file))
write.table(coding_cpm, file=data_file, sep="\t", quote = F, row.names = F)

# save top 50 (coding)
top50_coding <- head(coding_cpm,50)

data_file = paste(data_dir, "cpm_top50_coding_genes.tsv", sep="/")
print(paste("Writing top 50 coding genes:", data_file))
write.table(top50_coding, file=data_file, sep="\t", quote = F, row.names = F)

# heatmap (coding)
top50_coding$gene_name <- NULL
top50_coding <- top50_coding[,selected_samples]

colData1 = NULL
colData1$samples <- as.vector(selected_samples)
colData1 <- as.data.frame(colData1)

colData1$samples <- as.factor(colData1$samples)
rownames(colData1) <- colnames(top50_coding)

outfile = paste(plot_dir, "/", filename_preffix, "heatmap_tpm_top50_coding.jpeg", sep="")
print(paste("Saving file:", outfile))

pheatmap(top50_coding, 
         scale = 'row',
         annotation_col = colData1, 
         show_rownames = T, filename=outfile, width = 16, height = 20, legend = F)

# NON-CODING
noncoding_cpm <- as.data.frame(noncoding_dt)
rownames(noncoding_cpm) <- noncoding_cpm$gene_name
noncoding_cpm$gene_name <- NULL
# drop rows where sum of counts in the row is less than 10 per sample
noncoding_cpm <- noncoding_cpm[as.logical(rowSums(noncoding_cpm >= length(samples) * 10)), ] 
# normalization to cpm (counts per million)
noncoding_cpm <- apply(noncoding_cpm, 2, function(x) x/sum(as.numeric(x)) * 10^6)

# Sort (non-coding)
# 1 means rows (2 means columns, c(1, 2) means both), var - is a function computing variance across the data frame
V <- apply(noncoding_cpm, 1, var)
row_order = names(V[order(V, decreasing = T)])
noncoding_cpm <- noncoding_cpm[row_order,]

# Save dataframe (non-coding)
noncoding_cpm <- as.data.frame(noncoding_cpm)
noncoding_cpm$gene_name <- rownames(noncoding_cpm)
data_file = paste(data_dir, "cpm_non_coding_genes.tsv", sep="/")
print(paste("Writing CPM for non-coding genes:", data_file))
write.table(noncoding_cpm, file=data_file, sep="\t", quote = F, row.names = F)
noncoding_cpm$gene_name <- NULL

# save top 50 (non-coding)
top50_noncoding <- head(noncoding_cpm,50)
data_file = paste(data_dir, "cpm_top50_noncoding_genes.tsv", sep="/")
print(paste("Writing top 50 non-coding genes:", data_file))
write.table(top50_noncoding, file=data_file, sep="\t", quote = F, row.names = F)


# heatmap (non-coding)
top50_noncoding$gene_name <- NULL
top50_noncoding <- top50_noncoding[,selected_samples]

colData1 = NULL
colData1$samples <- as.vector(selected_samples)
colData1 <- as.data.frame(colData1)

colData1$samples <- as.factor(colData1$samples)
rownames(colData1) <- colnames(top50_noncoding)

print(paste("Saving file:", outfile))
outfile = paste(plot_dir, "/", filename_preffix, "heatmap_tpm_top50_noncoding.jpeg", sep="")
print(paste("Saving file:", outfile))
pheatmap(top50_noncoding, 
         scale = 'row',
         annotation_col = colData1, 
         show_rownames = T, filename = outfile, width = 16, height = 20)
