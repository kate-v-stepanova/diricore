library(data.table)
# library(DESeq2)
library(stats)
library(pheatmap)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
genome_id= "hg19"
bam_type = args[2]

# cpm - counts per million
BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = paste(BASE_DIR, project_id, sep="/")
# READS_DIR = paste(project_dir, "analysis/output/alignments/reads_per_gene", sep="/") # has to be NOT normalized!
READS_DIR = paste(project_dir, "analysis/output/alignments/reads_per_gene/tsv", bam_type, sep="/")
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
# 
# selected_samples = c("Mix_mC", "CoC_mC", "159_7_mC")
# selected_samples = c("CoC_GFP", "Mix_GFP", "MRC5_G_GFP")
# 
# MRC5_genes2 = c("THBS1", "COL1A2", "FN1", "AKR1B1", "ACTB", "LGALS1", "FAT1", "ALDOA", "ITGB1", "PTX3", "TNFRSF11B", "PLAU", "GAPDH", "NPTX1", "MMP1", "MYL9", "HMOX1", "SERPINE1", "LAMC1", "CANX", "HSP90B1", "TFRC", "ADAMTS1", "MKI67", "STAT1", "ATP5B", "GPI", "HSP90AB1", "WDR1", "LOXL2", "MZF1", "SEMA3A", "RPL11", "PLOD2", "CTNNB1", "HSPG2", "WLS", "PRKDC", "G3BP1", "FAM129B", "MAP1A", "RPL10A", "CPA4", "SEC61A1", "ATP2A2", "AXL", "RPS4Y1", "FBLN1", "APLP2", "DDOST", "NPTN", "MYOF", "SCARB2", "HNRNPD", "HSPA13", "ALDH18A1", "LMAN1", "PALLD", "PGD", "ARL6IP5", "CAPNS1", "CITED2", "PICALM", "YME1L1", "HSPA9", "RPL3", "UGGT1", "XPO1", "HK1", "VDAC3", "POMP", "EZR")
# genes_159 = c("AKR1B1", "LTBP1", "LDHA", "TIMP1", "ALDOA", "PXDN", "CD44", "VIM", "HGF", "GAPDH", "RND3", "FLNC", "NT5E", "NRP1", "CYR61", "COL12A1", "FST", "COL6A1", "HIF1A", "STAT1", "EFEMP1", "THY1", "CRIM1", "PLAT", "QSOX1", "DKK1", "ATP1A1", "KIAA1199", "HIST1H2AG", "SLC16A3", "SRGN", "UGDH", "CORO1C", "G3BP1", "CAPRIN1", "COL5A1", "LAMB1", "PLOD1", "RPL10A", "RPL38", "EIF4G2", "CPA4", "DDB1", "PSMB1", "MT-ND2", "ANPEP", "PRPF8", "PPAP2B", "FXYD5", "AXL", "TFPI2", "LAMP2", "EIF3B", "CD9", "DDOST", "CD164", "PDIA6", "SNRPB", "TIMP2", "TRIP12", "ESYT1", "SSR2", "THRAP3", "MCM7", "G3BP2", "RPN1", "CPD", "SERINC1", "H2AFZ", "MYL12A", "PSMA7", "PSMB6", "SNRNP200", "PICALM", "GNAI3", "NAV2", "OS9", "GNAS", "DEGS1", "MCAM", "CHD4")

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
# # drop rows where sum of counts in the row is less than 10 per sample
# coding_cpm <- coding_cpm[as.logical(rowSums(coding_cpm >= length(samples) * 10)), ] 

controls = c('159_7_mC', 'MRC5_G_GFP')
coding_cpm <- coding_cpm[!apply(coding_cpm[,controls], 1, function(row) all(row<=10)),]

# normalization to cpm (counts per million)
coding_cpm <- apply(coding_cpm, 2, function(x) x/sum(as.numeric(x)) * 10^6)

# Sort by variance (coding)
# 1 means rows (2 means columns, c(1, 2) means both), var - is a function computing variance across the data frame
V <- apply(coding_cpm, 1, var)
row_order = names(V[order(V, decreasing = T)])
# coding_cpm <- coding_cpm[sorted_names_coding,] # fold_changes_from_counts.r
coding_cpm <- coding_cpm[row_order,]
coding_cpm <- coding_cpm[MRC5_genes2,]
coding_cpm <- coding_cpm[genes_159,]
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
# top50_coding <- top50_coding[,selected_samples]
coding_cpm <- coding_cpm[,selected_samples]

colData1 = NULL
colData1$samples <- as.vector(selected_samples)
colData1 <- as.data.frame(colData1)

colData1$samples <- as.factor(colData1$samples)
# rownames(colData1) <- colnames(top50_coding)
rownames(colData1) <- colnames(coding_cpm)

outfile = paste(plot_dir, "/", filename_preffix, "heatmap_tpm_top50_coding.jpeg", sep="")
outfile = paste(plot_dir, "/", "159", "heatmap_cpm_selected_genes_coding.jpeg", sep="")
print(paste("Saving file:", outfile))
pheatmap(top50_coding, 
         scale = 'row',
         annotation_col = colData1, 
         show_rownames = T, filename=outfile, width = 16, height = 20, legend = F)

pheatmap(coding_cpm, 
         scale = 'row',
         annotation_col = colData1, 
         show_rownames = T, filename=outfile, width = 16, height = 20, legend = F)

# NON-CODING
noncoding_cpm <- as.data.frame(noncoding_dt)
rownames(noncoding_cpm) <- noncoding_cpm$gene_name
noncoding_cpm$gene_name <- NULL

# # drop rows where sum of counts in the row is less than 10 per sample
# noncoding_cpm <- noncoding_cpm[as.logical(rowSums(noncoding_cpm >= length(samples) * 10)), ] 

controls = c('159_7_mC', 'MRC5_G_GFP')
noncoding_cpm <- noncoding_cpm[!apply(noncoding_cpm[,controls], 1, function(row) all(row<=10)),]

# normalization to cpm (counts per million)
noncoding_cpm <- apply(noncoding_cpm, 2, function(x) x/sum(as.numeric(x)) * 10^6)

# Sort (non-coding)
# 1 means rows (2 means columns, c(1, 2) means both), var - is a function computing variance across the data frame
V <- apply(noncoding_cpm, 1, var)
row_order = names(V[order(V, decreasing = T)])
noncoding_cpm <- noncoding_cpm[row_order,]
noncoding_cpm <- noncoding_cpm[genes_159,]
noncoding_cpm <- noncoding_cpm[MRC5_genes2,]
noncoding_cpm <- noncoding_cpm[selected_genes,]


# Save dataframe (non-coding)
noncoding_cpm <- as.data.frame(noncoding_cpm)
noncoding_cpm <- noncoding_cpm[,selected_samples]
noncoding_cpm$gene_name <- rownames(noncoding_cpm)
data_file = paste(data_dir, "cpm_non_coding_genes.tsv", sep="/")
data_file = paste(data_dir, "MRC5_cpm_non_coding_genes.tsv", sep="/")
data_file = paste(data_dir, "159_cpm_non_coding_genes.tsv", sep="/")
print(paste("Writing CPM for non-coding genes:", data_file))
write.table(noncoding_cpm, file=data_file, sep="\t", quote = F, row.names = F)
noncoding_cpm$gene_name <- NULL

# save top 50 (non-coding)
top50_noncoding <- head(noncoding_cpm,50)
data_file = paste(data_dir, "cpm_top50_noncoding_genes.tsv", sep="/")
print(paste("Writing top 50 non-coding genes:", data_file))
write.table(top50_noncoding, file=data_file, sep="\t", quote = F, row.names = F)


# heatmap (non-coding)
# top50_noncoding$gene_name <- NULL
# top50_noncoding <- top50_noncoding[,selected_samples]
noncoding_cpm$gene_name <- NULL
noncoding_cpm <- noncoding_cpm[,selected_samples]

colData1 = NULL
colData1$samples <- as.vector(selected_samples)
colData1 <- as.data.frame(colData1)

colData1$samples <- as.factor(colData1$samples)
# rownames(colData1) <- colnames(top50_noncoding)
rownames(colData1) <- colnames(noncoding_cpm)

print(paste("Saving file:", outfile))
outfile = paste(plot_dir, "/", filename_preffix, "heatmap_cpm_selected_genes.jpeg", sep="")
outfile = paste(plot_dir, "/", "159_heatmap_cpm_selected_genes.jpeg", sep="")
print(paste("Saving file:", outfile))
pheatmap(top50_noncoding, 
         scale = 'row',
         annotation_col = colData1, 
         show_rownames = T, filename = outfile, width = 16, height = 20)
pheatmap(noncoding_cpm, 
         scale = 'row',
         annotation_col = colData1, 
         show_rownames = T, filename = outfile, width = 16, height = 20)