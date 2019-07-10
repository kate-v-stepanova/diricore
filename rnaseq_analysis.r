library(data.table)
library(DESeq2)
library(stats)
library(pheatmap)
library(ggplot2)

BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = paste(BASE_DIR, project_id, sep="/")
INDIR = paste(project_dir, "analysis/output/kallisto_quant", sep="/")
READS_DIR = paste(project_dir, "analysis/output/alignments/reads_per_gene", sep="/") # has to be not normalized!
contrasts_file = paste(project_dir, "analysis/input/metadata/rpf_density_contrasts.tsv", sep="/")


data_dir = paste(project_dir, "analysis/output/rna_seq", sep="/")
plot_dir = paste(project_dir, "analysis/output/figures/rna_seq", sep="/")

print("Creating output directories")
dir.create(plot_dir)
dir.create(data_dir)

print("Correlation matrix")
tpm = NULL
samples = c()
for (sample_dir in list.dirs(INDIR)) {
  if (sample_dir != INDIR) {
    samplename = basename(sample_dir)
    samples = c(samples, samplename)
    counts_file = paste(sample_dir, "abundance.tsv", sep="/")
    counts <- fread(counts_file, header = T, sep = "\t", select=c("target_id", "tpm"))
    counts$target_id <- vapply(strsplit(counts$target_id,"\\|"), `[`, 1, FUN.VALUE=character(1))
    counts$tpm <- as.numeric(counts$tpm)
    setnames(counts, "tpm", samplename)
    if (is.null(tpm)) {
      tpm = counts
    } else {
      tpm = merge(tpm, counts, by="target_id")
    } 
  }
}
tpm <- as.data.frame(tpm)
rownames(tpm) <- tpm$target_id
tpm$target_id <- NULL

countData <- as.matrix(tpm)
correlationMatrix <- cor(countData)
outfile = paste(plot_dir, "correlation_matrix.jpeg", sep="/")
print(paste("Writing correlation heatmap to:", outfile))
pheatmap(correlationMatrix, filename=outfile)


all_counts = NULL
for (counts_file in list.files(READS_DIR, pattern="*.out.tab")) {
  samplename = basename(counts_file)
  samplename = gsub("_ReadsPerGene.out.tab", "", counts_file)
  counts <- fread(paste(READS_DIR, counts_file, sep="/"), header = F, sep = "\t")
  counts <- tail(counts,-4)
  setnames(counts, c("transcript", samplename, "pos", "neg"))
  counts$pos <- NULL
  counts$neg <- NULL
  if (is.null(all_counts)) {
    all_counts = counts
  } else {
    all_counts = merge(all_counts, counts, by="transcript")
  } 
}
counts_df <- as.data.frame(all_counts)
rownames(counts_df) <- counts_df$transcript
counts_df$transcript <- NULL

countData <- as.matrix(counts_df)
colData = NULL
colData$v1 <- as.vector(samples)
colData <- as.data.frame(colData)
colData$v1 <- as.factor(colData$v1)

rownames(colData) = as.vector(samples)
designFormula <- colData
# not normalized counts
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = as.formula(designFormula))
dds <- dds[ rowSums(counts(dds)) > 1, ] # remove genes with 0 counts
dds <- DESeq(dds) # for normalization

DESeq2::plotMA(object = dds, ylim = c(-5, 5))

contrasts_dt <- fread(contrasts_file, header=F, sep="\t")
colnames(contrasts_dt) <- c("sample", "control", "color")
contrasts_dt$color <- NULL

print("Differential expression")
all_de_results = NULL
all_de_results$a <- counts_df$`1_16h_21perO2`
for(i in 1:nrow(contrasts_dt)) {
  row <- contrasts_dt[i,]
  contrast_name <- paste(row$sample, "__vs__", row$control, sep="")
  de_outfile = paste(data_dir, "/diff_expr_", contrast_name, ".tsv", sep="")
  DEresults = results(dds, contrast = c("v1", row$sample, row$control))
  DEresults <- DEresults[order(DEresults$pvalue),]
  de_results <- as.data.frame(DEresults)
  de_results$transcript <- rownames(de_results)
  print(paste("Writing file:", de_outfile))
  write.table(de_results, file=de_outfile, sep="\t", quote = F, row.names = F)
  de_results$contrast <- contrast_name
  # if (is.null(all_de_results)) {
  #   all_de_results = de_results
  # } else {
  #   all_de_results <- rbind(all_de_results, de_results)  
  # }
  all_de_results[[contrast_name]] = de_results
}
all_de_results$a <- NULL
#my_ma_plot <- 
de_plot =  ggplot(as.data.frame(de_results), aes(x=baseMean, y=log2FoldChange)) + ylim(-5, 5) + geom_point(size = 1)
outfile = paste(plot_dir, "diff_expr_all.jpeg", sep="/")
ggsave(outfile, de_plot)
print("Created file:", outfile)
# p <- DESeq2::plotMA(object = dds, ylim = c(-5, 5)) # alpha - threshold p-values, main - title

print("p-value distribution")
pplot = ggplot(data = as.data.frame(all_de_results), aes(x = pvalue)) + geom_histogram(bins = 100)
outfile = paste(plot_dir, "p_values.jpeg", sep="/")
ggsave(outfile, pplot)
print(paste("Created file:", outfile))

print("PCA")
rld <- rlog(dds)
my_pca_plot = plotPCA(rld, ntop = 500, intgroup = 'v1')
my_pca_data = plotPCA(rld, ntop = 500, intgroup = 'v1', returnData=T)
outfile = paste(data_dir, "pca_rlog.tsv", sep="/")
write.table(my_pca_data, outfile, quote = F, sep="\t", row.names = F)
outplot = paste(plot_dir, "pca_rlog.jpeg", sep = "/")
ggsave(outplot, my_pca_plot)


#countsNormalized <- counts(dds, normalized = TRUE)
# # select top 500 most variable genes
#selectedGenes <- names(sort(apply(countsNormalized, 1, var), decreasing = TRUE)[1:500])
#DESeq2::plotPCA(countsNormalized[selectedGenes,], col = as.numeric(colData$v1))

