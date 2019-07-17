library(data.table)
library(DESeq2)
library(stats)
library(pheatmap)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
genome_id=args[2]

BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = paste(BASE_DIR, project_id, sep="/")
READS_DIR = paste(project_dir, "analysis/output/alignments/reads_per_gene", sep="/") # has to be NOT normalized!
contrasts_file = paste(project_dir, "analysis/input/metadata/rpf_density_contrasts.tsv", sep="/")
gene_names = paste(BASE_DIR, "static", genome_id, "gene_names.txt", sep="/")
ma_dir = paste(project_dir, "analysis/output/ma_plot", sep="/")
ma_plot_dir = paste(project_dir, "analysis/output/figures/ma_plot", sep="/")

data_dir = paste(project_dir, "analysis/output/fold_change", sep="/")
plot_dir = paste(project_dir, "analysis/output/figures/fold_change", sep="/")

print("Creating output directories")
dir.create(plot_dir)
dir.create(data_dir)
dir.create(ma_dir)
dir.create(ma_plot_dir)

# GENE NAMES
gene_names_dt <- as.data.frame(fread(gene_names, sep="|", header=F))
colnames(gene_names_dt) <- c("transcript_id", "gene_id", "gene_name")
#gene_names_dt$transcript_id <- NULL # discard
gene_names_dt <- gene_names_dt[!duplicated(gene_names_dt$gene_id), ]
rownames(gene_names_dt) <- gene_names_dt$gene_id

# CONTRASTS
contrasts_dt <- fread(contrasts_file, header=F, sep="\t")
colnames(contrasts_dt) <- c("sample", "control", "color")
contrasts_dt$color <- NULL

# COUNTS - from STAR
all_counts = NULL
samples = c()
for (counts_file in list.files(READS_DIR, pattern="*.out.tab")) {
  samplename = basename(counts_file)
  samplename = gsub("_ReadsPerGene.out.tab", "", counts_file)
  samples = c(samples, samplename)
  counts <- fread(paste(READS_DIR, counts_file, sep="/"), header = F, sep = "\t")
  counts <- tail(counts,-4)
  setnames(counts, c("gene_id", samplename, "pos", "neg"))
  counts$pos <- NULL
  counts$neg <- NULL
  if (is.null(all_counts)) {
    all_counts = counts
  } else {
    all_counts = merge(all_counts, counts, by="gene_id")
  } 
}
counts_df <- as.data.frame(all_counts)
rownames(counts_df) <- counts_df$gene_id
counts_df$gene_id <- NULL


countData <- as.matrix(counts_df)
countData <- countData[!apply(countData[,samples], 1, function(row) all(row<=30)),]
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


# CODING AND NON-CODING GENES
all_de_results = NULL
all_de_results$a <- counts_df[[samplename]]
only_coding_dt = NULL
only_coding_dt$a <- counts_df[[samplename]]
list_of_contrasts = c()
for(i in 1:nrow(contrasts_dt)) {
  row <- contrasts_dt[i,]
  contrast_name <- paste(row$sample, "__vs__", row$control, sep="")
  list_of_contrasts <- c(list_of_contrasts, contrast_name)
  DEresults = results(dds, contrast = c("v1", row$sample, row$control))
  DEresults <- DEresults[order(DEresults$pvalue),]
  de_results <- as.data.frame(DEresults)
  
  # Replace gene_id with gene_name
  de_results$gene_id <- rownames(de_results)
  de_results = merge(de_results, gene_names_dt, by="gene_id", all.x = T)
  
  # drop non-coding genes (without name)
  coding_dt <- de_results[!is.na(de_results$gene_name),]
  # duplicated genes concatenate with the transcript_id
  coding_dt[duplicated(coding_dt$gene_name),]$gene_name <- paste(coding_dt[duplicated(coding_dt$gene_name),]$gene_name, " - ", coding_dt[duplicated(coding_dt$gene_name),]$gene_id, sep="")
  
  coding_dt$transcript_id <- NULL
  coding_dt <- as.data.frame(coding_dt)
  rownames(coding_dt) <- coding_dt$gene_name
  only_coding_file = paste(ma_dir, "/ma_plot_coding_", contrast_name, ".tsv", sep="")
  print(paste("Writing file:", only_coding_file))
  write.table(coding_dt, file=only_coding_file, sep="\t", quote = F, row.names = F)
  
  
  de_results$gene_trans <- de_results$gene_name
  de_results[!is.na(de_results$gene_name),]$gene_name <- paste(de_results[!is.na(de_results$gene_name),]$gene_name, " - ", de_results[!is.na(de_results$gene_name),]$gene_id, sep="")
  de_results[!duplicated(de_results$gene_trans),]$gene_name <- de_results[!duplicated(de_results$gene_trans),]$gene_trans
  de_results[is.na(de_results$gene_name),]$gene_name <- de_results[is.na(de_results$gene_name),]$gene_id
  
  de_results$transcript_id <- NULL
  de_results$gene_trans <- NULL
  de_results <- as.data.frame(de_results)
  rownames(de_results) <- de_results$gene_name
  
  de_outfile = paste(ma_dir, "/ma_plot_all_", contrast_name, ".tsv", sep="")
  print(paste("Writing file:", de_outfile))
  write.table(de_results, file=de_outfile, sep="\t", quote = F, row.names = F)
  
  all_de_results[[contrast_name]] <- de_results
  only_coding_dt[[contrast_name]] <- coding_dt
  
  outfile = paste(ma_plot_dir, "/ma_plot_all_", contrast_name, ".jpeg", sep="")
  print(paste("Saving plot:", outfile))
  my_ma_plot <-  ggplot(as.data.frame(de_results), aes(x=baseMean, y=log2FoldChange)) + ylim(-5, 5) + geom_point(size = 1)
  ggsave(outfile, my_ma_plot)
  
  outfile = paste(ma_plot_dir, "/ma_plot_coding_", contrast_name, ".jpeg", sep="")
  print(paste("Saving plot:", outfile))
  my_ma_plot <-  ggplot(as.data.frame(coding_dt), aes(x=baseMean, y=log2FoldChange)) + ylim(-5, 5) + geom_point(size = 1)
  ggsave(outfile, my_ma_plot)
}
all_de_results$a <- NULL
only_coding_dt$a <- NULL

## ALL FOLD CHANGES INTO ONE TABLE
coding_fc <- NULL
all_fc_df <- NULL
for (contrast_name in list_of_contrasts) {
  # ALL
  con_df <- all_de_results[[contrast_name]]$log2FoldChange
  con_df <- as.data.frame(con_df)
  rownames(con_df) <- all_de_results[[contrast_name]]$gene_name
  con_df$con_df <- NULL
  con_df[[contrast_name]] <- all_de_results[[contrast_name]]$log2FoldChange
  if (is.null(all_fc_df)) {
    all_fc_df <- con_df
  } else {
    all_fc_df[[contrast_name]] <- con_df[[contrast_name]]
  }
  
  # CODING
  con_df <- NULL
  con_df <- only_coding_dt[[contrast_name]]$log2FoldChange
  con_df <- as.data.frame(con_df)
  con_df$con_df <- NULL
  rownames(con_df) <- only_coding_dt[[contrast_name]]$gene_name
  con_df[[contrast_name]] <- only_coding_dt[[contrast_name]]$log2FoldChange
  if (is.null(coding_fc)) {
    coding_fc <- con_df
  } else {
    coding_fc[[contrast_name]] <- con_df[[contrast_name]]
  }
}

# SORT AND WRITE FC OF ALL GENES
sorted_names <- names(sort(apply(all_fc_df, 1, var), decreasing = T)) # 1 means rows (2 means columns, c(1, 2) means both), 
                                                                      # var - is a function computing variance across the data frame
all_fc_df <- all_fc_df[sorted_names,]
all_fc_df$gene_name <- rownames(all_fc_df)
col_order <- c("gene_name", list_of_contrasts)
all_fc_df <- all_fc_df[, col_order]
fc_data_file = paste(data_dir, "/FC_all_genes.tsv", sep="")
print(paste("Writing FC data file:", fc_data_file))
write.table(all_fc_df, file=fc_data_file, sep="\t", quote = F, row.names = F)


# SORT AND WRITE FC OF CODING GENES
sorted_names <- names(sort(apply(coding_fc, 1, var), decreasing = T)) # 1 means rows (2 means columns, c(1, 2) means both), 
                                                                      # var - is a function computing variance across the data frame
coding_fc <- coding_fc[sorted_names,]
coding_fc$gene_name <- rownames(coding_fc)
col_order <- c("gene_name", list_of_contrasts)
coding_fc <- coding_fc[, col_order]
coding_data_file = paste(data_dir, "/FC_coding_genes.tsv", sep="")
print(paste("Writing FC data file:", coding_data_file))
write.table(coding_fc, file=coding_data_file, sep="\t", quote = F, row.names = F)

# TOP 50 OF ALL GENES
all_fc_df <- head(all_fc_df,50)
all_fc_df$gene_name <- rownames(all_fc_df)
col_order <- c("gene_name", list_of_contrasts)
all_fc_df <- all_fc_df[, col_order]
fc_data_file = paste(data_dir, "/FC_top_50_all.tsv", sep="")
print(paste("Writing FC data file:", fc_data_file))
write.table(all_fc_df, file=fc_data_file, sep="\t", quote = F, row.names = F)
all_fc_df$gene_name <- NULL

colData1 = NULL
colData1$contrasts <- as.vector(list_of_contrasts)
colData1 <- as.data.frame(colData1)

colData1$contrasts <- as.factor(colData1$contrasts)
rownames(colData1) <- colnames(all_fc_df)

print(paste("Saving file:", outfile))
outfile = paste(plot_dir, "FC_top50_all_genes.jpeg", sep="/")
print(paste("Saving file:", outfile))
pheatmap(all_fc_df, 
         scale = 'row',
         annotation_col = colData1, 
         show_rownames = T, filename = outfile, width = 16, height = 20)



# TOP 50 OF CODING GENES
coding_fc <- head(coding_fc,50)
coding_fc$gene_name <- rownames(coding_fc)
col_order <- c("gene_name", list_of_contrasts)
coding_fc <- coding_fc[, col_order]
fc_data_file = paste(data_dir, "/FC_top_50_coding.tsv", sep="")
print(paste("Writing FC data file:", fc_data_file))
write.table(coding_fc, file=fc_data_file, sep="\t", quote = F, row.names = F)
coding_fc$gene_name <- NULL

colData1 = NULL
colData1$contrasts <- as.vector(list_of_contrasts)
colData1 <- as.data.frame(colData1)

colData1$contrasts <- as.factor(colData1$contrasts)
rownames(colData1) <- colnames(coding_fc)

print(paste("Saving file:", outfile))
outfile = paste(plot_dir, "FC_top50_coding_genes.jpeg", sep="/")
print(paste("Saving file:", outfile))
pheatmap(coding_fc, 
         scale = 'row',
         annotation_col = colData1, 
         show_rownames = T, filename = outfile, width = 16, height = 20)
