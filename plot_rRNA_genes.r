#!/usr/bin/env Rscript

trna=F
args = commandArgs(trailingOnly=TRUE)
if (length(args) <2 ) {
    stop(paste("Usage example:", args[0], "14522_B6a", "trna/rrna", "[number_of_genes_per_plot]"))
} else {
    dataset_id=args[1]
    trna = args[2]
    if (trna == "trna") {
        trna=T
    } else {
        trna=F
    }
    if (length(args) == 3) {
        number_of_genes_per_plot=as.integer(args[3])
    } else {
        # this is only for trna genes, we have too many of them
        number_of_genes_per_plot=100
    }
}

require(ggplot2)
require(data.table)

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR=paste(BASE_DIR, dataset_id, sep="/")
if (trna) {
    INDIR=paste(PROJECT_DIR, "analysis/output/trna_fragments", sep="/")
    OUTDIR=paste(PROJECT_DIR, "analysis/output/figures/trna_fragments", sep="/")
} else {
    INDIR=paste(PROJECT_DIR, "analysis/output/rrna_fragments", sep="/")
    OUTDIR=paste(PROJECT_DIR, "analysis/output/figures/rrna_fragments", sep="/")
}

sample_files = list.files(INDIR, pattern="*_reads_per_gene.txt")

for (filename in sample_files) {
    print(paste("Reading", filename))
    sample_file = paste(INDIR, filename, sep="/")
    sample_name = gsub("_reads_per_gene.txt", "", filename)
    outfile = paste(OUTDIR, "/", sample_name, ".pdf", sep="")
    rrna_dt <- fread(sample_file, header=TRUE)
    if (trna) { 
        rrna_plot <- ggplot(rrna_dt, aes(x=gene_name, y=gene_counts))
        rrna_plot <- rrna_plot + geom_bar(stat='identity', width=.5, fill="tomato3") + labs(title=paste("Reads per gene for sample:", sample_name), subtitle=paste(nrow(rrna_dt), " genes. Total reads:", sum(rrna_dt$gene_counts))) + ylab("Reads") + xlab(NULL)
        rrna_plot <- rrna_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        ggsave(outfile, rrna_plot, width=20)
    } else {
        rrna_plot <- ggplot(rrna_dt, aes(x=gene_name, y=gene_counts/1000000))
        rrna_plot <- rrna_plot + geom_bar(stat='identity', width=.5, fill="tomato3") + labs(title=paste("Reads per gene for sample:", sample_name), subtitle=paste(nrow(rrna_dt), " genes. Total reads:", sum(rrna_dt$gene_counts))) + ylab("Milion reads") + xlab(NULL)

        ggsave(outfile, rrna_plot)
    }
    print(paste("Created plot", outfile))

}
