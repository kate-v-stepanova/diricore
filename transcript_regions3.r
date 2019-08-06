#!/usr/bin/env Rscript

library(riboWaltz)
library(tools)
library(ggplot2)
library(data.table)
library(dplyr)


args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
genome=args[2]
bam_type = "hq_unique"
bam_types = c("hq", "all", "all_unique")
if (length(args) >= 3) {
  if (args[3] %in% bam_types) {
    bam_type = args[3]
  }
}

base_dir = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = paste(base_dir, project_id, sep="/")
data_dir = paste(project_dir, "analysis/output/transcript_regions", bam_type, sep="/")
input_file = paste(data_dir, "reads_per_region.tsv", sep="/")

plot_dir = paste(project_dir, "analysis/output/figures/transcript_regions", bam_type, sep="/")
plot_file = paste(plot_dir, "/", project_id, "_transcript_regions.pdf", sep="")
percentage_plot = paste(plot_dir, "/", project_id, "_transcript_regions_percentage.pdf", sep="")

print("Creating output directory")
dir.create(plot_dir)

print(paste("Reading input file:", input_file))
dt <- read.csv(input_file, sep="\t")

print(paste("Creating plot:", plot_file))
colors = c("gray85", "gray70", "gray50", "gray32", "gray10")
p <- ggplot(dt, aes(x=sample, y=reads, fill=factor(region, levels=c("5' UTR", "Start", "CDS", "3' UTR")))) +
  geom_bar(stat='identity',color = "white",width = 0.65, alpha = 0.9, size = 0.025 * 25) +
  scale_fill_manual(name = "", values = colors)
ggsave(plot_file, device = NULL, plot=p)

print(paste("Creating plot with percentages:", percentage_plot))
p <- ggplot(dt, aes(x=sample, y=percentage, fill=factor(region, levels=c("5' UTR", "Start", "CDS", "3' UTR")))) +
  geom_bar(stat='identity',color = "white",width = 0.65, alpha = 0.9, size = 0.025 * 25) +
  scale_fill_manual(name = "", values = colors)
ggsave(percentage_plot, device = NULL, plot=p)
