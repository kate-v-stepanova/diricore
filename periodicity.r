#!/usr/bin/env Rscript

library(riboWaltz)
library(tools)
library(ggplot2)
library(data.table)

print("DO NOT FORGET to sort bamfiles with the command: samtools sort sample.bam -o sample.sorted.bam")

args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
genome=args[2]

project_id = "3192"
base_dir = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = paste(base_dir, project_id, sep="/")
static_dir = paste(base_dir, "static", genome, sep="/")
transcriptome_fasta = paste(static_dir, "transcripts.fa", sep="/")
contrasts_file = paste(project_dir, "analysis/input/metadata/rpf_density_contrasts.tsv", sep="/")

gtf_file = paste(static_dir, "gencode.annotation.gtf", sep="/")

bamfolder = paste(project_dir, "analysis/input/periodicity_bam", sep="/")
plot_dir = paste(project_dir, "analysis/output/figures/periodicity", sep="/")
data_dir = paste(project_dir, "analysis/output/periodicity", sep="/")

heatmaps = paste(plot_dir, "heatmap", sep="/")
frames_path = paste(plot_dir, "frames", sep="/")
metaprofile_path = paste(plot_dir, "metaprofile", sep="/")
#codon_usage_path = paste(plot_dir, "codon_usage", sep="/")

print("Creating output directories")
dir.create(heatmaps)
dir.create(frames_path)
dir.create(metaprofile_path)

samplenames = c()
for (filename in list.files(bamfolder)) {
  samplename = file_path_sans_ext(filename)
  samplenames <- c(samplenames, samplename)
}

print("Loading input files into memory")
annotation_dt = create_annotation(gtf_file)
reads_list <- bamtolist(bamfolder = bamfolder, annotation = annotation_dt)

print("P-site offset")
psite_offset <- psite(reads_list, plot=F, plot_dir=plot_dir)
reads_psite_list <- psite_info(reads_list, psite_offset)

print("Read ends heatmap")
for (filename in list.files(bamfolder)) {
  print(filename)
  samplename = file_path_sans_ext(filename)
  ends_heatmap <- rends_heat(reads_list, annotation_dt, sample = samplename, cl = 85, utr5l = 25, cdsl = 40, utr3l = 25)
  output_file = paste(data_dir, "/", samplename, ".heatmap.tsv", sep="")
  write.table(ends_heatmap[["dt"]], file=output_file, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  plot = ends_heatmap[["plot"]] + scale_colour_gradient2(low="#bfbfbf", mid="#333333", high="black")
  ggsave(paste(heatmaps, "/", samplename, ".jpeg", sep=""), plot=plot, width=15, height=7, units="in")
}

print("P-site frames")
for (filename in list.files(bamfolder, pattern="*.bam")) {
  samplename = file_path_sans_ext(filename)
  output_file = paste(frames_path, "/", samplename, ".jpeg", sep="")
  frames <- frame_psite(reads_psite_list, sample = samplename, region = "all")
  ggsave(output_file, device = NULL, plot=frames[["plot"]], width=10, height=7, units="in")
}

print("P-site frames length")
for (filename in list.files(bamfolder, pattern="*.bam")) {
  samplename = file_path_sans_ext(filename)
  print(samplename)
  output_file = paste(frames_path, "/", samplename, ".length", ".jpeg", sep="")
  frames_stratified <- frame_psite_length(reads_psite_list, sample = samplename, region = "all", cl = 90)
  ggsave(output_file, device = NULL, plot=frames_stratified[["plot"]], width=10, height=7, units="in")
}

print("Metaprofile")
for (filename in list.files(bamfolder, pattern="*.bam")) {
  samplename = file_path_sans_ext(filename)
  output_file = paste(metaprofile_path, "/", samplename, ".jpeg", sep="")
  metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt, sample = samplename,
                                           utr5l = 20, cdsl = 40, utr3l = 20,
                                           plot_title = "auto")
  
  ggsave(output_file, device = NULL, plot=metaprofile[["plot"]], width=10, height=7, units="in")
}

# print("Codon usage barplot")
# # Codon usage barplot
# for (filename in list.files(bamfolder, pattern="*.bam")) {
#   samplename = file_path_sans_ext(filename)
#   print(samplename)
#   output_file = paste(codon_usage_path, "/", samplename, ".jpeg", sep="")
#   print(output_file)
#   codon_usage_barplot <- codon_usage_psite(reads_psite_list, annotation_dt, sample = samplename,
#                                            fastapath = transcriptome_fasta,
#                                            fasta_genome = FALSE,
#                                            frequency_normalization = FALSE) 
#   ggsave(output_file, device = NULL, plot=codon_usage_barplot[["plot"]], width=10, height=7, units="in")
# }
# 
# print("Codon usage pairwise (contrasts)")
# # Codon usage pairwise (contrasts)
# contrasts_dt <- read.table(file = contrasts_file, sep = '\t', header = FALSE)
# 
# for(i in 1:nrow(contrasts_dt)) {
#   row <- contrasts_dt[i,]
#   sample1 = paste(row["V1"], ".hqmapped.sorted", sep="")
#   sample2 = paste(row["V2"], ".hqmapped.sorted", sep="")
#   comparison_list <- list()
#   comparison_list[[sample1]] <- reads_psite_list[[sample1]]
#   comparison_list[[sample2]] <- reads_psite_list[[sample2]]
#   output_file = paste(codon_usage_path, "/", sample1, "_vs_", sample2, ".jpeg", sep="")
#   codon_usage_scatter <- codon_usage_psite(comparison_list, annotation_dt, 
#                                       sample = c(sample1, sample2),
#                                       fastapath = transcriptome_fasta,
#                                       gtfpath = gtf_file,
#                                       frequency_normalization = TRUE)
#   ggsave(output_file, device = NULL, plot=codon_usage_scatter[["plot"]], width=10, height=7, units="in")
# }

