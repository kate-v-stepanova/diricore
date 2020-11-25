#!/usr/bin/env Rscript

### R - 3.6.2
library(data.table)
library(riboWaltz)
library(tools)
library(ggplot2)
# INSTALL_opts = c('--no-lock')

print("DO NOT FORGET to sort bamfiles with the command: samtools sort sample.bam -o sample.sorted.bam")

args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
genome="hg19" #args[2]
# genome = 'mm9'

bam_type = "all_unique"
if (length(args) >= 3) {
  bam_type = args[3] # can be: hq, hq_unique, all, all_unique
}

base_dir = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = paste(base_dir, project_id, sep="/")
static_dir = paste(base_dir, "static", genome, sep="/")
transcriptome_fasta = paste(static_dir, "transcripts_single_header.fa", sep="/")
genome_fasta = paste(static_dir, "genome.fa", sep="/")
contrasts_file = paste(project_dir, "analysis/input/metadata/rpf_density_contrasts.tsv", sep="/")

gtf_file = paste(static_dir, "gencode.annotation.gtf", sep="/")
# gtf_file = paste(static_dir, "Mt_tRNA.gtf", sep="/")
#gtf_file = paste(static_dir, "protein_coding.annotation.gtf", sep="/")
# gtf_file = paste(static_dir, "Saccharomyces_cerevisiae.R64-1-1.98.gtf", sep="/")

bamfolder = paste(project_dir, "analysis/input/periodicity_bam", bam_type, sep="/")
plot_dir = paste(project_dir, "analysis/output/figures/periodicity", bam_type, sep="/")
data_dir = paste(project_dir, "analysis/output/periodicity", bam_type, sep="/")

heatmaps = paste(plot_dir, "heatmap", sep="/")
frames_path = paste(plot_dir, "frames", sep="/")
metaprofile_path = paste(plot_dir, "metaprofile", sep="/")
length_path = paste(plot_dir, "read_lengths", sep="/")
psite_offset_dir = paste(plot_dir, "psite_offset", sep="/")
psite_region_dir = paste(plot_dir, "psite_region", sep="/")
psite_region_data_dir = paste(data_dir, "psite_region", sep="/")
#codon_usage_path = paste(plot_dir, "codon_usage", sep="/")

print("Creating output directories")
dir.create(heatmaps, recursive=T, showWarnings = F)
dir.create(frames_path, recursive=T, showWarnings = F)
dir.create(metaprofile_path, recursive=T, showWarnings = F)
dir.create(psite_offset_dir, recursive=T, showWarnings = F)
dir.create(psite_region_dir, recursive=T, showWarnings = F)
dir.create(psite_region_data_dir, recursive=T, showWarnings = F)
dir.create(paste(data_dir, "psites", sep="/"), recursive = T, showWarnings = F)
dir.create(paste(plot_dir, "psites", sep="/"), recursive = T, showWarnings = F)
dir.create(paste(data_dir, "asites", sep="/"), recursive = T, showWarnings = F)
dir.create(paste(plot_dir, "asites", sep="/"), recursive = T, showWarnings = F)
dir.create(length_path, recursive=T, showWarnings = F)

samplenames = c()
for (filename in list.files(bamfolder)) {
  samplename = file_path_sans_ext(filename)
  samplenames <- c(samplenames, samplename)
}

print("Loading input files into memory")
annotation_dt = create_annotation(gtf_file)
annotation_dt = annotation_dt[annotation_dt$l_cds != 0]
# write.table(annotation_dt, file="/icgc/dkfzlsdf/analysis/OE0532/static/mm9/annotation_dt_ribowaltz.tsv", sep="\t", quote=F, row.names=F, col.names = T)

# reads_list <- bamtolist(bamfolder = bamfolder, annotation = annotation_dt, transcript_align = F)
reads_list <- bamtolist(bamfolder = bamfolder, annotation = annotation_dt, transcript_align = T)
# remove non-coding
for (sample in samplenames) {
  reads_list[[sample]] <- reads_list[[sample]][reads_list[[sample]]$cds_stop != 0]
}

print('P-site stats')
ribo_stats = paste(data_dir, '/ribowaltz_stats.tsv', sep="") # created manually from the output of reads_list()
# header: sample  input_reads not_in_annotation   negative_strand output_reads    cds_reads   inframe_reads
# example: /icgc/dkfzlsdf/analysis/OE0532/18927/analysis/output/periodicity/all_unique/ribowaltz_stats.tsv
stats_df = as.data.frame(fread(ribo_stats, sep="\t", header = T, fill=T))
stats_df$cds_reads = 0
stats_df$inframe_reads = 0
for (sample in samplenames) {
  cur_psite = reads_psite_list[[sample]]
  cds_df = cur_psite[complete.cases(cur_psite),]
  cds_reads = nrow(cds_df)
  inframe_df = cds_df[cds_df$psite_from_start %% 3 == 0,]
  inframe_reads = nrow(inframe_df)
  stats_df[stats_df$sample == sample,]$cds_reads = cds_reads
  stats_df[stats_df$sample == sample,]$inframe_reads = inframe_reads
}
print(paste('Writing stats:', ribo_stats))
write.table(stats_df, file=ribo_stats, quote = F, sep="\t", row.names = F)

print("P-site offset")
psite_offset <- psite(reads_list, extremity="auto", plot=F, plot_dir=psite_offset_dir)
reads_psite_list <- psite_info(reads_list, psite_offset)
# Getting rid of transcripts which are not in fasta (otherwise RiboWaltz fails)
# Error: subscript contains invalid names
# sequences <- Biostrings::readDNAStringSet(transcriptome_fasta, format = "fasta", use.names = TRUE)
# name_seq <- names(sequences)
name_seq = unique(annotation_dt$transcript)
length(name_seq)
name_tr <- unique(reads_psite_list[[samplename]]$transcript)
length(name_tr)
table(name_tr %in% name_seq)
for (filename in list.files(bamfolder)) {
  samplename = file_path_sans_ext(filename)
  print(samplename)
  data_for_current_sample = reads_psite_list[[samplename]] #
  data_for_current_sample = data_for_current_sample[data_for_current_sample$transcript %in% name_seq]
  data_for_current_sample = data_for_current_sample[data_for_current_sample$transcript %in% name_tr]
  reads_psite_list1 = NULL
  reads_psite_list1[[samplename]] = data_for_current_sample
  ## if get this error:
  ##  -> "Error in .Call2("solve_user_SEW", refwidths, start, end, width, translate.negative.coord,  : 
  ##  ->   solving row 318282: 'allow.nonnarrowing' is FALSE and the supplied start (0) is < 1""
  ##    -> do this:
  # reads_psite_list1[[samplename]] = reads_psite_list1[[samplename]][-c(318282)]
  
  # p-site
  psite_codons <- codon_usage_psite(reads_psite_list1, annotation_dt, sample = samplename,
                                    fastapath = genome_fasta,
                                    fasta_genome = T,
                                    # frequency_normalization = FALSE,
                                    gtfpath = gtf_file)
  # data
  datafile = paste(data_dir, "/psites/", samplename, "_psite_codon_counts.tsv", sep="")
  print(paste("Writing:", datafile))
  write.table(psite_codons$dt, file=datafile, sep="\t", quote=FALSE, row.names = F, col.names = T)
  # plot
  imgfile = paste(plot_dir, "/psites/", samplename, "_psite_codon_counts.jpeg", sep="")
  print(paste("Saving plot:", imgfile))
  ggsave(imgfile, device = NULL, plot=psite_codons$plot, width=15, height=7, units="in")
}


# print("A-site")
# # remove from reads list
# asite_offset <- psite(reads_list, plot=F, extremity="auto")
# reads_asite_list <- psite_info(reads_list, asite_offset, 
#                                site='asite', 
#                                fastapath=genome_fasta, fasta_genome = T,
#                                gtfpath = gtf_file)
# sequences <- Biostrings::readDNAStringSet(transcriptome_fasta, format = "fasta", use.names = TRUE)
# name_seq <- names(sequences)
# length(name_seq)
# name_tr <- unique(reads_psite_list[[samplename]]$transcript)
# length(name_tr)
# table(name_tr %in% name_seq)
# for (filename in list.files(bamfolder)) {
#   samplename = file_path_sans_ext(filename)
#   print(samplename)
# 
#   # a-site
#   asite_codons <- codon_usage_psite(reads_asite_list, annotation_dt, sample = samplename,
#                                     fastapath = genome_fasta,
#                                     frequency_normalization = FALSE, site="asite")
#   # data
#   datafile = paste(data_dir, "/asites/", samplename, "_asite_codon_counts.tsv", sep="")
#   print(paste("Writing:", datafile))
#   write.table(asite_codons$dt, file=datafile, sep="\t", quote=FALSE, row.names = F, col.names = T)
#   # plot
#   imgfile = paste(plot_dir, "/asites/", samplename, "_asite_codon_counts.jpeg", sep="")
#   print(paste("Saving plot:", imgfile))
#   ggsave(imgfile, device = NULL, plot=asite_codons$plot, width=15, height=7, units="in")
# }

print("Read ends heatmap")
for (filename in list.files(bamfolder)) {
  print(filename)
  samplename = file_path_sans_ext(filename)
  ends_heatmap <- rends_heat(reads_list, annotation_dt, sample = samplename, cl = 85, utr5l = 25, cdsl = 40, utr3l = 25)
  output_file = paste(data_dir, "/", samplename, ".heatmap.tsv", sep="")
  write.table(ends_heatmap[["dt"]], file=output_file, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  plot = ends_heatmap[["plot"]] + scale_colour_gradient2(low="#bfbfbf", mid="#333333", high="black")
  plot_file = paste(heatmaps, "/", samplename, ".jpeg", sep="")
  print(paste('Writing', plot_file))
  ggsave(plot_file, plot=plot, width=15, height=7, units="in")
}
  
print("P-site frames")
for (filename in list.files(bamfolder, pattern="*.bam")) {
  samplename = file_path_sans_ext(filename)
  print(samplename)
  output_file = paste(frames_path, "/", samplename, ".jpeg", sep="")
  frames <- frame_psite(reads_psite_list, sample = samplename, region = "all")
  print(paste("Writing", output_file))
  ggsave(output_file, device = NULL, plot=frames[["plot"]], width=10, height=7, units="in")
}

print("P-site frames length")
for (filename in list.files(bamfolder, pattern="*.bam")) {
  samplename = file_path_sans_ext(filename)
  print(samplename)
  frames_stratified <- frame_psite_length(reads_psite_list, sample = samplename, region = "all", cl = 90)
  output_file = paste(frames_path, "/", samplename, ".length", ".jpeg", sep="")
  print(paste("Writing", output_file))
  ggsave(output_file, device = NULL, plot=frames_stratified[["plot"]], width=10, height=7, units="in")
}

print("Metaprofile")
for (filename in list.files(bamfolder, pattern="*.bam")) {
  samplename = file_path_sans_ext(filename)
  print(samplename)
  ## Error in parse(text = x) : <text>:1:4: unexpected input
  ## 1: 159_
  ##        ^
  ## Means: sample name cannot start with a number.
  ## To solve: rename files in analysis/input/periodicity_bam
  ## Example: for f in $(ls *); do fn=${f#"231_"}; echo "mv $f $fn"; done
  ## And get again reads_psite_list
  metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt, sample = samplename,
                                   utr5l = 20, cdsl = 40, utr3l = 20,
                                   plot_title = samplename)
  output_file = paste(metaprofile_path, "/", samplename, ".jpeg", sep="")
  print(paste("Writing", output_file))
  ggsave(output_file, device = NULL, plot=metaprofile[[paste("plot_", samplename, sep="")]], width=10, height=7, units="in")
}

print("Read lengths distribution")
for (sample in samplenames) {
  example_length_dist <- rlength_distr(reads_list, sample = sample)
  outfile = paste(length_path, "/", sample, ".jpeg", sep="")
  print(paste("Writing", outfile))
  ggsave(outfile, device = NULL, plot=example_length_dist[["plot"]])
}
print('Psite regions')
for (sample in samplenames) {
  psite_region <- region_psite(reads_psite_list, annotation_dt, sample = sample)
  plot_file = paste(psite_region_dir, "/", sample, ".pdf", sep="")
  data_file = paste(psite_region_data_dir, "/", sample, ".tsv", sep="")
  print(paste('Writing data:', data_file))
  write.table(psite_region$dt, file=data_file, sep="\t", quote=FALSE, row.names = F, col.names = T)
  print(paste('Writing image file: ', plot_file))
  ggsave(plot_file, psite_region$plot)
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

