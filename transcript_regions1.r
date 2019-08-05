#!/usr/bin/env Rscript

library(riboWaltz)
library(tools)
library(ggplot2)
library(data.table)
library(dplyr)


args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
genome=args[2]

base_dir = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = paste(base_dir, project_id, sep="/")
static_dir = paste(base_dir, "static", genome, sep="/")
transcriptome_fasta = paste(static_dir, "transcripts.fa", sep="/")
gtf_file = paste(static_dir, "gencode.annotation.gtf", sep="/")

bamfolder = paste(project_dir, "analysis/input/periodicity_bam", sep="/")
data_dir = paste(project_dir, "analysis/output/transcript_regions", sep="/")

if (!file.exists(bamfolder)) {
  stop(paste("No input files in ", bamfolder))
}

print("Creating output directory")
dir.create(data_dir)

print("Loading input files into memory. This might take a while")
annotation_dt = create_annotation(gtf_file)
reads_list <- bamtolist(bamfolder = bamfolder, annotation = annotation_dt)


samplenames = c()
for (filename in list.files(bamfolder)) {
  samplename = file_path_sans_ext(filename)
  samplenames <- c(samplenames, samplename)
}

print("Identifying the regions")
for (sample in samplenames) {
  current_dt <- reads_list[[sample]]
  current_dt <- current_dt[, .(transcript, end5, end3, cds_start, cds_stop)]
  current_dt <- current_dt %>%
    mutate(region = case_when(end5 > cds_start & end5 <= cds_stop  ~ 'CDS',
                              end5 <= cds_start & end5 >= cds_start ~ 'Start',
                              end5 < cds_start ~ "5' UTR",
                              end3 > cds_stop ~ "3' UTR",
                              TRUE ~ 'Undefined'))  
  dropped_dt <- current_dt[!(current_dt$cds_start == 0 & current_dt$cds_stop == 0),]
  filename = paste(data_dir, "/", sample, "_reads_list.tsv", sep="")
  print(paste("Writing:", filename))
  write.table(dropped_dt, file=filename, sep="\t", row.names = F, col.names = T, quote = F) 
  
  # ## The difference with my script is that they count non-coding genes as 3' ??
  # ## TODO: to test this, drop rows where cds_start and cds_stop is 0
  # #d<-d[!(d$A=="B" & d$E==0),]
  # dropped_dt <- current_dt[!(current_dt$cds_start == 0 & current_dt$cds_stop == 0),]
  # filename = paste(p_sites_per_region, "/", sample, "_dropped_reads_custom_list.tsv", sep="")
  # print(paste("Writing:", filename))
  # write.table(current_dt, file=filename, sep="\t", row.names = F, col.names = T, quote = F) 
  # ## Conclusion: there definitely is a difference.
}
print("Now aggregate samples in python script!")