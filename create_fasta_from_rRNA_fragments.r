#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop(paste("Usage example:", args[0], "14548"))
} else if (length(args)==1) {
  dataset_id=args[1]
}

trna=F
if (length(args)<2) {
   stop(paste("Usage example:", args[0], "14548", "trna/rrna"))
} else {
   dataset_id=args[1]
   trna=args[2]
   if (trna == "trna") {
      trna=T
   } else {
      trna=F
   }
}

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR=paste(BASE_DIR, dataset_id, sep="/")
if (trna) {
    input_dir=paste(PROJECT_DIR, "analysis/output/trna_fragments", sep="/")
} else {
    input_dir=paste(PROJECT_DIR, "analysis/output/rrna_fragments", sep="/")
}

for (filename in list.files(input_dir, pattern="grouped_")) {
    sample_name=strsplit(filename, '[.]')[[1]][1]
    sample_name=gsub("grouped_dem_", "", sample_name)
    sample_name=gsub("_umi_extracted_", "", sample_name)
    output_file=paste(input_dir, "/", sample_name, ".fasta", sep="")
    input_file=paste(input_dir, filename, sep="/")
    grouped_table <- read.table(input_file, header=TRUE)
    grouped_table = grouped_table[, c("Group", "Counts_grouped")]
    unique_rows = unique(grouped_table)
    print(unique_rows)
    print(filename)
    unique_rows$fasta <- paste(">", unique_rows$Counts_grouped, "_", unique_rows$Group, sep="")
    unique_rows = unique_rows[, c("fasta", "Group")]
    write.table(unique_rows, output_file, sep="\n", quote=F, row.names=F, col.names=F)
    print(paste("Created file", output_file))
}

