#!/usr/bin/env Rscript

# parsing arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop(paste("Usage example:", args[0], "14548"))
} else if (length(args)==1) {
  dataset_id=args[1]
}

# constant paths
BASE_PATH="/icgc/dkfzlsdf/analysis/OE0532"
METADATA_PATH="analysis/input/metadata"

dataset_path=paste(BASE_PATH, dataset_id, sep="/")
infile=paste(dataset_path, "/", dataset_id, "_meta.tsv", sep="")
outfile=paste(dataset_path, METADATA_PATH, "bc_file.txt", sep="/")

# reading .tsv file
bc <- read.table(infile, sep="\t", header=TRUE)
bc <- bc[, c("SAMPLE_ID", "BARCODE")]
bc <- bc[with(bc, !grepl("Undetermined", BARCODE)),]
bc <- as.data.frame.matrix(bc)
bc$SAMPLE_ID <- gsub(" ", "_", bc$SAMPLE_ID)

# writing bc_file
write.table(bc, file=outfile, sep="\t", col.names=F, row.names=F, quote=F)
print(paste("Created file: ", outfile))


