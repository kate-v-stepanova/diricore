#!/usr/bin/env Rscript

trna=F
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop(paste("Usage example:", args[0], "14548", "trna/rrna"))
} else {
   dataset_id=args[1]
   trna=args[2]=="trna"
}
#if (length(args)==0) {
#  stop(paste("Usage example:", args[0], "14548"))
#} else if (length(args)==1) {
#   dataset_id=args[1]
#}

BASE_DIR="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_DIR=paste(BASE_DIR, dataset_id, sep="/")

if (trna) {
    input_dir=paste(PROJECT_DIR, "analysis/output/trna_fragments", sep="/")
} else {
    input_dir=paste(PROJECT_DIR, "analysis/output/rrna_fragments", sep="/")
}

require(data.table)
require(Biostrings)
require(ggplot2)

for (input_file in list.files(input_dir, pattern="*_top_tRNA_seqs.txt")) {
    print(paste("Processing", input_file))
    output_file=paste(input_dir, paste("grouped_", strsplit(input_file, '[.]')[[1]][1], ".txt", sep=""), sep="/")
    filename=paste(input_dir, input_file, sep="/")
    if (file.size(filename)==0) {
        print(paste("WARNING: file is empty. Skipping:", filename))
    } else {
        dt <- fread(filename, col.names=c('Counts', 'Seq'))

        #  Group the sequences using pairwiseAlignment( more flexible)
        print("Group the sequences using pairwise Alignment")
        seqs <- DNAStringSet(unique(dt$Seq))
        seqs <- seqs[order(width(seqs))]
        mylist <- list()
        while (length(seqs) > 0){
             p <- pairwiseAlignment(seqs, seqs[[length(seqs)]])
             id <- nmatch(p) / width(seqs) # calculate percentaje of identity
             g <- nmismatch(p) == 0 & id > 0.7
             mylist[as.character(seqs[[length(seqs)]])] <- list(as.character(seqs[g]))
             seqs <- seqs[!g]
        }  
        mygroups <- data.table(Group = unlist(sapply(1:length(mylist), function(x) rep(names(mylist[x]), length(mylist[[x]])))), Seq = unlist(mylist))

        # Add group info to counts table
        print("Add group info to counts table")
        mysum <- merge(dt,mygroups, by='Seq')
        mysum[,Counts_grouped:=sum(Counts),by=.(Group)]
        write.table(mysum[order(-Counts_grouped)], output_file, sep='\t', quote=F,row.names=F)
        print(paste("Created file", output_file))
   }
}
