#!/usr/bin/env Rscript

input_file="./data/output/rrna_fragments/run_UMIs_top_rRNA_seqs.txt"
output_file="data/output/rrna_fragments/run_UMIs_top_rRNA_seqs_grouped.txt"

require(data.table)
require(Biostrings)
require(ggplot2)

dt <- fread(input_file,col.names=c('Counts','Seq'))


#  Group the sequences using pairwiseAlignment( more flexible)
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
mysum <- merge(dt,mygroups, by='Seq')
mysum[,Counts_grouped:=sum(Counts),by=.(Group)]
write.table(mysum[order(-Counts_grouped)], output_file, sep='\t', quote=F,row.names=F)
