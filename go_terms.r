library(DESeq2)
library(gProfileR)
library(data.table)

pro_file = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/pro_asp_leu/pro__geneid_top1000.tsv"
asp_file = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/pro_asp_leu/asp__geneid_top1000.tsv"
leu_file = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/pro_asp_leu/leu__geneid_top1000.tsv"

pro_out = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/go_terms/pro.tsv"
asp_out = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/go_terms/asp.tsv"
leu_out = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/go_terms/leu.tsv"

pro_out2 = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/go_terms/pro_without_genes.tsv"
asp_out2 = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/go_terms/asp_without_genes.tsv"
leu_out2 = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/go_terms/leu_without_genes.tsv"

pro_genes = fread(pro_file, sep="\t")
asp_genes = fread(asp_file, sep="\t")
leu_genes = fread(leu_file, sep="\t")


proResults <- gprofiler(query = pro_genes$gene_name,
                        organism = 'hsapiens', 
                        src_filter = 'GO', 
                        significant = F,
                        hier_filtering = 'moderate')
aspResults <- gprofiler(query = asp_genes$gene_name, 
                        organism = 'hsapiens', 
                        src_filter = 'GO', 
                        significant = F,
                        hier_filtering = 'moderate')
leuResults <- gprofiler(query = leu_genes$gene_name, 
                        organism = 'hsapiens', 
                        src_filter = 'GO', 
                        significant = F,
                        hier_filtering = 'moderate')

pro = proResults[order(proResults$p.value), c(3:4, 7, 10, 12, 14)]
asp = aspResults[order(aspResults$p.value), c(3:4, 7, 10, 12, 14)]
leu = leuResults[order(leuResults$p.value), c(3:4, 7, 10, 12, 14)]

write.table(pro, file=pro_out, sep="\t", quote = F, row.names = F)
write.table(asp, file=asp_out, sep="\t", quote = F, row.names = F)
write.table(leu, file=leu_out, sep="\t", quote = F, row.names = F)

pro2 = proResults[order(proResults$p.value), c(3:4, 7, 10, 12)]
asp2 = aspResults[order(aspResults$p.value), c(3:4, 7, 10, 12)]
leu2 = leuResults[order(leuResults$p.value), c(3:4, 7, 10, 12)]

write.table(pro2, file=pro_out2, sep="\t", quote = F, row.names = F)
write.table(asp2, file=asp_out2, sep="\t", quote = F, row.names = F)
write.table(leu2, file=leu_out2, sep="\t", quote = F, row.names = F)
