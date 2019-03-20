#!/bin/bash

gtf='./staticdata/human/gencode.v29.basic.annotation.gtf'
files=$(ls -1 ./data/output/tophat_out/*/accepted_hits.hqmapped_dedup.bam | tr '\n' ' ')
htseq-count -f bam -r pos -s yes -t exon -i gene_id --additional-attr=gene_name --additional-attr=gene_type -m union $files $gtf > ./data/output/htseq_counts_dedup.txt

# Add the header (NOTE the double quotes to pass var to sed)
# Note, the first step removes a trailing space
files=$(echo $files | sed 's/[[:blank:]]+$//')
header="Gene_id\tGene_name\tGene_type\t${files// /\\t}"
sed -i "1i$header" ./data/output/htseq_counts_dedup.txt

# repeat the same steps using all the HQ algs (before deduping)
files=$(ls -1 data/output/tophat_out/*/accepted_hits.hqmapped.bam | tr '\n' ' ')
htseq-count -f bam -r pos -s yes -t exon -i gene_id --additional-attr=gene_name --additional-attr=gene_type -m union $files $gtf > ./data/output/htseq_counts.txt

# Add the header (NOTE the double quotes to pass var to sed)
# Note, the first step removes a trailing space
files=$(echo $files | sed 's/[[:blank:]]+$//')
header="Gene_id\tGene_name\tGene_type\t${files// /\\t}"
sed -i "1i$header" ./data/output/htseq_counts.txt


