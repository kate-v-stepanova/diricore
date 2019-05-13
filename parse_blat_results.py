#!/usr/bin/env python

import sys
import os
import glob
import pandas as pd

if len(sys.argv) <= 2:
    print("project_id is missing! ")
    print("Usage {} 14464 trna/rrna".format(sys.argv[0]))
    exit(1)

project_id = sys.argv[1]
trna = sys.argv[2]
if trna == "trna":
    trna = True
    prefix = "tRNA"
else:
    trna = False
    prefix = "rRNA"


BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
DIRICORE_PATH = "/home/e984a/diricore"

project_dir = os.path.join(BASE_DIR, project_id)
rrna_genes_file = os.path.join(DIRICORE_PATH, "staticdata/human/{}_genes.txt".format(prefix))

if trna:
    input_dir = os.path.join(project_dir, "analysis/output/trna_fragments")
else:
    input_dir = os.path.join(project_dir, "analysis/output/rrna_fragments")
print(input_dir)
# these files are manually created from UCSC browser
blat_results = glob.glob(os.path.join(input_dir, "*.parsed_psl.txt"))
rrna_genes_df = pd.read_csv(rrna_genes_file, sep="\t")
for infile in blat_results:
    print("Processing ", infile)
    output_file = infile.replace('.parsed_psl.txt', '_reads_per_gene.txt')
    # a line looks like: 345202_TTCGCGCGGGTCGGGGGGCGGGGCGGACTGT 18 8 25 31 100.0 chr6_cox_hap2 - 4700634 4700651 18 None None
    # we only need our id (which is counts+sequence), then chromosome, start and end positions
    df = pd.read_csv(infile, sep=" ", header=None, usecols=[0,6,8,9])
    df.columns = ['counts_sequence', 'chrom', 'start', 'end']

    # keep only the first hits (the sequences can be aligned to multiple regions, and blat outputs all of them)
    # "counts_sequence" column contains sequence which was aligned and the number of reads which match to this sequence.
    # Since all the duplicated sequences were merged, this column can be used as an identifier
#    df = df.drop_duplicates(subset="counts_sequence", keep="first")
    if df.empty:
         print("No data loaded! Skipping sample")
         continue
    # now I want to separate sequence and the number of reads
    df[['counts', 'sequence']] = df.get('counts_sequence').str.split('_', expand=True)
    result = []
    # should be just 3 genes, so not a big deal for performance
    for row_id, gene in rrna_genes_df.iterrows():
        #gene_name = row.get('gene_name')
        # now we extract from df only the hits which belong to this gene
        # we check that the start and end position of the hit is between start and end position of the gene
        gene_df = df.loc[(df.get('chrom') == gene.get('chrom')) & (df.get('start') >= gene.get('start')) & (df.get('end') <= gene.get('end'))]
        # now we can remove duplicates. Assuming there is no overlap betweeen the regions of rRNA genes, we can safely remove the duplicates so that each sequence will only be counted once
        gene_df = gene_df.drop_duplicates(subset="counts_sequence", keep="first")
        gene_name = gene.get('gene_name')
        # now just count total number of reads
        gene_counts = gene_df.get('counts').astype(int).sum()
        if int(gene_counts) != 0:
            result.append({
                'gene_name': gene_name,
                'gene_counts': gene_counts
            })
    result_df = pd.DataFrame(result)
    result_df.drop_duplicates(subset='gene_name')
#    df['Total'] = df.groupby(['Fullname', 'Zip'])['Amount'].transform('sum')
    # result_df['gene_counts'] = result_df.groupby(['gene_name'])['gene_counts'].transform('sum')
    result_df.to_csv(output_file, sep="\t", index=False)
    print("Created file: ", output_file)
