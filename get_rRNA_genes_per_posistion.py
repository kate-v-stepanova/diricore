#!/usr/bin/env python

import sys
import os
import glob
import pandas as pd

if len(sys.argv) < 2:
    print("project_id is missing! ")
    print("Usage {} 14464 trna/rrna [hg19]".format(sys.argv[0]))
    exit(1)

project_id = sys.argv[1]
trna = sys.argv[2]
if trna == "trna":
    trna = True
    prefix = "tRNA"
else:
    trna = False
    prefix = "rRNA"

species='hg19'
if len(sys.argv) >= 4:
   species = sys.argv[3]

BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
DIRICORE_PATH = "/home/e984a/diricore"
STATIC_PATH = os.path.join(BASE_DIR, 'static', species)

project_dir = os.path.join(BASE_DIR, project_id)
rrna_genes_file = os.path.join(STATIC_PATH, "{}_genes.txt".format(prefix))

if trna:
    input_dir = os.path.join(project_dir, "analysis/output/trna_fragments")
    output_dir = os.path.join(project_dir, "analysis/output/trna_positions")
else:
    input_dir = os.path.join(project_dir, "analysis/output/rrna_fragments")
    output_dir = os.path.join(project_dir, "analysis/output/rrna_positions")

if not os.path.isdir(output_dir):
     os.makedirs(output_dir)

blat_results = glob.glob(os.path.join(input_dir, "*.parsed_psl.txt"))
rrna_genes_df = pd.read_csv(rrna_genes_file, sep="\t")
chromosomes = rrna_genes_df['chrom'].unique().tolist()
for infile in blat_results:
    print("Processing ", infile)
    sample = os.path.basename(infile).replace('.parsed_psl.txt', '')
    # a line looks like: 345202_TTCGCGCGGGTCGGGGGGCGGGGCGGACTGT 18 8 25 31 100.0 chr6_cox_hap2 - 4700634 4700651 18 None None
    # we only need our id (which is counts+sequence), then chromosome, start and end positions
    df = pd.read_csv(infile, sep=" ", header=None, usecols=[0,6,8,9])
    df.columns = ['counts_sequence', 'chrom', 'start', 'end']

    # filter results (keep only the hits that aligned to chromosomes of rRNA genes)
    df = df[df.get('chrom').isin(chromosomes)]
    print(len(df))
    # keep only the first hits (the sequences can be aligned to multiple regions, and blat outputs all of them)
    # "counts_sequence" column contains sequence which was aligned and the number of reads which match to this sequence.
    # Since all the duplicated sequences were merged, this column can be used as an identifier
#    df = df.drop_duplicates(subset="counts_sequence", keep="first")

    # now I want to separate sequence and the number of reads
    if len(df) > 0:
        df[['counts', 'sequence']] = df.get('counts_sequence').str.split('_', expand=True)
        # should be just 3 genes, so not a big deal for performance
        for row_id, row in rrna_genes_df.iterrows():
            print(row['gene_name'])            
            #gene_name = row.get('gene_name')
            # now we extract from df only the hits which belong to this gene
            # we check that the start and end position of the hit is between start and end position of the gene
            gene_df = df.loc[(df.get('chrom') == row.get('chrom')) & (df.get('start') >= row.get('start')) & (df.get('end') <= row.get('end'))]
            print(len(gene_df))
            # now we can remove duplicates. Assuming there is no overlap betweeen the regions of rRNA genes, we can safely remove the duplicates so that each sequence will only be counted once
            # gene_df = gene_df.drop_duplicates(subset="counts_sequence", keep="first")
            gene_name = row.get('gene_name')
            gene_df['gene'] = gene_name
            gene_df['counts'] = gene_df['counts'].astype('int')
            gene_df = gene_df.sort_values(by=['chrom', 'start', 'end'])
            gene_df = gene_df.groupby(['chrom','start', 'end', 'gene'],as_index=False).sum() #agg({'counts': 'sum'})
            output_file = os.path.join(output_dir, '{}_{}_reads_per_positions.txt'.format(sample, gene_name))
            if len(gene_df) != 0:
                gene_df.to_csv(output_file, sep="\t", index=False)
                print("Created file: ", output_file)
    else:
        print("Empty dataset. Skipping file: {}".format(infile))
