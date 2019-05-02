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
    df.columns = ['counts_sequence', 'gene', 'start', 'end']
    df[['counts', 'sequence']] = df.get('counts_sequence').str.split('_', expand=True)
    df['gene'] = df['gene'].str.split('__').str[-1]
    df['read_length'] = df['sequence'].str.len()
    df['reads_info'] = ''

    if len(df) > 0:
        # should be just 3 genes, so not a big deal for performance
        for row_id, row in rrna_genes_df.iterrows():
            # now we extract from df only the hits which belong to this gene
            # we check that the start and end position of the hit is between start and end position of the gene
            gene_df = df.loc[df.get('gene') == row.get('gene_name')] #  & (df.get('start') >= row.get('start')) & (df.get('end') <= row.get('end'))]
            gene_name = row.get('gene_name')
            gene_df['counts'] = gene_df['counts'].astype('int')
            gene_df = gene_df.sort_values(by=['gene', 'start', 'end'])
            start_df = gene_df.groupby(['start', 'read_length'], as_index=False).agg({'counts': 'sum'})
            start_positions = start_df['start'].unique()
            for start in start_positions:
                all_reads_info = []
                for row_id, row in start_df.loc[start_df['start'] == start].iterrows():
                     all_reads_info.append( "{} reads of length {}".format(row['counts'], row['read_length']))
                gene_df.loc[gene_df['start'] == start, 'reads_info'] = ",".join(all_reads_info)
            gene_df = gene_df.groupby(['start'], as_index=False).agg({'counts': 'sum', 'gene': 'first', 'reads_info': 'first'})
            output_file = os.path.join(output_dir, '{}_{}_reads_per_position.txt'.format(sample, gene_name))
            if len(gene_df) != 0:
                gene_df.to_csv(output_file, sep="\t", index=False)
                print("Created file: ", output_file)
    else:
        print("Empty dataset. Skipping file: {}".format(infile))
