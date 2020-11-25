import pandas as pd
import os
import sys
import glob

project_id = sys.argv[1]

INDIR = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ribo_diff/te".format(project_id)

print("Processing files in: {}".format(INDIR))
for infile in glob.glob("{}/*.txt".format(INDIR)):
    df = pd.read_csv(infile, sep="\t")
    df = df.round(3)
    # drop column with all NaNs (it appeared for some reason)
    df = df.dropna(axis=1, how="all")
    df = df.dropna()
    df.to_csv(infile, sep="\t", index=False)    
