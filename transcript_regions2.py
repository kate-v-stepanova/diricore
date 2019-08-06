import pandas as pd
import os
import sys

project_id = sys.argv[1]
bam_type = "hq_unique"
bam_types = ["all", "all_unique", "hq"]
if len(sys.argv) >= 3:
    if sys.argv[2] in bam_types:
        bam_type = sys.argv[2]

basedir = "/icgc/dkfzlsdf/analysis/OE0532/"
project_dir = os.path.join(basedir, project_id)
input_dir = os.path.join(project_dir, "analysis/output/transcript_regions", bam_type)
outfile = os.path.join(input_dir, "reads_per_region.tsv")

files = os.listdir(input_dir)
files = list(filter(lambda x: x.endswith("_reads_list.tsv"), files))
samples = [f.replace("_reads_list.tsv", "") for f in files]
result = {}
for inf in files:
   sample = inf.replace("_reads_list.tsv", "")
   print(inf)
   df = pd.read_csv(os.path.join(input_dir, inf), sep="\t", header=0)

   # check there are no NAs in other columns than psite_region
   # df1 = df.drop(['psite_region'], axis=1)
   # len(df1.dropna()) == len(df)
   df = df.dropna()
   result[sample] = {
     "3' UTR": len(df.loc[df["region"] == "3' UTR"]),
     "CDS": len(df.loc[df["region"] == "CDS"]),
     "Start": len(df.loc[df["region"] == "Start"]),
     "5' UTR": len(df.loc[df["region"] == "5' UTR"]),
   }


full_df = pd.DataFrame(columns=sorted(result.keys()), index=("5' UTR", "Start", "CDS", "3' UTR"))
for sample in sorted(result.keys()): 
    full_df[sample] = pd.DataFrame({sample: result[sample]})
full_df.to_csv("reads_per_region_from_custom.tsv")

perc_df = pd.DataFrame(columns=["sample", "region", "reads", "percentage"])
for region in ["5' UTR", "Start", "CDS", "3' UTR"]:
    data = []
    for sample in full_df.columns:
        data.append({"sample": sample, "region": region, "reads": full_df[sample][region], "percentage": float(full_df[sample][region]) / float(full_df[sample].sum()) * 100})
    df = pd.DataFrame(data)
    if len(perc_df) == 0:
       perc_df = df
    else:
       perc_df = perc_df.append(df)
    
print("Writing: ", outfile)
perc_df.to_csv(outfile, sep="\t", index=False)

