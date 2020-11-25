library(coRdon)

args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
genome=args[2]

BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = paste(BASE_DIR, project_id, sep="/")
indir = paste(project_dir, "analysis/output/alignments/toGenome", sep="/")

for (input_file in list.files(indir, pattern="*.fasta", full.names=T)) {
  reads_df <- readSet(file=input_file)
  codons_df <- codonTable(reads_df)
  codonCounts(codons_df)
}