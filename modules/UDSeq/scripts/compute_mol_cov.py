import gzip
import glob
import argparse
import os

parser = argparse.ArgumentParser("compute_mol_cov.py")

parser.add_argument("-i", "--input")
parser.add_argument("-o", "--output")

args = parser.parse_args()

files = glob.glob(os.path.join(args.input, "*coverage.bed.gz"))
size = 0
cov = 0
for file in files:
    #print(file)
    with gzip.open(file, mode = "rt") as f:
        for line in f:
            #print(line)
            line = line.strip().split("\t")
            size += 1
            cov += min(int(line[3]), 1)
        #print(f"Completed {file}, processed {size} bp with {cov} covered bases")

with open(args.output, "w") as f:
    f.write("total_mol_bases\tsnv_mol_bases\tpct_covered\n")
    f.write(f"{size}\t{cov}\t{cov/3137300923}\n")