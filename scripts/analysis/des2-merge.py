# curates tall table of des2 data

import csv
import argparse
from pathlib import Path
from tangle.exp import DESeq2Table

parser = argparse.ArgumentParser()
parser.add_argument("experiment_id")
parser.add_argument("output_dir")
parser.add_argument("des2_tsv_files", nargs="+")
args = parser.parse_args()

ANALYSIS_FIELD = "analysis_type"

rows = []

def to_float(v):
    if v:
        return float(v)
    return "nan"

for tsv_file in args.des2_tsv_files:
    analysis_type = Path(tsv_file).stem
    print(analysis_type)
    with open(tsv_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row[ANALYSIS_FIELD] = analysis_type
            row["experiment_id"] = args.experiment_id
            row["baseMean"] = to_float(row["baseMean"]) 
            row["log2FoldChange"] = to_float(row["log2FoldChange"])
            row["lfcSE"] = to_float(row["lfcSE"])
            row["stat"] = to_float(row["stat"])
            row["pvalue"] = to_float(row["pvalue"])
            row["padj"] = to_float(row["padj"])
            row["max_cv"] = to_float(row["max_cv"])
            row["mean_base"] = to_float(row["mean_base"])
            row["mean_testgroup"] = to_float(row["mean_testgroup"])
            rows.append(row)

des2_merged_tsv = args.output_dir+"/des2_tall.tsv"
DESeq2Table.write_tsv(des2_merged_tsv, rows)
