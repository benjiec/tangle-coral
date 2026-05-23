# curates tall table of des2 data

import csv
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("experiment_id")
parser.add_argument("output_dir")
parser.add_argument("des2_tsv_files", nargs="+")
args = parser.parse_args()

ANALYSIS_FIELD = "analysis_type"

entries = []

for tsv_file in args.des2_tsv_files:
    analysis_type = Path(tsv_file).stem
    print(analysis_type)
    with open(tsv_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row[ANALYSIS_FIELD] = analysis_type
            row["experiment_id"] = args.experiment_id
            entries.append(row)

des2_merged_tsv = args.output_dir+"/des2_tall.tsv"

with open(des2_merged_tsv, "w") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(entries[0].keys()))
    writer.writeheader()
    for entry in entries:
        writer.writerow(entry)
