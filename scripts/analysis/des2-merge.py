# curates tall table of des2 data, with a manifest

import re
import csv
import argparse
from pathlib import Path
from needle.seq import read_fasta_as_dict

parser = argparse.ArgumentParser()
parser.add_argument("output_dir")
parser.add_argument("faa_file")
parser.add_argument("des2_tsv_files", nargs="+")
args = parser.parse_args()

ANALYSIS_FIELD = "analysis_type"
seq_dict = read_fasta_as_dict(args.faa_file)
# strip off orifypy suffix on top of trinity suffix
seq_ids = { re.sub(r'_i\d+_ORF\.\d+$', '', k) for k in seq_dict.keys() }

entries = []

for tsv_file in args.des2_tsv_files:
    analysis_type = Path(tsv_file).stem
    print(analysis_type)
    with open(tsv_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row[ANALYSIS_FIELD] = analysis_type
            entries.append(row)

manifest_tsv = args.output_dir+"/sequence_list.tsv"
des2_merged_tsv = args.output_dir+"/des2_tall.tsv"

with open(manifest_tsv, "w") as f:
    # add an empty column to help Tableau understand this is a tab delimited file
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=["sequence_id", "empty"])
    writer.writeheader()
    for x in seq_ids:
        writer.writerow(dict(sequence_id=x, empty=""))

with open(des2_merged_tsv, "w") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(entries[0].keys()))
    writer.writeheader()
    for entry in entries:
        writer.writerow(entry)
