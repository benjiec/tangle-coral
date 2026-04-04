import csv
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("sequence_ko_tsv")
args = ap.parse_args()

with open(args.sequence_ko_tsv, "r") as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = [row for row in reader]

# filter by score
rows = [row for row in rows if row["score_threshold"] != "-"]
rows = [row for row in rows if float(row["dom_score"]) >= float(row["score_threshold"])]

with open(args.sequence_ko_tsv+"_filtered", "w") as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=list(rows[0].keys()))
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
