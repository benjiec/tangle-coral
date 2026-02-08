import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("tsv")
parser.add_argument("genome_accession")
parser.add_argument("--when", default=None)
args = parser.parse_args()

FIELD = "genome_accession"

entries = []
with open(args.tsv, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if args.when:
            if FIELD in row and row[FIELD] == args.when:
                row[FIELD] = args.genome_accession
        else:
            row[FIELD] = args.genome_accession
        entries.append(row)

with open(args.tsv, "w") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(entries[0].keys()))
    writer.writeheader()
    for entry in entries:
        writer.writerow(entry)
