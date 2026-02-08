import csv
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("tsv_file")
ap.add_argument("--delimiter", type=str, default=None)
args = ap.parse_args()

delimiter = "\t" if args.delimiter is None else args.delimiter

with open(args.tsv_file, "r") as f:
    reader = csv.DictReader(f, delimiter=delimiter)
    rows = list([row for row in reader])
    print(f"read {len(rows)} rows")
    if len(rows) > 0:
        print(f"headers: {list(rows[0].keys())}")
