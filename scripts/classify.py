import csv
import gzip
from pathlib import Path
from needle.classify import classify

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("hmm_file")
ap.add_argument("module_id")
ap.add_argument("--disable-cutoff-ga", action="store_true", default=False)
args = ap.parse_args()

if "Pfam-A.hmm" in args.hmm_file and args.disable_cutoff_ga:
    print("warning: Pfam-A.hmm selected, recommend leaving on cutff-ga")
cutoff_ga = not args.disable_cutoff_ga
print("using cutoff_ga:", cutoff_ga)

# always load ko thresholds, if we use Pfam this will just get ignored

score_threshold_dict = None
with gzip.open("data/ko_thresholds.gz", mode='rt') as f:
    reader = csv.DictReader(f, delimiter='\t')
    score_threshold_dict = {row['knum']: row['threshold'] for row in reader}

hmm_file_stem = Path(args.hmm_file).stem
proteins_faa = f"data/{args.module_id}_results/proteins.faa"
proteins_tsv = f"data/{args.module_id}_results/proteins.tsv"
output_tsv = f"data/{args.module_id}_results/classify_{hmm_file_stem}.tsv"

classify(args.hmm_file, proteins_faa, proteins_tsv, cutoff_ga, output_tsv, score_threshold_dict)
