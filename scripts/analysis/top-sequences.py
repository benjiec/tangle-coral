import re
import argparse
from collections import defaultdict
from tangle.exp import DESeq2Table, TranscriptGenesTable
from tangle.detected import DetectedTable
from tangle.models import CSVSource

ap = argparse.ArgumentParser()
ap.add_argument("des2_tall_fn")
ap.add_argument("--mean-threshold", default=200, type=float)
ap.add_argument("--l2fc-threshold", default=1, type=float)
ap.add_argument("--padj-threshold", default=0.05, type=float)
args = ap.parse_args()

source = CSVSource(DESeq2Table, args.des2_tall_fn)
des2_rows = source.values()
des2_rows = [row for row in des2_rows if row["sequence_type"] == "gene"]

filtered_rows = [r for r in des2_rows if (abs(r["log2foldchange"]) >= args.l2fc_threshold and \
                                          r["padj"] <= args.padj_threshold) or \
                                         max(r["mean_base"], r["mean_testgroup"]) >= args.mean_threshold]

for row in filtered_rows:
    sequence_id = row["sequence_id"]
    print(sequence_id)
