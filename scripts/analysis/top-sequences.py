import re
import argparse
from collections import defaultdict
from tangle.exp import DESeq2Table
from tangle.detected import DetectedTable
from tangle.models import CSVSource

ap = argparse.ArgumentParser()
ap.add_argument("des2_tall_fn")
ap.add_argument("--results-fn", help="If specified, ignore sequences already ran, from this file")
ap.add_argument("--l2fc-threshold", default=1, type=float)
ap.add_argument("--padj-threshold", default=0.05, type=float)
args = ap.parse_args()

prev_results = []
if args.results_fn:
    source = CSVSource(DetectedTable, args.results_fn)
    prev_results = source.values()

source = CSVSource(DESeq2Table, args.des2_tall_fn)
des2_rows = source.values()

prev_sequences = {r["query_accession"]:1 for r in prev_results}
filtered_rows = [r for r in des2_rows if abs(r["log2foldchange"]) >= args.l2fc_threshold and \
                                         r["padj"] <= args.padj_threshold]
sequence_ids = {r["sequence_id"]:1 for r in filtered_rows}
for s in sequence_ids.keys():
    print(s)
