import re
import argparse
from collections import defaultdict
from tangle.exp import DESeq2Table, TranscriptGenesTable
from tangle.detected import DetectedTable
from tangle.models import CSVSource

ap = argparse.ArgumentParser()
ap.add_argument("des2_tall_fn")
ap.add_argument("--transcript-genes-fn")
ap.add_argument("--transcript-proteins-fn")
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
des2_rows = [row for row in des2_rows if row["sequence_type"] == "gene"]

prev_sequences = {r["query_accession"]:1 for r in prev_results}
filtered_rows = [r for r in des2_rows if abs(r["log2foldchange"]) >= args.l2fc_threshold and \
                                         r["padj"] <= args.padj_threshold]
gene_ids = {r["sequence_id"]:1 for r in filtered_rows}
final_ids = gene_ids

if args.transcript_genes_fn:
    source = CSVSource(TranscriptGenesTable, args.transcript_genes_fn)
    transcript_genes = source.values()
    transcript_ids = {}
    for row in transcript_genes:
        if row["gene_id"] in gene_ids:
            transcript_ids[row["transcript_id"]] = 1
    final_ids = transcript_ids

    if args.transcript_proteins_fn:
        source = CSVSource(DetectedTable, args.transcript_proteins_fn)
        transcript_proteins = source.values()

        cds_ids = {}
        for row in transcript_proteins:
            if row["query_accession"] in transcript_ids:
                cds_ids[row["target_accession"]] = 1

        protein_ids = {}
        for row in transcript_proteins:
            if row["query_accession"] in cds_ids:
                protein_ids[row["target_accession"]] = 1
        final_ids = protein_ids

for s in final_ids.keys():
    print(s)
