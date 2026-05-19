import re
import argparse
from collections import defaultdict
from tangle.detected import DetectedTable
from tangle.models import CSVSource

ap = argparse.ArgumentParser()
ap.add_argument("sequence_detected_tsv")
args = ap.parse_args()

source = CSVSource(DetectedTable, args.sequence_detected_tsv)
rows = source.values()

by_gene_and_target = {}

for row in rows:
    protein_acc = row["query_accession"]
    target_acc = row["target_accession"]
    bitscore = row["bitscore"]
    gene_acc = re.sub(r"_i(\d+)_ORF\.(\d+)$", "", protein_acc)
    assert gene_acc != protein_acc and protein_acc.startswith(gene_acc)
    k = (gene_acc, target_acc)
    if k not in by_gene_and_target or by_gene_and_target[k]["bitscore"] < bitscore:
        by_gene_and_target[k] = row

by_gene = defaultdict(list)

for (gene_acc, _), row in by_gene_and_target.items():
    by_gene[gene_acc].append(row)

rows = []
for gene_acc, gene_rows in by_gene.items():
    gene_rows = sorted(gene_rows, key=lambda r: float(r["evalue"]))
    for i, gene_row in enumerate(gene_rows):
        gene_row["query_accession"] = gene_acc
        gene_row["custom_metric_name"] = "evalue-rank"
        gene_row["custom_metric_value"] = i+1
    rows.extend(gene_rows)

fn = args.sequence_detected_tsv+"_aggregated"
DetectedTable.write_tsv(fn, rows)
