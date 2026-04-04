import re
import csv
import argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("sequence_ko_tsv")
args = ap.parse_args()

with open(args.sequence_ko_tsv, "r") as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = [row for row in reader]

by_gene_and_hmm = {}

for row in rows:
    protein_acc = row["protein_accession"]
    hmm_acc = row["hmm_accession"]
    dom_score = row["dom_score"]
    gene_acc = re.sub(r"_i(\d+)_ORF\.(\d+)$", "", protein_acc)
    assert gene_acc != protein_acc and protein_acc.startswith(gene_acc)
    k = (gene_acc, hmm_acc)
    if k not in by_gene_and_hmm or by_gene_and_hmm[k]["dom_score"] < dom_score:
        by_gene_and_hmm[k] = row

by_gene = defaultdict(list)

for (gene_acc, _), row in by_gene_and_hmm.items():
    by_gene[gene_acc].append(row)

rows = []
for gene_acc, gene_rows in by_gene.items():
    gene_rows = sorted(gene_rows, key=lambda r: float(r["dom_evalue"]))
    for i, gene_row in enumerate(gene_rows):
        gene_row["protein_accession"] = gene_acc
        gene_row["dom_rank_for_protein"] = i+1
    rows.extend(gene_rows)

with open(args.sequence_ko_tsv+"_aggregated", "w") as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=list(rows[0].keys()))
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
