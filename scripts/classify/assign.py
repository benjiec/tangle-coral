import os
import csv
import argparse
from pathlib import Path
from scripts.defaults import DefaultPath
from needle.duckdb import load, CandidateClassifiedProteins

ap = argparse.ArgumentParser()
ap.add_argument("module_id")
ap.add_argument("--additional-genome-accession", type=str, default=None)
args = ap.parse_args()

protein_faa = f"data/{args.module_id}_results/proteins.faa"
output_dir = f"data/{args.module_id}_results/faa"
os.makedirs(output_dir, exist_ok=True)

multiple_protein_faas = []
if args.additional_genome_accession:
    if os.path.exists(args.additional_genome_accession):
        with open(args.additional_genome_accession, "r") as f:
            for acc in f.readlines():
                acc = acc.strip()
                multiple_protein_faas.append(DefaultPath.ncbi_genome_protein_faa(acc))
                print("additional protein sequences from", multiple_protein_faas[-1])
    else:
        multiple_protein_faas.append(DefaultPath.ncbi_genome_protein_faa(args.additional_genome_accession))
        print("additional protein sequences from", multiple_protein_faas[-1])

load(args.module_id)
candidate_proteins = CandidateClassifiedProteins()
ko_matches = candidate_proteins.ko_matches()
pfam_matches = candidate_proteins.pfam_matches()

print(ko_matches)
print(pfam_matches)

sorter = lambda d: (d['genome_accession'], d['protein_accession'], d['hmm_db'], d['dom_rank_for_protein'], d['hmm_accession'], d['hmm_start'], d['hmm_end'])

ko_matches = ko_matches.to_dict(orient='records')
ko_matches = sorted(ko_matches, key=sorter)
print(ko_matches[0])

pfam_matches = pfam_matches.to_dict(orient='records')
pfam_matches = sorted(pfam_matches, key=sorter)
print(pfam_matches[0])

fieldnames = list(ko_matches[0].keys())
output_ko = f"data/{args.module_id}_results/candidate_ko.tsv"
output_pf = f"data/{args.module_id}_results/candidate_pfam.tsv"

with open(output_ko, "w") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for rec in ko_matches:
        writer.writerow(rec)

with open(output_pf, "w") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for rec in pfam_matches:
        writer.writerow(rec)




"""
#
# first, non-reference genome, use default scoring threshold
#
assign_ko(classify_rows, ko_hmm_db_name, protein_faa, output_dir, score_to_threshold_ratio=0.9)

#
# then, reference genome, use more stringent scoring threshold
#

# use more stringent scoring threshold
assign_ko(classify_rows, ko_hmm_db_name, multiple_protein_faas, output_dir, score_to_threshold_ratio=1)
"""
