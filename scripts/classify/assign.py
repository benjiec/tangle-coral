import os
from pathlib import Path
from scripts.defaults import DefaultPath
from needle.classify import assign_ko, ClassifyTSV

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("ortholog_hmm_file")
ap.add_argument("pfam_hmm_file")
ap.add_argument("module_id")
ap.add_argument("--additional-genome-accession", type=str, default=None)
args = ap.parse_args()

classify_tsv = f"data/{args.module_id}_results/classify.tsv"

classify_rows = ClassifyTSV.from_tsv_to_rows(classify_tsv)
ortholog_hmm_db_name = Path(args.ortholog_hmm_file).stem
protein_faa = f"data/{args.module_id}_results/proteins.faa"
output_dir = f"data/{args.module_id}_results/faa"
domain_hmm_db_name = Path(args.pfam_hmm_file).stem

os.makedirs(output_dir, exist_ok=True)

#
# first, non-reference genome, use default scoring threshold
#
assign_ko(classify_rows, ortholog_hmm_db_name, protein_faa, output_dir, prefix_match = True, domain_hmm_db_name = domain_hmm_db_name, score_to_threshold_ratio=0.9)

#
# then, reference genome, use more stringent scoring threshold
#
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

# use more stringent scoring threshold
assign_ko(classify_rows, ortholog_hmm_db_name, multiple_protein_faas, output_dir, domain_hmm_db_name = domain_hmm_db_name, score_to_threshold_ratio=1)
