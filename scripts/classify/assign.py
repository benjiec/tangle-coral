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
proteins_faa = f"data/{args.module_id}_results/proteins.faa"
output_dir = f"data/{args.module_id}_results/faa"
domain_hmm_db_name = Path(args.pfam_hmm_file).stem

proteins_faa = [proteins_faa]
if args.additional_genome_accession:
    proteins_faa.append(DefaultPath.ncbi_genome_protein_faa(args.additional_genome_accession))

os.makedirs(output_dir, exist_ok=True)
assign_ko(classify_rows, ortholog_hmm_db_name, proteins_faa, output_dir, domain_hmm_db_name = domain_hmm_db_name)
