# writes old detected format to DetectedTable format

import argparse
import pandas as pd
from tangle.detected import DetectedTable
from tangle import unique_batch

ap = argparse.ArgumentParser()
ap.add_argument("tsv_file")
ap.add_argument("experiment_id")
ap.add_argument("target_database")
ap.add_argument("--output", "-o", help="Output file path (defaults to input)")
ap.add_argument("--pfam-mode", action="store_true", default=False)
args = ap.parse_args()

df = pd.read_csv(args.tsv_file, sep='\t')
output_path = args.output if args.output else args.tsv_file

df["detection_type"] = "sequence"
df["detection_method"] = "hmm"
df["batch"] = unique_batch()
df["query_database"] = args.experiment_id
df["query_accession"] = df["protein_accession"]
df["query_type"] = "protein"
df["query_start"] = df["protein_start"]
df["query_end"] = df["protein_end"]
df["target_database"] = args.target_database
df["target_accession"] = df["hmm_accession"]
df["target_type"] = "protein"
df["target_start"] = df["hmm_start"]
df["target_end"] = df["hmm_end"]
df["evalue"] = df["dom_evalue"]
df["bitscore"] = df["dom_score"]

if not args.pfam_mode:
  df["bitscore_threshold"] = df["score_threshold"]
  df["custom_metric_name"] = "evalue-rank"
  df["custom_metric_value"] = df["dom_rank_for_protein"]

df = df.drop(columns=['hmm_db', 'protein_end', 'dom_evalue_cond',
                      'hmm_accession', 'hmm_start', 'score_threshold', 'hmm_end', 'dom_evalue',
                      'protein_accession', 'protein_start', 'dom_score', 'genome_accession',
                      'dom_rank_for_protein'])

rows = df.to_dict(orient='records')
DetectedTable.write_tsv(output_path, rows)
