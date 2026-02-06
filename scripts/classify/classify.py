import os
import csv
import gzip
import tempfile
from pathlib import Path
from scripts.defaults import DefaultPath
from needle.seq import read_fasta_as_dict, write_fasta_from_dict
from needle.match import ProteinsTSV
from needle.classify import ClassifyTSV, classify

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("hmm_file")
ap.add_argument("module_id")
ap.add_argument("output_tsv")
ap.add_argument("--disable-cutoff-ga", action="store_true", default=False)
ap.add_argument("--genome-accession", type=str, default=None)
ap.add_argument("--fasta-file", type=str, default=None)
ap.add_argument("--filter-by", type=str, default=None)
ap.add_argument("--requires-prefix-match", action="store_true", default=False)
ap.add_argument("--cpu", type=int, default=None)
ap.add_argument("--max-rank", type=int, default=25)
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

proteins_faa = None
protein_genome_accession_dict = None

output_tsv = args.output_tsv

if args.fasta_file:
    if args.genome_accession is None:
        print("Classifying a protein FASTA file requires a genome accession")
        exit(-1)
    if not os.path.exists(args.fasta_file):
        print(f"Cannot find protein FASTA file {args.fasta_file}")
        exit(-1)
    proteins_faa = args.fasta_file
    proteins_fasta = read_fasta_as_dict(proteins_faa)
    protein_genome_accession_dict = { k: args.genome_accession for k in proteins_fasta.keys() }

elif args.genome_accession:
    proteins_faa = DefaultPath.ncbi_genome_protein_faa(args.genome_accession)
    if not os.path.exists(proteins_faa):
        print(f"Cannot find protein FASTA for genome {args.genome_accession}")
        exit(-1)
    proteins_fasta = read_fasta_as_dict(proteins_faa)
    protein_genome_accession_dict = { k: args.genome_accession for k in proteins_fasta.keys() }

else:    
    proteins_faa = f"data/{args.module_id}_results/protein_detected.faa"
    proteins_tsv = f"data/{args.module_id}_results/protein_detected.tsv"
    if not os.path.exists(proteins_faa) or not os.path.exists(proteins_tsv):
        print(f"Cannot find module assets in data/{args.module_id}_results")
        exit(-1)
    protein_tsv_rows = ProteinsTSV.from_tsv_to_rows(proteins_tsv)
    protein_genome_accession_dict = {
      row["protein_hit_id"]: row["genome_accession"]
      for row in protein_tsv_rows
    }

assert proteins_faa and protein_genome_accession_dict
print(f"classifying {proteins_faa}")

tmp_fn = None
if args.filter_by and os.path.exists(args.filter_by):
    print(f"filtering fasta using {args.filter_by}")

    proteins_fasta = read_fasta_as_dict(proteins_faa)
    classify_rows = ClassifyTSV.from_tsv_to_rows(args.filter_by)
    proteins = {row["protein_accession"] for row in classify_rows}
    proteins_fasta = {k:v for k,v in proteins_fasta.items() if k in proteins}

    tmpf = tempfile.NamedTemporaryFile(delete=False, suffix=".faa", mode="w")
    tmpf.close()
    write_fasta_from_dict(proteins_fasta, tmpf.name)
    proteins_faa = tmpf.name
    tmp_fn = tmpf.name

classify(args.hmm_file, proteins_faa, cutoff_ga, output_tsv, protein_genome_accession_dict, score_threshold_dict,
         requires_prefix_match = args.requires_prefix_match, cpu = args.cpu, max_rank = args.max_rank)

if tmp_fn:
    os.remove(tmp_fn)
