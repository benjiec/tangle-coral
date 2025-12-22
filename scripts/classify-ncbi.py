import csv
import gzip
from defaults import DefaultPath
from needle.seq import read_fasta_as_dict
from needle.classify import classify

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("hmm_file")
ap.add_argument("genome_accession")
ap.add_argument("output_tsv")
ap.add_argument("--disable-cutoff-ga", action="store_true", default=False)
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

proteins_faa = DefaultPath.ncbi_genome_protein_faa(args.genome_accession)
proteins_fasta = read_fasta_as_dict(proteins_faa)
protein_genome_accession_dict = { k: args.genome_accession for k in proteins_fasta.keys() }

classify(args.hmm_file, proteins_faa, cutoff_ga, args.output_tsv, protein_genome_accession_dict, score_threshold_dict)
