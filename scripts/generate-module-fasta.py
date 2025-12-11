import os
import argparse
import subprocess
from defaults import DefaultPath
from needle.ortholog import load_modules
from needle.hmm import parse_hmmsearch_domtbl

parser = argparse.ArgumentParser(
    description="Generate FASTA file for KEGG module"
)
parser.add_argument("module_id", help="Module ID")
args = parser.parse_args()

print("Loading KEGG module and definition")
modules = load_modules(
  'data/modules.tsv',
  'data/ko.tsv',
  'data/module_ko.tsv',
  'data/ko.fasta.gz',
  module_id=args.module_id.upper()
)

assert len(modules.keys()) == 1
module = modules[args.module_id.upper()]

module_prefix = f"data/{args.module_id.lower()}"
module_ko_faa_fn = f"{module_prefix}_ko.faa"
domtbl_path = f"{module_prefix}_pfam.domtbl"
pfam_id_file = f"{module_prefix}_pfam_ids.txt"
module_hmm_fn = f"{module_prefix}_query.hmm"
module_stat_fn = f"{module_prefix}_query.hmm"
module_faa_fn = f"{module_prefix}_query.faa"

print("Curating KO consensus sequences")
module.create_consensus_aa_fasta(module_ko_faa_fn)

print("Scanning against Pfam")
cmd = ["hmmscan", "--cut_ga", "--domtblout", domtbl_path, DefaultPath.pfam_hmm(), module_ko_faa_fn]
subprocess.run(cmd, check=True, capture_output=True)

hmmscan_rows = parse_hmmsearch_domtbl(domtbl_path)
pfam_accessions = [row["target_accession"] for row in hmmscan_rows]
pfam_accessions = list(set(pfam_accessions))

with open(pfam_id_file, "w") as f:
    for acc in pfam_accessions:
        f.write(acc+"\n")

print("Fetch identified Pfam HMMs")
cmd = ["hmmfetch", "-f", "-o", module_hmm_fn, DefaultPath.pfam_hmm(), pfam_id_file]
subprocess.run(cmd, check=True, capture_output=True)

print("Generating final Pfam consensus query FASTA")

cmd = ["hmmemit", "-c", "-o", module_faa_fn, module_hmm_fn]
subprocess.run(cmd, check=True, capture_output=True)

cmd = ["hmmstat", module_hmm_fn]
res = subprocess.run(cmd, check=True, capture_output=True)
name_to_id = {}

for line in str(res.stdout).split("\\n"):
    if line.startswith("#"):
        continue
    tokens = line.split()
    if len(tokens) >= 3:
        name_to_id[tokens[1]] = tokens[2]

updated_lines = []
with open(module_faa_fn, "r") as f:
    lines = f.readlines()

for line in lines:
    line = line.strip()
    if line.startswith(">"):
        assert line.endswith("-consensus")
        name = line[1:].replace("-consensus", "")
        updated_lines.append(">"+name_to_id[name])
    else:
        updated_lines.append(line)

with open(module_faa_fn, "w") as f:
    for line in updated_lines:
        f.write(line+"\n")

print("Cleanup")
os.remove(module_ko_faa_fn)
os.remove(domtbl_path)
os.remove(pfam_id_file)
