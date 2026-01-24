import os
import argparse
from scripts.defaults import DefaultPath
from needle.duckdb import load, PutativeProteins
from needle.seq import read_fasta_as_dict, write_fasta_from_dict

parser = argparse.ArgumentParser()
parser.add_argument("module")
parser.add_argument("ko")
args = parser.parse_args()

load(args.module.lower())
protein_genome_accs = PutativeProteins(args.ko.upper()).protein_genome_accessions()

assigned_dict = read_fasta_as_dict(DefaultPath.module_ko_assigned_proteins(args.module, args.ko))
seq_dict = read_fasta_as_dict(DefaultPath.module_detected_proteins(args.module))

genome_accs = list(set([genome_acc for _, genome_acc in protein_genome_accs]))
for acc in genome_accs:
    fn = DefaultPath.ncbi_genome_protein_faa(acc)
    if os.path.exists(fn):
        d = read_fasta_as_dict(DefaultPath.ncbi_genome_protein_faa(acc))
        seq_dict.update(d)

found_seq = {}

for protein_acc, _ in protein_genome_accs:
    if protein_acc in assigned_dict:
        # print(f"skip assigned protein {protein_acc}")
        pass
    elif protein_acc not in seq_dict:
        print(f"cannot find {protein_acc} sequence")
    else:
        found_seq[protein_acc] = seq_dict[protein_acc]

output_fn = DefaultPath.module_ko_putative_proteins(args.module, args.ko)
write_fasta_from_dict(found_seq, output_fn)
