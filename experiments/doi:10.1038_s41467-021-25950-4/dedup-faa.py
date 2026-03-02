import re
import argparse
from collections import defaultdict
from needle.seq import read_fasta_as_dict, write_fasta_from_dict

parser = argparse.ArgumentParser()
parser.add_argument("faa_file")
args = parser.parse_args()

seq_dict = read_fasta_as_dict(args.faa_file)

by_gene = defaultdict(list)

for k, v in seq_dict.items():
    gene_id = re.sub(r'_i\d+_ORF\.\d+$', '', k)

    found = False
    for saved_sequence_id, saved_aa in by_gene[gene_id]:
        if v in saved_aa:
            found = True
            # print(f"{k} contained in {saved_sequence_id}, ignore")
        elif saved_aa in v:
            found = True
            # print(f"{saved_sequence_id} contained in {k}, swap")
            by_gene[gene_id] = [t for t in by_gene[gene_id] if t[0] != saved_sequence_id]
            by_gene[gene_id].append((k, v))

    if not found:
        by_gene[gene_id].append((k, v))

entries = {}
for gene, group in by_gene.items():
    for k, v in group:
        entries[k] = v
print(len(entries))

write_fasta_from_dict(entries, args.faa_file+"_dedup")
