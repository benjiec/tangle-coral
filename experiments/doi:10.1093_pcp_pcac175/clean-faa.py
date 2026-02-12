import csv
from needle.seq import read_fasta_as_dict, write_fasta_from_dict

data_fn = "alldata.csv"
protein_faa = "uniprotkb_taxonomy_id_2951_AND_reviewed_2026_02_12.fasta"
to_faa = "proteins.faa"

protein_sequences = read_fasta_as_dict(protein_faa)
new_dict = {}

for acc,seq in protein_sequences.items():
    new_acc = acc.split("|")[1]
    new_dict[new_acc] = seq

write_fasta_from_dict(new_dict, to_faa)
