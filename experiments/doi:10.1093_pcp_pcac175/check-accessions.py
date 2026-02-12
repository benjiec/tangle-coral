import csv
from needle.seq import read_fasta_as_dict

data_fn = "alldata.csv"
protein_faa = "uniprotkb_taxonomy_id_2951_AND_reviewed_2026_02_12.fasta"

protein_sequences = read_fasta_as_dict(protein_faa)
fasta_accessions = {p.split("|")[1] for p in protein_sequences.keys()}

with open(data_fn, "r") as f:
    reader = csv.DictReader(f, delimiter=",")
    accessions = [row["accession"] for row in reader]

checked = {}
found = set()
for acc in accessions:
    if acc in checked:
        print(f"{acc} duplicated")
    checked[acc] = 1
    if acc not in fasta_accessions:
        print(f"{acc} not found")
    else:
        found.add(acc)
print(f"total {len(found)} found")
