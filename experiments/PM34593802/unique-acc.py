import argparse
from pathlib import Path
from needle.seq import read_fasta_as_dict, write_fasta_from_dict

parser = argparse.ArgumentParser()
parser.add_argument("fna")
args = parser.parse_args()

sequences = read_fasta_as_dict(args.fna)
prefix = Path(args.fna).stem.split(".")[0]

updated = {f"{prefix}_{k}":v for k,v in sequences.items()}
write_fasta_from_dict(updated, args.fna)

print(len(sequences), prefix)
