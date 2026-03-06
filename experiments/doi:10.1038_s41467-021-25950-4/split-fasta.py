import os
import math
import random
import argparse
from needle.seq import read_fasta_as_dict, write_fasta_from_dict

parser = argparse.ArgumentParser()
parser.add_argument("fna")
parser.add_argument("n", type=int)
args = parser.parse_args()

sequences = read_fasta_as_dict(args.fna)

accs = [x for x in sequences.keys()]
random.shuffle(accs)

pre, ext = os.path.splitext(args.fna)
step = math.ceil(len(accs) / int(args.n))

j = 0
for i in range(0, len(accs), step):
    print(i)
    selected = {k:sequences[k] for k in accs[i:i+step]}
    write_fasta_from_dict(selected, f"{pre}_{j}{ext}")
    j += 1
