import itertools
from needle.seq import read_fasta_as_dict, write_fasta_from_dict

sequences = read_fasta_as_dict("ssa01.faa")

prefix_len = len("GLGU01070895.1")

kf = lambda k: k[:prefix_len]

selected = {}
for acc, seq in sequences.items():
    transcript = kf(acc)

    if transcript not in selected:
        selected[transcript] = seq
    elif len(seq) > len(selected[transcript]):
        selected[transcript] = seq

write_fasta_from_dict(selected, "ssa01.selected.faa")

