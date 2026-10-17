import sys
import json
import argparse
from tangle.sequence import read_fasta_as_dict
from tangle import open_file_to_write

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_x")
    parser.add_argument("fasta_y")
    parser.add_argument("json_output")
    parser.add_argument("--x", nargs="*", type=str)
    args = parser.parse_args(argv)

    seqs_x = read_fasta_as_dict(args.fasta_x)
    seqs_y = read_fasta_as_dict(args.fasta_y)
    if args.x:
        seqs_x = {k:v for k,v in seqs_x.items() if k in args.x}

    jobs = []

    for x_name, x_sequence in seqs_x.items():
        for y_name, y_sequence in seqs_y.items():
            job = dict(name=f"{x_name} {y_name}",
                       modelSeeds=[],
                       sequences=[dict(proteinChain=dict(sequence=x_sequence, count=1)),
                                  dict(proteinChain=dict(sequence=y_sequence, count=1))],
                       dialect="alphafoldserver",
                       version=1)
            jobs.append(job)

    with open_file_to_write(args.json_output, "wt") as f:
        f.write(json.dumps(jobs))

if __name__ == "__main__":
    sys.exit(main())
