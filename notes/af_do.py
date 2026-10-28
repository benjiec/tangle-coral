import sys
import json
import argparse
from tangle.sequence import read_fasta_as_dict
from tangle import open_file_to_write

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("json_output")
    args = parser.parse_args(argv)

    seqs = read_fasta_as_dict(args.fasta)
    jobs = []

    for x_name, x_sequence in seqs.items():
        job = dict(name=f"{x_name}",
                   modelSeeds=[],
                   sequences=[dict(proteinChain=dict(sequence=x_sequence, count=1))],
                   dialect="alphafoldserver",
                   version=1)
        jobs.append(job)

    with open_file_to_write(args.json_output, "wt") as f:
        f.write(json.dumps(jobs))

if __name__ == "__main__":
    sys.exit(main())
