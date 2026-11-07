import re
import sys
import json
import argparse
from tangle.sequence import read_fasta_as_dict
from tangle import open_file_to_write, open_file_to_read

# input file should be
#   candidate_faa - fasta file of all sequences
#   selection_fn - each line is comma separated list of targets

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_faa")
    parser.add_argument("selection_fn")
    parser.add_argument("json_output")
    args = parser.parse_args(argv)

    candidates = read_fasta_as_dict(args.candidate_faa, preserve_full_accession=True)
    selections = []
    error = False

    with open_file_to_read(args.selection_fn) as f:
        for line in f:
            seq_names = line.strip().split(",")
            for sn in seq_names:
                if sn not in candidates:
                    print(f"{sn} not in candidate fasta file")
                    error = True
            selections.append(seq_names)

    jobs = []
    if error is False:
        for seq_names in selections:
            sequences=[dict(proteinChain=dict(sequence=candidates[sn], count=1)) for sn in seq_names]
            names = " ".join([re.sub(r"\W","_",sn) for sn in seq_names])
            job = dict(name=names,
                       modelSeeds=[],
                       sequences=sequences,
                       dialect="alphafoldserver",
                       version=1)
            jobs.append(job)

        with open_file_to_write(args.json_output, "wt") as f:
            f.write(json.dumps(jobs))

if __name__ == "__main__":
    sys.exit(main())
