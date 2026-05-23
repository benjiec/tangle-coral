# curates transcript manifest

import argparse
from tangle.sequence import read_fasta_as_dict
from tangle.manifest import ManifestTable

parser = argparse.ArgumentParser()
parser.add_argument("experiment_id")
parser.add_argument("output_dir")
parser.add_argument("transcripts_fna_file")
args = parser.parse_args()

if Path(args.transcripts_fna_file+".gz").exists():
    seq_dict = read_fasta_as_dict(args.transcripts_fna_file+".gz")
else:
    seq_dict = read_fasta_as_dict(args.transcripts_fna_file)

manifest_tsv = args.output_dir+"/transcript_list.tsv"

sequences = []
for k,v in seq_dict.items():
    sequences.append(dict(
        sequence_accession=k,
        sequence_database=args.experiment_id,
        sequence_type="transcript",
        sequence_source="paper",
        sequence_length=len(v))
    )
ManifestTable.write_tsv(manifest_tsv, sequences)
