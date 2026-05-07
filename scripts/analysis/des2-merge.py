# curates tall table of des2 data, with a manifest

import re
import csv
import argparse
from pathlib import Path
from tangle.sequence import read_fasta_as_dict
from tangle.manifest import ManifestTable

parser = argparse.ArgumentParser()
parser.add_argument("experiment_id")
parser.add_argument("output_dir")
parser.add_argument("faa_file")
parser.add_argument("des2_tsv_files", nargs="+")
args = parser.parse_args()

ANALYSIS_FIELD = "analysis_type"
seq_dict = read_fasta_as_dict(args.faa_file)
# strip off orifypy suffix on top of trinity suffix
seq_ids = { re.sub(r'_i\d+_ORF\.\d+$', '', k) for k in seq_dict.keys() }

entries = []

for tsv_file in args.des2_tsv_files:
    analysis_type = Path(tsv_file).stem
    print(analysis_type)
    with open(tsv_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row[ANALYSIS_FIELD] = analysis_type
            row["experiment_id"] = args.experiment_id
            entries.append(row)

manifest_tsv = args.output_dir+"/sequence_list.tsv"
des2_merged_tsv = args.output_dir+"/des2_tall.tsv"

sequences = []
for x in seq_ids:
    sequences.append(dict(sequence_accession=x, sequence_database=args.experiment_id, sequence_type="transcript", sequence_source="paper"))
ManifestTable.write_tsv(manifest_tsv, sequences)

with open(des2_merged_tsv, "w") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(entries[0].keys()))
    writer.writeheader()
    for entry in entries:
        writer.writerow(entry)
