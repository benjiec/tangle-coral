# writes old detected format to DetectedTable format

import argparse
import pandas as pd
from pathlib import Path
from tangle.exp import TranscriptCountsTable, DESeq2Table
from tangle import unique_batch

ap = argparse.ArgumentParser()
ap.add_argument("experiment_dir")
ap.add_argument("experiment_id")
args = ap.parse_args()

sequence_data_tsv = Path(args.experiment_dir) / "sequence_data.tsv"
deseq2_tsv = Path(args.experiment_dir) / "des2_tall.tsv"

df = pd.read_csv(sequence_data_tsv, sep='\t')
df["experiment_id"] = args.experiment_id
rows = df.to_dict(orient='records')
TranscriptCountsTable.write_tsv(str(sequence_data_tsv)+".new", rows)

df = pd.read_csv(deseq2_tsv, sep='\t')
df["experiment_id"] = args.experiment_id
rows = df.to_dict(orient='records')
DESeq2Table.write_tsv(str(deseq2_tsv)+".new", rows)
