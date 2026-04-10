import os
import csv
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("data_dir")
parser.add_argument("output_dir")
args = parser.parse_args()


def get_entries(metadata, genome_accession, quant_fn):
    entries = []

    with open(quant_fn, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            d = metadata.copy()
            d["genome_accession"] = genome_accession
            d["sequence_id"] = row["Name"]
            d["count"] = row["NumReads"]
            d["tpm"] = row["TPM"]

            # convert to protein name
            if d["sequence_id"].startswith("lcl|"):
                d["sequence_id"] = "_".join(d["sequence_id"].split("_cds_")[1].split("_")[:-1])

            entries.append(d)

    return entries


def process_dir(rootdirn, genome_accession, fn_to_metadata):
    root = Path(rootdirn)
    subdirs = [x for x in root.iterdir() if x.is_dir()]
    entries = []

    for subdir in subdirs:
        quant_file = subdir / "quant.sf"
        if os.path.exists(quant_file):
            print(subdir.stem, quant_file)
            entries.extend(get_entries(fn_to_metadata[subdir.stem], genome_accession, quant_file))

    return entries


md = {
  "SRR9331965": dict(timepoint=1, cohort="WT1", sample="27_50.9"),
  "SRR9331964": dict(timepoint=1, cohort="WT1", sample="27_50.8"),
  "SRR9331963": dict(timepoint=1, cohort="WT1", sample="27_50.7"),
  "SRR9331962": dict(timepoint=1, cohort="WT1", sample="27_50.6"),
  "SRR9331961": dict(timepoint=1, cohort="WT1", sample="27_50.5"),
  "SRR9331960": dict(timepoint=1, cohort="WT1", sample="27_50.1"),
  "SRR9331959": dict(timepoint=1, cohort="SS8", sample="27_40.9"),
  "SRR9331958": dict(timepoint=1, cohort="SS8", sample="27_40.8"),
  "SRR9331957": dict(timepoint=1, cohort="SS8", sample="27_40.7"),
  "SRR9331956": dict(timepoint=1, cohort="SS8", sample="27_40.6"),
  "SRR9331955": dict(timepoint=1, cohort="SS8", sample="27_40.5"),
  "SRR9331954": dict(timepoint=1, cohort="SS5", sample="27_25.9"),
  "SRR9331953": dict(timepoint=1, cohort="SS5", sample="27_25.8"),
  "SRR9331952": dict(timepoint=1, cohort="SS5", sample="27_25.7"),
  "SRR9331951": dict(timepoint=1, cohort="SS5", sample="27_25.6"),
  "SRR9331950": dict(timepoint=1, cohort="SS5", sample="27_25.5"),
  "SRR9331949": dict(timepoint=1, cohort="SS5", sample="27_25.1"),
  "SRR9331948": dict(timepoint=1, cohort="SS3", sample="27_14.8"),
  "SRR9331947": dict(timepoint=1, cohort="SS3", sample="27_14.7"),
  "SRR9331946": dict(timepoint=1, cohort="SS3", sample="27_14.6"),
  "SRR9331945": dict(timepoint=1, cohort="SS3", sample="27_14.5"),
  "SRR9331944": dict(timepoint=1, cohort="SS3", sample="27_14.4"),
  "SRR9331943": dict(timepoint=1, cohort="SS3", sample="27_14.3"),
}

all_entries = []
all_entries.extend(process_dir(args.data_dir+"/aten-quants", "doi:10.1126/sciadv.aba2498-a_tenuis", md))
all_entries.extend(process_dir(args.data_dir+"/c_goreaui-quants", "doi:10.1126/sciadv.aba2498-c_goreaui", md))

with open(args.output_dir+"/sequence_data.tsv", "w") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(all_entries[0].keys()))
    writer.writeheader()
    for entry in all_entries:
        writer.writerow(entry)
