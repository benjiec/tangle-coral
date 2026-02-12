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
  "SRR33692232": dict(cohort="34C", timepoint=0, sample="1"),
  "SRR36399503": dict(cohort="34C", timepoint=0, sample="2"),
  "SRR36399537": dict(cohort="34C", timepoint=0, sample="3"),
  "SRR36399877": dict(cohort="34C", timepoint=0, sample="4"),
  "SRR36401143": dict(cohort="34C", timepoint=3, sample="3"),
  "SRR36401144": dict(cohort="34C", timepoint=3, sample="2"),
  "SRR36401145": dict(cohort="34C", timepoint=3, sample="1"),
  "SRR36417296": dict(cohort="34C", timepoint=12, sample="3"),
  "SRR36417297": dict(cohort="34C", timepoint=12, sample="2"),
  "SRR36417298": dict(cohort="34C", timepoint=12, sample="1"),
  "SRR36417818": dict(cohort="34C", timepoint=24, sample="3"),
  "SRR36417819": dict(cohort="34C", timepoint=24, sample="2"),
  "SRR36417820": dict(cohort="34C", timepoint=24, sample="1"),
  "SRR36490025": dict(cohort="34C", timepoint=48, sample="3"),
  "SRR36490026": dict(cohort="34C", timepoint=48, sample="2"),
  "SRR36490027": dict(cohort="34C", timepoint=48, sample="1"),
  "SRR36494936": dict(cohort="34C", timepoint=96, sample="3"),
  "SRR36494937": dict(cohort="34C", timepoint=96, sample="2"),
  "SRR36494938": dict(cohort="34C", timepoint=96, sample="1"),
  "SRR36494939": dict(cohort="34C", timepoint=72, sample="3"),
  "SRR36494940": dict(cohort="34C", timepoint=72, sample="2"),
  "SRR36494941": dict(cohort="34C", timepoint=72, sample="1"),
  "SRR36495115": dict(cohort="27C", timepoint=192, sample="2"),
  "SRR36495116": dict(cohort="27C", timepoint=192, sample="1"),
  "SRR36495117": dict(cohort="34C", timepoint=192, sample="2"),
  "SRR36495118": dict(cohort="34C", timepoint=192, sample="1"),
}

all_entries = process_dir(args.data_dir+"/quants", "symbiodinium_linucheae_doi_10.1093_ismejo_wraf268_ssa01", md)

with open(args.output_dir+"/sequence_data.tsv", "w") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(all_entries[0].keys()))
    writer.writeheader()
    for entry in all_entries:
        writer.writerow(entry)
