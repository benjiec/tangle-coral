import re
import os
import csv
import argparse
import itertools
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("data_dir")
parser.add_argument("output_dir")
args = parser.parse_args()

from collections import defaultdict

def tximport(metadata, quant_fn):
    aggregated = defaultdict(lambda: {"count": 0.0, "tpm": 0.0, "iso_count": 0})

    with open(quant_fn, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seq_id = re.sub(r'_i\d+$', '', row["Name"])
            aggregated[seq_id]["count"] += float(row["NumReads"])
            aggregated[seq_id]["tpm"] += float(row["TPM"])
            aggregated[seq_id]["iso_count"] += 1

    entries = []
    for seq_id, stats in aggregated.items():
        d = metadata.copy()
        d["sequence_id"] = seq_id
        d["count"] = stats["count"]
        # Standard tximport sums TPMs for the gene level
        d["tpm"] = stats["tpm"] 
        # Optional: store how many isoforms were collapsed
        d["isoform_count"] = stats["iso_count"]
        entries.append(d)

    return entries


def process_dir(rootdirn, fn_to_metadata):
    root = Path(rootdirn)
    subdirs = [x for x in root.iterdir() if x.is_dir()]
    entries = []

    for subdir in subdirs:
        quant_file = subdir / "quant.sf"
        if os.path.exists(quant_file):
            print(subdir.stem, quant_file)
            entries.extend(tximport(fn_to_metadata[subdir.stem], quant_file))

    return entries


md = {
  "SRR6255820": "Orbicella faveolata, Symbiodinium A3, control, A",
  "SRR6255822": "Orbicella faveolata, Symbiodinium A3, control, B",
  "SRR6255825": "Orbicella faveolata, Symbiodinium A3, control, C",
  "SRR6255844": "Orbicella faveolata, Symbiodinium A3, treatment, A",
  "SRR6255862": "Orbicella faveolata, Symbiodinium A3, treatment, B",
  "SRR6255864": "Orbicella faveolata, Symbiodinium A3, treatment, C",
  "SRR6255880": "Pseudodiploria clivosa, Breviolum faviinorum, control, A",
  "SRR6255882": "Pseudodiploria clivosa, Breviolum faviinorum, control, B",
  "SRR6255888": "Pseudodiploria clivosa, Breviolum faviinorum, control, C",
  "SRR6255887": "Pseudodiploria clivosa, Breviolum faviinorum, treatment, A",
  "SRR6255886": "Pseudodiploria clivosa, Breviolum faviinorum, treatment, B",
  "SRR6255890": "Pseudodiploria clivosa, Breviolum faviinorum, treatment, C",
  "SRR6255889": "Siderastrea radians, Breviolum B5, control, A",
  "SRR6256312": "Siderastrea radians, Breviolum B5, control, B",
  "SRR6256323": "Siderastrea radians, Breviolum B5, control, C",
  "SRR6256322": "Siderastrea radians, Breviolum B5, treatment, A",
  "SRR6256324": "Siderastrea radians, Breviolum B5, treatment, B",
  "SRR6256325": "Siderastrea radians, Breviolum B5, treatment, C",
}

symb_md = {
  fn: dict(
    timepoint=0,
    cohort=desc.split(", ")[2],
    sample=desc.split(", ")[3],
    genome_accession="doi:10.1038_s41467-021-25950-4_"+desc.split(", ")[1].lower().replace(" ", "_")
  )
  for fn,desc in md.items()
}

host_md = {
  fn: dict(
    timepoint=0,
    cohort=desc.split(", ")[2],
    sample=desc.split(", ")[3],
    genome_accession="doi:10.1038_s41467-021-25950-4_"+desc.split(", ")[0].lower().replace(" ", "_")
  )
  for fn,desc in md.items()
}

host_entries = process_dir(args.data_dir+"/host_quants", host_md)
symb_entries = process_dir(args.data_dir+"/symb_quants", symb_md)
all_entries = host_entries + symb_entries

kf = lambda entry: entry["genome_accession"]
all_entries = sorted(all_entries, key=kf)
for acc, group in itertools.groupby(all_entries, kf):
    group = list(group)
    print(acc, len(group))
    avg_count = sum([float(x["count"]) for x in group]) / len(group)
    print(acc, avg_count)

with open(args.output_dir+"/sequence_data.tsv", "w") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(all_entries[0].keys()))
    writer.writeheader()
    for entry in all_entries:
        writer.writerow(entry)
