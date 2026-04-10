import re
import os
import csv
import argparse
import itertools
import pandas as pd
from pathlib import Path
from pytximport import tximport

parser = argparse.ArgumentParser()
parser.add_argument("data_dir")
parser.add_argument("output_dir")
args = parser.parse_args()

from collections import defaultdict

def load_sf(metadata, quant_fn):

    tx_to_gene = []
    with open(quant_fn, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = re.sub(r'_i\d+$', '', row["Name"])
            tx_to_gene.append(dict(transcript_id=row["Name"], gene_id=gene_id))
    tx_to_gene_map = pd.DataFrame(tx_to_gene)

    txi = tximport(
        [quant_fn],
        data_type="salmon",
        transcript_gene_map=tx_to_gene_map,
        counts_from_abundance="length_scaled_tpm"
    )

    counts_df = pd.DataFrame(
        txi.X.T,
        index=txi.var_names,   # These are your 'gene_id's (stripped Trinity names)
        columns=["corrected_counts"]
    )

    tpm_df = pd.DataFrame(
        txi.obsm['abundance'].T, # Transpose to get Genes as Rows
        index=txi.var_names,
        columns=[f"TPM" for s in txi.obs_names]
    )

    assert (counts_df.index == tpm_df.index).all(), "Index mismatch: Counts and TPM rows do not align"
   
    combined_df = pd.concat([counts_df, tpm_df], axis=1)
    combined_df.index.name = 'sequence_id'
    combined_df = combined_df.reset_index()
    for mk, mv in metadata.items():
        combined_df[mk] = mv
    combined_df["count"] = combined_df["corrected_counts"].astype(float).round().astype(int).round().astype(int)

    return combined_df


def process_dir(rootdirn, fn_to_metadata):
    root = Path(rootdirn)
    subdirs = [x for x in root.iterdir() if x.is_dir()]
    dataframes = []
    for subdir in subdirs:
        quant_file = subdir / "quant.sf"
        if os.path.exists(quant_file):
            print(subdir.stem, quant_file)
            dataframes.append(load_sf(fn_to_metadata[subdir.stem], quant_file))
    return dataframes


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

host_dataframes = process_dir(args.data_dir+"/host_quants", host_md)
symb_dataframes = process_dir(args.data_dir+"/symb_quants", symb_md)
all_dataframes = host_dataframes + symb_dataframes
tall_df = pd.concat(all_dataframes, axis=0, ignore_index=True)
tall_df.to_csv(args.output_dir+"/sequence_data.tsv", sep="\t", index=False)
