import re
import csv
import argparse
import itertools
import pandas as pd
import numpy as np
from pathlib import Path
from pytximport import tximport
from tangle.exp import SequenceCountsTable, TranscriptGenesTable

parser = argparse.ArgumentParser()
parser.add_argument("experiment_id")
parser.add_argument("data_dir")
parser.add_argument("output_dir")
args = parser.parse_args()


def load_quants_by_gene(metadata: dict, quant_fn: str) -> pd.DataFrame:
    """
    Processes a single Salmon quant.sf file. Aggregates transcript metrics to
    gene level (Trinity convention). Returns a dataframe containing corrected
    counts, TPM, abundance-weighted effective length, and abundance-weighted
    physical length.
    """

    # Phase 1: Parse and calculate transcript-level weights
    tx_to_gene = []
    raw_records = []
    
    with open(quant_fn, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tx_id = row["Name"]
            if not tx_id.startswith(metadata["species_prefix"]):
                continue

            gene_id = re.sub(r'_i\d+$', '', tx_id)
            tx_to_gene.append({"transcript_id": tx_id, "gene_id": gene_id})
            raw_records.append({
                "gene_id": gene_id,
                "Length": float(row["Length"]),
                "TPM": float(row["TPM"])
            })
    
    tx_to_gene_map = pd.DataFrame(tx_to_gene)
    raw_df = pd.DataFrame(raw_records)

    # Phase 2: Compute abundance-weighted physical length
    raw_df['tpm_times_length'] = raw_df['TPM'] * raw_df['Length']
    grouped = raw_df.groupby('gene_id')
    
    tpm_sum = grouped['TPM'].sum()
    tpm_times_length_sum = grouped['tpm_times_length'].sum()
    
    # Avoid division by zero for unexpressed genes by falling back to arithmetic mean
    weighted_phys_len = np.where(
        tpm_sum > 0,
        tpm_times_length_sum / tpm_sum,
        grouped['Length'].mean()
    )
    
    phys_len_df = pd.DataFrame({
        f"weighted_physical_length": weighted_phys_len
    }, index=tpm_sum.index)

    # Phase 3: Execute pytximport for effective metrics
    txi = tximport(
        [quant_fn],
        data_type="salmon",
        transcript_gene_map=tx_to_gene_map,
        counts_from_abundance="length_scaled_tpm"
    )

    # Phase 4: Construct and align core dataframes (Genes to Rows via .T)
    counts_df = pd.DataFrame(
        txi.X.T, 
        index=txi.var_names, 
        columns=[f"count"]
    )
    
    tpm_df = pd.DataFrame(
        txi.obsm['abundance'].T, 
        index=txi.var_names, 
        columns=[f"tpm"]
    )
    
    eff_len_df = pd.DataFrame(
        txi.obsm['length'].T, 
        index=txi.var_names, 
        columns=[f"effective_length"]
    )
   
    # Map Phase 2 data explicitly to the pytximport variable index order
    phys_len_df = phys_len_df.reindex(txi.var_names)
    phys_len_df.columns = [f"weighted_physical_length"]

    # Phase 5: Consolidated Merge and Metadata Injection
    combined_df = pd.concat([counts_df, tpm_df, eff_len_df, phys_len_df], axis=1)
    
    combined_df.index.name = 'sequence_id'
    combined_df = combined_df.reset_index()
    combined_df['sequence_type'] = 'gene'
    
    for mk, mv in metadata.items():
        combined_df[mk] = mv
        
    combined_df["count"] = combined_df[f"count"].round().astype(int)

    return tx_to_gene_map, combined_df


def load_quants_by_transcript(metadata: dict, quant_fn: str) -> pd.DataFrame:
    """
    Processes a single Salmon quant.sf file at the transcript level.
    Bypasses tximport to retain individual isoform metrics (_i versions).
    """

    # Phase 1: Read raw Salmon quantification matrix directly
    raw_df = pd.read_csv(quant_fn, sep="\t")
    raw_df = raw_df[raw_df['Name'].str.startswith(metadata["species_prefix"])]

    # Phase 2: Extract, isolate, and rename target metrics
    # Salmon schema: Name | Length | EffectiveLength | TPM | NumReads
    transcript_df = pd.DataFrame(index=raw_df["Name"])

    transcript_df["count"] = raw_df["NumReads"].values
    transcript_df["tpm"] = raw_df["TPM"].values
    transcript_df["effective_length"] = raw_df["EffectiveLength"].values
    transcript_df["weighted_physical_length"] = raw_df["Length"].values

    # Phase 5: Index and Metadata Injection

    transcript_df.index.name = 'sequence_id'
    combined_df = transcript_df.reset_index()
    combined_df['sequence_type'] = 'transcript'

    for mk, mv in metadata.items():
        combined_df[mk] = mv

    combined_df["count"] = combined_df["count"].round().astype(int)

    return combined_df


def process_dir(rootdirn, fn_to_metadata):
    root = Path(rootdirn)

    tx_to_gene = []
    gene_dfs = []
    transcript_dfs = []

    for k,v in fn_to_metadata.items():
        quant_path = root / v["transcriptome"] / k / "quant.sf"
        print(quant_path)

        t2g, df = load_quants_by_gene(v, quant_path)
        tx_to_gene.append(t2g)
        gene_dfs.append(df)
        df = load_quants_by_transcript(v, quant_path)
        transcript_dfs.append(df)

    return tx_to_gene, gene_dfs, transcript_dfs


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
}

symb_md = {
  fn: dict(
    timepoint=0,
    cohort=desc["cohort"],
    sample=desc["sample"],
    transcriptome="a_tenuis_c_goreaui",
    species_prefix="c_goreaui",
    genome_accession="x_c_goreaui"
  )
  for fn, desc in md.items()
}

host_md = {
  fn: dict(
    timepoint=0,
    cohort=desc["cohort"],
    sample=desc["sample"],
    transcriptome="a_tenuis_c_goreaui",
    species_prefix="a_tenuis",
    genome_accession="x_a_tenuis"
  )
  for fn, desc in md.items()
}

host_tx_to_gene, host_gene_dfs, host_transcript_dfs = process_dir(args.data_dir, host_md)
symb_tx_to_gene, symb_gene_dfs, symb_transcript_dfs = process_dir(args.data_dir, symb_md)

all_tx_to_gene = host_tx_to_gene + symb_tx_to_gene
df = pd.concat(all_tx_to_gene, axis=0, ignore_index=True)
df = df.drop_duplicates()
rows = df.to_dict(orient="records")
for row in rows:
    row["experiment_id"] = args.experiment_id
TranscriptGenesTable.write_tsv(args.output_dir+"/transcript_genes.tsv", rows)

gene_dfs = host_gene_dfs + symb_gene_dfs
gene_df = pd.concat(gene_dfs, axis=0, ignore_index=True)
rows = gene_df.to_dict(orient="records")
for row in rows:
    row["experiment_id"] = args.experiment_id
    del row["transcriptome"]
    del row["species_prefix"]
SequenceCountsTable.write_tsv(args.output_dir+"/gene_counts.tsv", rows)

transcript_dfs = host_transcript_dfs + symb_transcript_dfs
transcript_df = pd.concat(transcript_dfs, axis=0, ignore_index=True)
rows = transcript_df.to_dict(orient="records")
for row in rows:
    row["experiment_id"] = args.experiment_id
    del row["transcriptome"]
    del row["species_prefix"]
SequenceCountsTable.write_tsv(args.output_dir+"/transcript_counts.tsv", rows)
