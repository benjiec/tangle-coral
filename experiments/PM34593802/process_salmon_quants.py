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
    transcriptome=desc.split(", ")[1].lower().replace(" ", "_"),
    genome_accession="x_"+desc.split(", ")[1].lower().replace(" ", "_")
  )
  for fn,desc in md.items()
}

host_md = {
  fn: dict(
    timepoint=0,
    cohort=desc.split(", ")[2],
    sample=desc.split(", ")[3],
    transcriptome=desc.split(", ")[0].lower().replace(" ", "_"),
    genome_accession="x_"+desc.split(", ")[0].lower().replace(" ", "_")
  )
  for fn,desc in md.items()
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
SequenceCountsTable.write_tsv(args.output_dir+"/gene_counts.tsv", rows)

transcript_dfs = host_transcript_dfs + symb_transcript_dfs
transcript_df = pd.concat(transcript_dfs, axis=0, ignore_index=True)
rows = transcript_df.to_dict(orient="records")
for row in rows:
    row["experiment_id"] = args.experiment_id
    del row["transcriptome"]
SequenceCountsTable.write_tsv(args.output_dir+"/transcript_counts.tsv", rows)
