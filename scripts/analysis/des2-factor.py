import pandas as pd
import argparse
from itertools import combinations
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

ap = argparse.ArgumentParser()
ap.add_argument("data_tsv")
ap.add_argument("cohort_metadata_tsv")
ap.add_argument("output_dir")
ap.add_argument("--min-count", type=int, default=10)
ap.add_argument("--factor", type=str, required=True)
ap.add_argument("--category", type=str, required=True)
args = ap.parse_args()

def compute_by_category(cond_name, full_df, condition_df):

    data_df = pd.merge(full_df, condition_df, how="inner", on="cohort")
    print(f"DES2 analysis on {cond_name}: {full_df.shape}->{data_df.shape}")
    data_df['cohort_timepoint_sample'] = data_df['sample'].astype(str) + '/' + data_df['timepoint'].astype(str) + '/' + data_df['cohort'].astype(str)

    # split into counts and metadata dataframes
    metadata_df = data_df[['cohort_timepoint_sample', 'condition']].drop_duplicates()
    counts_df = data_df[['cohort_timepoint_sample', 'sequence_id', 'count']]

    # set both to be indexed by sample
    metadata_df = metadata_df.set_index('cohort_timepoint_sample')
    counts_df = counts_df.set_index('cohort_timepoint_sample')

    # convert to wide with sequence_id as column, sample as row
    counts_df = counts_df.pivot(columns='sequence_id', values='count')

    # filter out genes that have less than min_count read counts in total
    sequence_to_keep = counts_df.columns[counts_df.sum(axis=0) >= args.min_count]
    counts_df = counts_df[sequence_to_keep]
    counts_df = counts_df.fillna(0.0).round().astype(int)

    # get them in the same index order
    metadata_df = metadata_df.loc[counts_df.index]

    print("metadata", metadata_df.shape)
    print("counts", counts_df.shape)

    inference = DefaultInference(n_cpus=2)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design="~condition",
        refit_cooks=True,
        inference=inference,
    )

    dds.deseq2()

    # print(dds)
    # print(dds.var["dispersions"])
    # print(dds.varm["LFC"])

    unique_values = metadata_df['condition'].unique()
    pairs = list(combinations(unique_values, 2))

    for pair in pairs:
        fn = "_".join(pair)
        fn = f"{args.output_dir}/deseq2_{cond_name}_{fn}.tsv"
        print("generating", fn)

        contrast = ["condition"]+list(pair)
        ds = DeseqStats(dds, contrast=contrast, inference=inference, quiet=True)
        ds.summary()
        ds.results_df.to_csv(fn, sep='\t')


data_df = pd.read_csv(args.data_tsv, delimiter='\t')
metadata_df = pd.read_csv(args.cohort_metadata_tsv, delimiter='\t')

# within each category, we are comparing every possible combination of factors

category_values = list(metadata_df[metadata_df['term'] == args.category]['value'].drop_duplicates())

for category_value in category_values:
    print(f"{args.category}={category_value}")
    cohorts = metadata_df[(metadata_df['term'] == args.category) & (metadata_df['value'] == category_value)]['cohort'].unique()
    condition_df = metadata_df[metadata_df['cohort'].isin(cohorts)]
    condition_df = condition_df[condition_df['term'] == args.factor].copy()
    condition_df['condition'] = condition_df['value']
    condition_df = condition_df[['cohort', 'condition']].drop_duplicates()
    cond_name = f"{args.category}_{category_value}_{args.factor}"
    compute_by_category(cond_name, data_df, condition_df)
