import re
import pandas as pd
import argparse
from itertools import combinations
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

ap = argparse.ArgumentParser()
ap.add_argument("data_tsv")
ap.add_argument("output_dir")
ap.add_argument("--control-sequence", type=str, default=None, action="append")
ap.add_argument("--genome-accession", type=str, default=None)
ap.add_argument("--min-count", type=int, default=10)
ap.add_argument("--cohort", type=str, default=None)
ap.add_argument("--timepoint", type=str, default=None)
args = ap.parse_args()

cohort = args.cohort
timepoint = args.timepoint
if cohort and timepoint:
    raise Exception("--cohort and --timepoint: please specify only one")
if not cohort and not timepoint:
    raise Exception("--cohort and --timepoint: please specify exactly one")

tall_df = pd.read_csv(args.data_tsv, delimiter='\t')

if args.genome_accession:
    print(f"filtering by {args.genome_accession}")
    tall_df = tall_df[tall_df['genome_accession'] == args.genome_accession]

cond_name = None

if cohort is not None:
    tall_df = tall_df[tall_df['cohort'] == cohort].copy()
    tall_df['condition'] = tall_df['timepoint']
    cond_name = f"cohort_{cohort}_timepoint"
else:
    tall_df = tall_df[tall_df['timepoint'].astype(str) == timepoint].copy()
    tall_df['condition'] = tall_df['cohort']
    cond_name = f"timepoint_{timepoint}_cohort"

if args.control_sequence:
    for ctrl in args.control_sequence:
        ctrl = re.sub(r'\W', '_', ctrl)
        cond_name = f"ctrl_{ctrl}_{cond_name}"

if args.genome_accession:
    cond_name = f"{args.genome_accession}_{cond_name}"

# generate unique sample-timepoint column
tall_df['cohort_timepoint_sample'] = tall_df['sample'].astype(str) + '/' + tall_df['timepoint'].astype(str) + '/' + tall_df['cohort'].astype(str)

# split into counts and metadata dataframes
metadata_df = tall_df[['cohort_timepoint_sample', 'condition']].drop_duplicates()
counts_df = tall_df[['cohort_timepoint_sample', 'sequence_id', 'count']]

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

stable_controls = None
if args.control_sequence:
    assert type(args.control_sequence) in (list, tuple)
    stable_controls = args.control_sequence
    print("using controls", stable_controls)

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata_df,
    design="~condition",
    refit_cooks=True,
    inference=inference,
    control_genes=stable_controls
)

dds.deseq2()

# print(dds)
# print(dds.var["dispersions"])
# print(dds.varm["LFC"])

unique_values = metadata_df['condition'].unique()
pairs = list(combinations(unique_values, 2))

for pair in pairs:
    fn = "_".join([str(x) for x in pair])
    fn = f"{args.output_dir}/deseq2_{cond_name}_{fn}.tsv"
    print("generating", fn)

    contrast = ["condition"]+list(pair)
    ds = DeseqStats(dds, contrast=contrast, inference=inference, quiet=True)
    ds.summary()
    ds.results_df.to_csv(fn, sep='\t')
