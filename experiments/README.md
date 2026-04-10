# Analysis Tasks

## Processing RNAseq Data and Quantification

See https://github.com/benjiec/pile


## Process Quants

Each directory must have its own script to process quants into
`sequence_data.tsv`.

You can put multiple genomes in a single dataset file, i.e.
`sequence_data.tsv`. This is preferred if multiple species are present in a
sample, as DESeq2 can normalize across samples and more data helps with
normalization. However, make sure, in the mapping workflow, that sequence IDs
are unique between species.

```
python3 experiments/doi:10.1126_sciadv.aba2498/process_salmon_quants.py \
  experiments/doi:10.1126_sciadv.aba2498 data/exp_results/doi:10.1126_sciadv.aba2498
```

## Classify Transcripts

Use heap-py and tangle-py commands from ../README.md to classify transcripts to
KOs and Pfam domains. For KO, remember to run the assignment step afterward.

These commands should generate

  * `sequence_ko.tsv`: sequence to KO mapping
  * `sequence_pfam.tsv`: sequence to Pfam mapping


## Run DESeq2

This section generates

  * `sequence_list.tsv`: from des2-merge script below
  * `des2_tall.tsv`: des2 data, from des2-merge script below


DESeq2 analysis on RNAseq or proteomics results can be done using two scripts.
First, the `des2-simple.py` script does simple comparison of every two types of
cohorts at a given timepoint, or every two timepoints at the same cohort. For
example,

```
python3 scripts/analysis/des2-simple.py \
  --timepoint 1 --min-count 5 \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv data/exp_results/doi:10.1126_sciadv.aba2498

python3 scripts/analysis/des2-simple.py \
  --cohort SS8 --min-count 5 \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv data/exp_results/doi:10.1126_sciadv.aba2498
```

Use the `--include-timepoint 0` flag with the `--cohort` flag to limit to
specific timepoints, when comparing timepoint data within a cohort.

There is also a `des2-specific.py` script to compare any two cohort and/or
timepoints

```
python3 scripts/analysis/des2-specific.py --cohort1 34C --timepoint1 0 --cohort2 27C --timepoint2 192 --min-count 5 \
  data/exp_results/doi:10.1093_ismejo_wraf268/sequence_data.tsv data/exp_results/doi:10.1093_ismejo_wraf268
```

Use the following to merge multiple DES2 TSV files into a single tall TSV file,
first argument is output directory. IMPORTANT: for RNAseq data, make sure the
mapped file from Salmon, or other tool, has been converted to refer to protein
sequence IDs in the proteins.faa file.

```
PYTHONPATH=. python3 scripts/analysis/des2-merge.py \
  data/exp_results/doi:10.1126_sciadv.aba2498 \
  data/exp_results/doi:10.1126_sciadv.aba2498/proteins.faa \
  data/exp_results/doi:10.1126_sciadv.aba2498/deseq2_*.tsv
```
