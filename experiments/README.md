# Analysis Tasks

## Processing RNAseq Data, Quantification

See https://github.com/benjiec/tangle-pile

In experiment specific directories under the pile directory structure, generate
quantification files as well as proteome (i.e. .faa) files (possibly more than
one, for different species in the same experiment).

See PM34593802/README.md for an example workflow.


## Running DESeq2

DESeq2 analysis on RNAseq or proteomics results can be done using two scripts.
First, the `des2-simple.py` script does simple comparison of every two types of
cohorts at a given timepoint, or every two timepoints at the same cohort. For
example,

```
coral-py coral/scripts/analysis/des2-simple.py \
  --timepoint 1 --min-count 5 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`

coral-py coral/scripts/analysis/des2-simple.py \
  --cohort SS8 --min-count 5 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`
```

Use the `--include-timepoint 0` flag with the `--cohort` flag to limit to
specific timepoints, when comparing timepoint data within a cohort.

There is also a `des2-specific.py` script to compare any two cohort and/or
timepoints

```
coral-py coral/scripts/analysis/des2-specific.py --cohort1 34C --timepoint1 0 --cohort2 27C --timepoint2 192 --min-count 5 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`
```

Use the following to merge multiple DES2 TSV files into a single tall TSV file,
first argument is output directory.

```
coral-py coral/scripts/analysis/des2-merge.py \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508` \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/deseq2_*.tsv
```


## Select Transcripts for Further Processing

Use the following to get a list of top transcripts. The last line specifies a
file of past results so we don't re-classify sequence IDs already classified.

```
coral-py coral/scripts/analysis/top-sequences.py \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/des2_tall.tsv \
  --l2fc-threshold 1 --padj-threshold 0.05 \
  --results-fn `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_fs.tsv
```

And to do the above, then use results to filter a fasta, add

```
  | tangle-py tangle/scripts/fasta-emit.py \
      `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/transcripts.fna.gz \
      --prefix-with-underscore -
```
