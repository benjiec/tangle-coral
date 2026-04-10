# Analysis Tasks

## Processing RNAseq Data and Quantification

See https://github.com/benjiec/tangle-pile

In experiment specific directories under the pile directory structure, generate
quantification files as well as proteome (i.e. .faa) files (possibly more than
one, for different species in the same experiment).


## Process Quants

Each directory must have its own script to process quants into
`sequence_data.tsv`.

You can put multiple genomes in a single dataset file, i.e.
`sequence_data.tsv`. This is preferred if multiple species are present in a
sample, as DESeq2 can normalize across samples and more data helps with
normalization. However, make sure, in the mapping workflow, that sequence IDs
are unique between species.

```
coral-py coral/experiments/PM32426508/process_salmon_quants.py \
  experiments/PM32426508 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`
```

## Classify Transcripts

Use heap-py and tangle-py commands from ../README.md to classify transcripts to
KOs and Pfam domains. For KO, remember to run the assignment step afterward.

As inputs, use the amino acid sequences of ORFs predicted by TransDecoder -- or
pile proteome files.

This section should generate

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
coral-py coral/scripts/analysis/des2-simple.py \
  --timepoint 1 --min-count 5 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/sequence_data.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`

coral-py coral/scripts/analysis/des2-simple.py \
  --cohort SS8 --min-count 5 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/sequence_data.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`
```

Use the `--include-timepoint 0` flag with the `--cohort` flag to limit to
specific timepoints, when comparing timepoint data within a cohort.

There is also a `des2-specific.py` script to compare any two cohort and/or
timepoints

```
coral-py coral/scripts/analysis/des2-specific.py --cohort1 34C --timepoint1 0 --cohort2 27C --timepoint2 192 --min-count 5 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_data.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`
```

Use the following to merge multiple DES2 TSV files into a single tall TSV file,
first argument is output directory. IMPORTANT: for RNAseq data, make sure the
mapped file from Salmon, or other tool, has been converted to refer to protein
sequence IDs in the proteins.faa file.

```
coral-py coral/scripts/analysis/des2-merge.py \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508` \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/proteins.faa \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/deseq2_*.tsv
```
