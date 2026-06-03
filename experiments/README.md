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

Use the following to get a list of top transcripts. The `--results-fn` line
specifies a file of past results so we don't re-classify sequence IDs already
classified. The `--transcript-genes-fn` and `--transcript-proteins-fn` lines
convert gene IDs to transcript IDs and protein IDs respectively using the
lookup files; ignoring one or all of them will leave the output IDs at the last
available conversion or just gene ids.

```
coral-py coral/scripts/analysis/top-sequences.py \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/des2_tall.tsv \
  --l2fc-threshold 0.1 --padj-threshold 0.05 --mean-threshold 200 \
```

And you can pipe the results from above to filter a protein fasta, for further
classification

```
  | tangle-py tangle/scripts/fasta-emit.py --prefix-with-underscore \
      `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/proteins.faa.gz - > top.faa
```


## Clustering of top transcripts

Add transcripts to a pool of sequences that have database name in the accessions

```
tangle-py tangle/scripts/fasta-pool.py \
  --database EXP_PM34593802 \
  top.faa \
  --append-to all_interesting_sequences.faa
```

Cluster using MMSeqs

```
docker run --rm \
  -v .:/work \
  ghcr.io/soedinglab/mmseqs2 \
  easy-cluster /work/all_interesting_sequences.faa /work/cluster /tmp \
  --cov-mode 0 -c 0.8 --min-seq-id 0
```

Prepare a cluster TSV that can be loaded into BigQuery

```
tangle-py tangle/scripts/demux-mmseq-clusters.py \
  --clustering-description glbtx \
  --parameters "c0.8 sid0" \
  --cluster-type aa \
  cluster_cluster.tsv clusters.tsv
```

Look at an alignment by cluster name

```
tangle-py tangle/scripts/cluster-align.py \
  cluster_all_seqs.fasta out.faa \
  ea23f31d37
```

The resulting alignment file `out.faa` can be visaulized at https://alignmentviewer.org/
