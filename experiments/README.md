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

Cluster using MMSeqs, for AA sequences

```
docker run --rm \
  -v .:/work \
  ghcr.io/soedinglab/mmseqs2 \
  easy-cluster /work/top_pooled.faa /work/cluster /tmp \
  --cov-mode 0 -c 0.8 --min-seq-id 0 -s 4
```

Prepare a cluster TSV that can be loaded into BigQuery

```
tangle-py tangle/scripts/demux-mmseq-clusters.py \
  --clustering-description glbtx \
  --parameters "c_m=0,c=0.8,s_id=0,s=4" \
  --cluster-type aa \
  cluster_cluster.tsv clusters.tsv
```

Look at an alignment by cluster name

```
tangle-py tangle/scripts/cluster-align.py \
  --clustering-description glbtx \
  --parameters "c_m=0,c=0.8,s_id=0,s=4" \
  --cluster-type aa \
  cluster_all_seqs.fasta out.faa \
  ea23f31d37
```

The resulting alignment file `out.faa` can be visaulized at https://alignmentviewer.org/

3Di search for KOs

```
docker run --platform linux/amd64 --rm \
  -v ./query_db_top_tx:/query \
  -v ./target_db_ko:/target \
  -v ./res:/res ghcr.io/steineggerlab/foldseek \
  search \
  /query/final_db \
  /target/final_db \
  /res/result_db \
  /tmp \
  --target-search-mode 1 \
  --alignment-type 2 \
  -e 0.001 \
  -s 9.5 \
  --cov-mode 0 -c 0.0

docker run --platform linux/amd64 --rm \
  -v ./query_db_top_tx:/query \
  -v ./target_db_ko:/target \
  -v ./res:/res ghcr.io/steineggerlab/foldseek \
  convertalis \
  /query/final_db \
  /target/final_db \
  /res/result_db \
  /res/results.tsv \
  --format-output "query,target,evalue,bits,qstart,qend,tstart,tend"

# Generates file that can be uploaded
heap-py heap/scripts/foldseek.py \
  --input-from-foldseek 3di/res-cm0-c0/results.tsv \
  --keep-all-results \
  --query-database-name EXP_PM34593802 \
  --query-type protein \
  --target-database-name 3di-ko \
  --target-type protein \
  _ _ _ sequence_3di_ko.tsv

# Load to BQ

tangle-py tangle/scripts/bq-schema.py tangle.detected > tangle_detected.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_detected \
  sequence_3di_ko.tsv \
  ./tangle_detected.schema.json 
```

The following is not very useful, so it's here as notes. Don't waste time on
this at the moment.

3Di clustering using FoldSeek. First cd into the directory with the 3Di
database, then (note the parameter differences)

```
docker run --platform linux/amd64 --rm \
  -v .:/work \
  ghcr.io/steineggerlab/foldseek cluster \
  /work/final_db /work/cluster_db /tmp \
  --target-search-mode 1 \
  --cov-mode 0 -c 0.8 -s 5 --min-seq-id 0

docker run --platform linux/amd64 --rm -v .:/work ghcr.io/steineggerlab/foldseek createtsv \
    /work/final_db \
    /work/final_db \
    /work/cluster_db \
    /work/top_cluster_3di.txt
```

Prepare a cluster TSV that can be loaded into BigQuery

```
tangle-py tangle/scripts/demux-mmseq-clusters.py \
  --member-database EXP_PM34593802 \
  --clustering-description glbtx \
  --parameters "c_m=0,c=0.8,s_id=0,s=5" \
  --cluster-type 3di \
  top_cluster_3di.txt clusters_3di.tsv
```
