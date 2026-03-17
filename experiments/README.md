# Analysis Tasks

## Processing RNAseq Data and Quantification

See https://github.com/benjiec/pile


## Process Quants

Each directory must have its own script to process quants into `sequence_data.tsv`.

You can put multiple genomes in a single dataset file, i.e.
`sequence_data.tsv`. This is preferred if multiple species are present in a
sample, as DESeq2 can normalize across samples, based on more data per sample.
However, make sure, in the mapping workflow, that sequence IDs are unique
between species.

```
python3 experiments/doi:10.1126_sciadv.aba2498/process_salmon_quants.py \
  experiments/doi:10.1126_sciadv.aba2498 data/exp_results/doi:10.1126_sciadv.aba2498
```

Helpers to update the genome accession to something more formal

```
python3 scripts/data/update-genome-accession.py \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv \
  GCA_014633955.1 --when doi:10.1126/sciadv.aba2498-a_tenuis
python3 scripts/data/update-genome-accession.py \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv \
  GCA_947184155.2 --when doi:10.1126/sciadv.aba2498-c_goreaui
```


## Classify Transcripts

```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --disable-cutoff --genome-accession _ \
  --fasta-file experiments/10.1126_sciadv.aba2498/c_goreaui.faa \
  kegg-downloads/ko.hmm _ experiments/10.1126_sciadv.aba2498/c_goreaui_ko.tsv

PYTHONPATH=. python3 scripts/classify/classify.py \
  --disable-cutoff --genome-accession _ \
  --fasta-file experiments/10.1126_sciadv.aba2498/aten.faa \
  kegg-downloads/ko.hmm _ experiments/10.1126_sciadv.aba2498/aten_ko.tsv

cat experiments/10.1126_sciadv.aba2498/c_goreaui_ko.tsv experiments/10.1126_sciadv.aba2498/aten_ko.tsv > experiments/10.1126_sciadv.aba2498/sequence_ko.tsv
```

Then filter the classified file to only keep classifications above the KO
threshold. For example,

```
python3 scripts/analysis/assign-ko.py experiments/10.1126_sciadv.aba2498/sequence_ko.tsv
mv experiments/10.1126_sciadv.aba2498/sequence_ko.tsv_filtered experiments/10.1126_sciadv.aba2498/sequence_ko.tsv
```

You can also directly classify/identify Pfam domains from the proteome

```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --genome-accession _ \
  --fasta-file data/exp_results/doi:10.1126_sciadv.aba2498/proteins.faa \
  pfam-downloads/Pfam-A.hmm _ data/exp_results/doi:10.1126_sciadv.aba2498/sequence_pfam.tsv \
   --cpu 4
```

Then, copy the .ko.tsv and .faa files to
data/exp_results/10.1126_sciadv.aba2498 directory. Usually, following these
conventions

  * `sequence_data.tsv`: data file
  * `sequence_ko.tsv`: sequence to KO mapping
  * `sequence_pfam.tsv`: sequence to Pfam mapping
  * `sequence_list.tsv`: from des2-merge script below
  * `des2_tall.tsv`: des2 data, from des2-merge script below


## Run DESeq2

First, identify an internal control gene, EIF5B has been shown to work well for
algae and likely coral as well.

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
