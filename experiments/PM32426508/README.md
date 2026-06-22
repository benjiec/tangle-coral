
## DESeq2

Gene level

```
coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession x_a_tenuis \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession x_c_goreaui \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`
```

Transcript level

```
coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession x_a_tenuis \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/transcript_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession x_c_goreaui \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/transcript_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`
```

Merge

```
coral-py coral/scripts/analysis/des2-merge.py EXP_PM32426508 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508` \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM32426508`/deseq2_*.tsv
```
