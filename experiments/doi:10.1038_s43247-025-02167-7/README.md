
```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --disable-cutoff-ga --max-rank 10 \
  --genome-accession x2 \
  --fasta-file peptides.faa kegg-downloads/ko.hmm x2 peptide_ko.tsv
```

How to specify genome accession? phylogeny, etc?
