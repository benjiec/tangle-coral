Raw seq reads were downloaded from NCBI using BioProject code referenced in
paper. 

Decoy and transcript fasta files for C. goreaui and A. tenius were constructed.
Salmon was used to map reads onto these transcripts.

DESeq2 performs count normalization by default.

```
python3 scripts/analysis/des2-simple.py --timepoint 1 --min-count 5 \
  --genome-accession GCA_014633955.1 \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv data/exp_results/doi:10.1126_sciadv.aba2498

python3 scripts/analysis/des2-simple.py --timepoint 1 --min-count 5 \
  --genome-accession GCA_947184155.2 \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv data/exp_results/doi:10.1126_sciadv.aba2498

PYTHONPATH=. python3 scripts/analysis/des2-merge.py \
  data/exp_results/doi:10.1126_sciadv.aba2498 \
  data/exp_results/doi:10.1126_sciadv.aba2498/proteins.faa \
  data/exp_results/doi:10.1126_sciadv.aba2498/deseq2_*.tsv
```
