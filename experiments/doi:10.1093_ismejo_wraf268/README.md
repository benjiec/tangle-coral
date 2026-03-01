
## Transcript file

ssa01.fna file is downloaded from NCBI: https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=GLGU01

SRA files were downloaded from NCBI. Salmon was used to map SRA files onto the
ssa01.fna.

## Protein / ORF detection

ssa01.fna is nucleic transcripts - we need proteins for our analysis. The
following was used to detect proteins in ssa01.fna.

```
pip3 install orfipy
orfipy ssa01.fna --pep ssa01.faa --min 200 --procs 4 --start ATG
```

200 min was used because the paper says they used a minimum 200 for transcript
length, during trinity assembly.

the longest protein for each transcript was then selected using the
select_orfs.py script, this resulted in the final proteins.faa file in the
data/exp_results directory. This .faa file was then classified by KO and Pfam
HMMs.

Comparing Pfam identification on longest ORFs vs all ORFs: 47130 Pfam matches
on longest ORFs, and 48702 on all ORFs or 3% higher. Making an assumption that
this is not a big deal and going with longest ORFs to simplify processing.

## DESeq2

```
python3 scripts/analysis/des2-simple.py --timepoint 192 --min-count 5 \
  data/exp_results/doi:10.1093_ismejo_wraf268/sequence_data.tsv data/exp_results/doi:10.1093_ismejo_wraf268

python3 scripts/analysis/des2-simple.py --cohort 34C --min-count 5 \
  --include-timepoint 0 \
  data/exp_results/doi:10.1093_ismejo_wraf268/sequence_data.tsv data/exp_results/doi:10.1093_ismejo_wraf268

python3 scripts/analysis/des2-specific.py --cohort1 34C --timepoint1 0 --cohort2 27C --timepoint2 192 --min-count 5 \
  data/exp_results/doi:10.1093_ismejo_wraf268/sequence_data.tsv data/exp_results/doi:10.1093_ismejo_wraf268

PYTHONPATH=. python3 scripts/analysis/des2-merge.py \
  data/exp_results/doi:10.1093_ismejo_wraf268 \
  data/exp_results/doi:10.1093_ismejo_wraf268/proteins.faa \
  data/exp_results/doi:10.1093_ismejo_wraf268/deseq2_*.tsv
```
