
Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC8484447/

## Trinity assembled transcripts and NGS reads

Zip file of data was downloaded from Dryad repository using accession
10.5061/dryad.k3j9kd57b. The zip file contains Trinity assembled .fna files by
species. The assemblies require some processing, see below.

Paper says

"In order to generate a de novo reference metatranscriptome for each holobiont
species, the cleaned and filtered reads from the replicates from control and
treatment fragments per each species were pooled and assembled using Trinity
(v2.4.0)72, generating three metatranscriptomes (O. faveolata, S. radians, and
P. clivosa)."

The SRA entries are

  * SRR6255820: Orbicella, control, A
  * SRR6255822: Orbicella, control, B
  * SRR6255825: Orbicella, control, C
  * SRR6255844: Orbicella, treatment, A
  * SRR6255862: Orbicella, treatment, B
  * SRR6255864: Orbicella, treatment, C
  * SRR6255880: Pseudodiploria, control, A
  * SRR6255882: Pseudodiploria, control, B
  * SRR6255888: Pseudodiploria, control, C
  * SRR6255887: Pseudodiploria, treatment, A
  * SRR6255886: Pseudodiploria, treatment, B
  * SRR6255890: Pseudodiploria, treatment, C
  * SRR6255889: Siderastrea, control, A
  * SRR6256312: Siderastrea, control, B
  * SRR6256323: Siderastrea, control, C
  * SRR6256322: Siderastrea, treatment, A
  * SRR6256324: Siderastrea, treatment, B
  * SRR6256325: Siderastrea, treatment, C	

DESeq2 results from the paper can be downloaded from supplemental data 1 from
the paper [https://pmc.ncbi.nlm.nih.gov/articles/PMC8484447/]. They can be used
as a reference.


## Transcript clustering

First rewrite the .fna.gz files to have a prefix for each entry, so we ensure
all entry names are unique across the three holobionts. E.g. do the following
for the 6 files (3 hosts, 3 symbionts)

```
PYTHONPATH=. python3 experiments/doi:10.1038_s41467-021-25950-4/unique-acc.py \
  data/exp_results/doi:10.1038_s41467-021-25950-4/pseudodiploria_symb.fna.gz
```

Use Pile to quantify, starting with the unique-fied transcripts.


## Annotation of KO and Pfam

Use `orfipy` to translate and compute ORFs for each of the 6 .fna clustered
sequences, e.g.

```
orfipy orbicella_host.fna.gz_rep_seq.fna.gz --pep orbicella_host.faa --min 150 --procs 4 --start ATG
```

Combine the 6 .faa output files into a single `proteins.faa`.

Perform KO and Pfam detection on protein fasta sequence, using tangle-py
commands. Don't forget to also filter KO classifications further to those
meeting KO HMM thresholds.

Then, aggregate hits on ORFs of isoforms of genes by gene. Basically, we want
classification outputs to use the same accessions as those used in Salmon and
DESeq2.

```
coral-py coral/scripts/analysis/classify-by-transcript.py sequence_ko_full.tsv
mv sequence_ko_full.tsv_aggregated sequence_ko.tsv
rm sequence_ko_full.tsv_aggregated

coral-py coral/scripts/analysis/classify-by-transcript.py sequence_pfam_full.tsv
mv sequence_pfam_full.tsv_aggregated sequence_pfam.tsv
```


## DESeq2

```
coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession doi:10.1038_s41467-021-25950-4_breviolum_b5 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_data.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession doi:10.1038_s41467-021-25950-4_breviolum_faviinorum \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_data.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession doi:10.1038_s41467-021-25950-4_symbiodinium_a3 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_data.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession doi:10.1038_s41467-021-25950-4_orbicella_faveolata \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_data.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession doi:10.1038_s41467-021-25950-4_pseudodiploria_clivosa \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_data.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession doi:10.1038_s41467-021-25950-4_siderastrea_radians \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_data.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-merge.py \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802` \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/proteins.faa \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/deseq2_*.tsv
```
