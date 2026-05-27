
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
coral-py coral/experiments/PM34593802/unique-acc.py \
  pseudodiploria_symb.fna.gz
```

Create a Pile directory for this project's transcriptomes

```
mkdir -p <pile_parent_dir>/PM34593802/transcriptomes
```

Place the 6 files into 6 separate transcriptome directories. Run clustering
using Pile. Then make sure each directory has a `transcripts.fna` file
(clustered) and `transcripts.unclustered.fna` file. They can be with `.gz`
suffixers. The directories are

```
breviolum_b5
breviolum_faviinorum
orbicella_faveolata
pseudodiploria_clivosa
siderastrea_radians
symbiodinium_a3
```

See Pile README. Run TransDecoder, then convert the protein GFF3 file into
Tangle detected format, as `detected.tsv` in each of the 6 transcriptome dirs. E.g.

```
PILE_WORKSPACE=PM34593802 pile-py pile/transdecoder_to_detected.py \
  orbicella_faveolata
```

Combine the 6 .tsv detected files into a single `transcript_proteins.tsv`.
Combine the 6 .faa output files into a single `proteins.faa.gz`, and all 6
.fna.gz output files into a single `transcripts.fna.gz`.

These three combined files should then be moved to `tangle-py
tangle/scripts/defaults.py -m area_experiment PM34593802`.

Generate transcript manifest, like this

```
coral-py coral/scripts/analysis/transcript-manifest.py EXP_PM34593802 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802` \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/transcripts.fna.gz
```

## Quantification

See Pile README on downloaded SRA assets. Then use Pile to quantify. See
`quant.sh` for exact commands for each transcriptome.

Process the quants using the following script, to generate a `gene_counts.tsv`
file and a `transcript_genes.tsv` file.

```
coral-py coral/experiments/PM34593802/process_salmon_quants.py EXP_PM34593802 \
  <pile workbase dir>/quants \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`
```


## Classification

Perform KO and Pfam detection on protein fasta sequence, using heap-py
commands, on GCloud. Don't forget to also filter KO classifications further to
those meeting KO HMM thresholds.


## DESeq2

Gene level

```
coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession x_breviolum_b5 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession x_breviolum_faviinorum \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession x_symbiodinium_a3 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession x_orbicella_faveolata \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession x_pseudodiploria_clivosa \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`

coral-py coral/scripts/analysis/des2-simple.py --timepoint 0 --min-count 50 \
  --genome-accession x_siderastrea_radians \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/gene_counts.tsv \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`
```

Merge

```
coral-py coral/scripts/analysis/des2-merge.py EXP_PM34593802 \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802` \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/deseq2_*.tsv
```
