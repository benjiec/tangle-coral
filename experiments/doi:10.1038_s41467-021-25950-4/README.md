
Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC8484447/

## Data

Zip file of data was downloaded from Dryad repository using accession
10.5061/dryad.k3j9kd57b. The zip file contains Trinity assembled .fna files
used for mapping.

Paper says

"In order to generate a de novo reference metatranscriptome for each holobiont
species, the cleaned and filtered reads from the replicates from control and
treatment fragments per each species were pooled and assembled using Trinity
(v2.4.0)72, generating three metatranscriptomes (O. faveolata, S. radians, and
P. clivosa)."

Raw sequencing data were downloaded from SRA using sratools, e.g.

```
fasterq-dump SRR6255822 --progress
```

The following SRA entries were used:

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

DESeq2 results were downloaded from supplemental data 1 from the paper
[https://pmc.ncbi.nlm.nih.gov/articles/PMC8484447/]. We re-mapped the reads and
re-generated the DESeq2 results, but are using the author's DESeq2 results for
reference and sanity check.


## Read mapping

For now just using host and algae transcriptomes.

First rewrite the .fna.gz files to have a prefix for each entry, so we ensure
all entry names are unique across the three holobionts. E.g. do the following
for the 6 files (3 hosts, 3 symbionts)

```
PYTHONPATH=. python3 experiments/doi:10.1038_s41467-021-25950-4/unique-acc.py \
  data/exp_results/doi:10.1038_s41467-021-25950-4/pseudodiploria_symb.fna.gz 
```

Create Salmon index, then use Salmon to map reads to transcriptomes, e.g.

```
salmon index -t pseudodiploria_symb.fna -i pseudodiploria_symb.salmon_index
salmon quant -i pseudodiploria_symb.salmon_index \
  -l A --validateMappings -o symb_quants/SRR6255880 -p 2\
  -1 SRR6255880_1.fastq -2 SRR6255880_2.fastq
```

See `mk_map_cmds.py` for a list of Salmon commands used on the SRA files.

`process_salmon_quants.py` was used to summarize the quantifications into
`sequence_data.full.tsv`. This TSV is rather large, and slow to load.


## Mapping to KO and Pfam

Use `orfipy` to translate and compute ORFs for each of the 6 .fna files, e.g.

```
orfipy orbicella_host.fna.gz --pep orbicella_host.faa --min 150 --procs 4 --start ATG
```

Combine the 6 .faa output files into a single `proteins.faa`, then perform KO
and Pfam detection on this file.

```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --cpu 2 --disable-cutoff-ga --genome-accession _ \
  --fasta-file experiments/doi:10.1038_s41467-021-25950-4/proteins.faa \
  kegg-downloads/ko.hmm _ data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_ko.orf.tsv

PYTHONPATH=. python3 scripts/classify/classify.py \
  --cpu 2 --genome-accession _ \
  --fasta-file experiments/doi:10.1038_s41467-021-25950-4/proteins.faa \
  pfam-downloads/Pfam-A.hmm _ data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_pfam.orf.tsv
```

TODO map to KO and Pfam
TODO remove _ORF suffix from classification .tsv files

## DESeq2

TODO find internal control genes
TODO run DESeq2 by each genome, total 6

```
python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --control-sequence SYMB_CONTROL \
  --genome-accession doi:10.1038_s41467-021-25950-4_breviolum_b5 \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --control-sequence SYMB_CONTROL \
  --genome-accession doi:10.1038_s41467-021-25950-4_breviolum_faviinorum \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --control-sequence SYMB_CONTROL \
  --genome-accession doi:10.1038_s41467-021-25950-4_symbiodinium_a3 \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --control-sequence HOST_CONTROL \
  --genome-accession doi:10.1038_s41467-021-25950-4_orbicella_faveolata \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --control-sequence HOST_CONTROL \
  --genome-accession doi:10.1038_s41467-021-25950-4_pseudodiploria_clivosa \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --control-sequence HOST_CONTROL \
  --genome-accession doi:10.1038_s41467-021-25950-4_siderastrea_radians \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4
```

