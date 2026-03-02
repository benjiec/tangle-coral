
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

Raw sequencing data were downloaded from SRA using sratools, e.g.

```
fasterq-dump SRR6255822 --progress
```

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

Cluster each of the 6 transcripts using MMSeqs. This step crushes/aggregates
genes on different alleles but provides more statistical power for differential
analysis through this aggregation.

```
scripts/cluster/mmseqs-cluster-trinity-transcripts \
  data/exp_results/doi:10.1038_s41467-021-25950-4/orbicella_host.fna.gz
scripts/cluster/mmseqs-cluster-trinity-transcripts \
  data/exp_results/doi:10.1038_s41467-021-25950-4/orbicella_symb.fna.gz

scripts/cluster/mmseqs-cluster-trinity-transcripts \
  data/exp_results/doi:10.1038_s41467-021-25950-4/pseudodiploria_host.fna.gz
scripts/cluster/mmseqs-cluster-trinity-transcripts \
  data/exp_results/doi:10.1038_s41467-021-25950-4/pseudodiploria_symb.fna.gz

scripts/cluster/mmseqs-cluster-trinity-transcripts \
  data/exp_results/doi:10.1038_s41467-021-25950-4/siderastrea_host.fna.gz
scripts/cluster/mmseqs-cluster-trinity-transcripts \
  data/exp_results/doi:10.1038_s41467-021-25950-4/siderastrea_symb.fna.gz
```


## Read mapping

Create Salmon index, then use Salmon to map reads to transcriptomes, e.g.

```
salmon index -t pseudodiploria_symb.fna.gz_rep_seq.fna -i pseudodiploria_symb.salmon_index
salmon quant -i pseudodiploria_symb.salmon_index \
  -l A --validateMappings -o symb_quants/SRR6255880 -p 2\
  -1 SRR6255880_1.fastq -2 SRR6255880_2.fastq
```

See `mk_map_cmds.py` for a list of Salmon commands used on the SRA files.

`process_salmon_quants.py` was used to summarize the quantifications into
`sequence_data.full.tsv`. This TSV is rather large, and slow to load. This
script not only rolls up the quantifications into a nice table, but also does
aggregation of both counts and TPM by "gene" in case some isoforms are still
present after mmseqs clustering.


## Mapping to KO and Pfam

Use `orfipy` to translate and compute ORFs for each of the 6 .fna files, e.g.

```
orfipy orbicella_host.fna.gz_rep_seq.fna.gz --pep orbicella_host.faa --min 150 --procs 4 --start ATG
```

Combine the 6 .faa output files into a single `proteins.faa`.

Perform KO and Pfam detection on protein fasta sequence

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
TODO deduplicate by gene, then remove _ORF so entries match sequence_list.tsv
TODO run scripts/analysis/assign-ko.py on KO TSV


## DESeq2

TODO quantify reads then run DESeq2

```
python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --genome-accession doi:10.1038_s41467-021-25950-4_breviolum_b5 \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --genome-accession doi:10.1038_s41467-021-25950-4_breviolum_faviinorum \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --genome-accession doi:10.1038_s41467-021-25950-4_symbiodinium_a3 \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --genome-accession doi:10.1038_s41467-021-25950-4_orbicella_faveolata \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --genome-accession doi:10.1038_s41467-021-25950-4_pseudodiploria_clivosa \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

python3 scripts/analysis/des2-simple.py --timepoint 0 --min-count 5 \
  --genome-accession doi:10.1038_s41467-021-25950-4_siderastrea_radians \
  data/exp_results/doi:10.1038_s41467-021-25950-4/sequence_data.tsv data/exp_results/doi:10.1038_s41467-021-25950-4

PYTHONPATH=. python3 scripts/analysis/des2-merge.py \
  data/exp_results/doi:10.1038_s41467-021-25950-4 \
  data/exp_results/doi:10.1038_s41467-021-25950-4/proteins.faa \
  data/exp_results/doi:10.1038_s41467-021-25950-4/deseq2_*.tsv
```
