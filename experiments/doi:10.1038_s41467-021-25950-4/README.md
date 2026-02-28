
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

Salmon map reads to transcriptomes, e.g.

```
salmon index -t pseudodiploria_symb.fna -i pseudodiploria_symb.salmon_index
salmon quant -i pseudodiploria_symb.salmon_index \
  -l A --validateMappings -o quants/SRR6255880 -p 2\
  -1 SRR6255880_R1.fastq -2 SRR6255880_R2.fastq
```

TODO map reads onto .fna.gz files
TODO summarize quantifications into sequence_data.tsv file with genome accessions


## Mapping to KO and Pfam

TODO translate and select ORFs, create proteins.faa
TODO map to KO and Pfam


## DESeq2

TODO find internal control genes
TODO run DESeq2 by each genome, total 6
