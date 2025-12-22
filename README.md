# Needle

Needle is a suite of tools to curate and compare sequences of proteins by
pathway, for understudied organisms. Cryptic gene structures, unconventional
splicing models, unfinished genomes, understudied proteomes, are some of the
challenges that make applying traditional workflows (e.g. gene prediction and
classification) difficult and/or unreliable. Needle focuses on working around
these challenges to detect presence of proteins for key pathways, and curates
information in a phylogentic aware manner to enable comparative analysis within
and across species.


## Setup

Install NCBI Docker image

```
docker pull ncbi/blast
```

Install MMSeqs2 Docker image

```
docker pull ghcr.io/soedinglab/mmseqs2
```

Setup SwissProt DB for MMSeqs2

```
scripts/mmseqs-swissprot-setup
```

Install Muscle Docker image

```
docker pull pegi3s/muscle
```

Install `HMMer` package. E.g. on MacOS run `brew install hmmer`.

Create Python virtualenv

```
python3 -m venv .venv
source .venv/bin/activate
pip3 install -r requirements.txt
```

## Initial Data

There are some initial data already in `data` directory. Below are instructions
to re-create them.

### Download a list of KO (KEGG Ortholog) numbers and names

```
curl https://rest.kegg.jp/list/ko -o data/ko.txt
echo "Ortholog ID\tOrtholog Name" | cat - data/ko.txt > data/ko.tsv
rm data/ko.txt
```

### Download a list of KEGG modules

```
curl https://rest.kegg.jp/list/module -o data/modules.txt
echo "Module ID\tModule Name" | cat - data/modules.txt > data/modules.tsv
rm data/modules.txt
```

### Fetch KO numbers for all the modules

The following creates `data/module_ko.tsv`

```
python3 scripts/fetch-kegg-module-ko.py
```

### Download KEGG KO profile HMMs

Download the HMM profiles from `https://www.genome.jp/ftp/db/kofam/`. The
`profiles.tar.gz` file is large, so this may take awhile.

Concatenate all the .hmm files together, e.g.

```
cat profiles/*.hmm > kegg_downloads/ko.hmm
```

Also, download the `ko_list.gz` file from the above location into
`data/ko_thresholds.gz`. This file contains scoring criteria for using the
HMMs.

### Download HMM profiles from Pfam

Download
`https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`. After
uncompress into `pfam-downloads` directory, run the following to create a
searchable HMM database

```
hmmpress pfam-downloads/Pfam-A.hmm
```

### Prepare List of Genomes

There is a list of Coral genomes in `data/genomes_coral.txt`, and a list of
Symbiodinium (algae) genomes in `data/genomes_algae.txt`. Run the following
command to generate `data/genomes.tsv`, which includes genome name and taxonomy
information.

```
PYTHONPATH=. python3 scripts/fetch-genomes.py data/genomes_coral.txt
```


## Workflow and Scripts

The general workflow looks like the following

  * Select a module to use; e.g. a module can be a pathway
  * Detect Proteins using ortholog profile HMMs: tblastn+HMM based refinement
  * Classify Proteins: hmmscan
  * Cluster Proteins: MMSeqs2
  * Generate MSAs: Muscle and hmmalign
  * Finish Protein Sequences
  * Improve Phylogene-aware family profiles

The following instructions use M00009 KEGG module, describing the TCA cycle, as
an example. Replace this number with other module IDs as appropriate.


### Generating Query .hmm for a KEGG Module

Use `hmmfetch` to create a smaller HMM database for the KOs of a module. Use
the following naming convention but change the module ID: `data/m00009_ko.hmm`.


### Detect Orthologs from Genomes

The following script puts outputs in `data/m00009_results` directory

```
./scripts/search-genome m00009 GCF_002042975.1
```


Or if you have a list of genome accessions in a file, e.g. `genomes.txt`, then do

```
./scripts/search-genomes m00009 genomes.txt
```

### Compare / Sanity check NCBI proteins against Needle detected proteins

Use the following script to compare, for a given HMM model, how NCBI annotated
proteins (i.e. in `protein.faa` and `genomic.gff`) compare against protein
found by Needle.

```
PYTHONPATH=. python3 scripts/compare-gff-with-match.py --best-hmm data/m00009_ko.hmm GCF_002042975.1 data/m00009_results/proteins.tsv \
  --output-file <filename>
```

### Classify Proteins

The following two commands will classify detected proteins first by KEGG ortholog, then Pfam domains.

```
PYTHONPATH=. python3 scripts/classify.py --disable-cutoff-ga data/m00009_ko.hmm m00009 data/m00009_results/classify.tsv
PYTHONPATH=. python3 scripts/classify.py pfam-downloads/Pfam-A.hmm m00009 data/m00009_results/classify.tsv
```

Classifications of annotated proteins submitted to NBCI can be classified in
the same way, and added to the same output TSV, using the following two
commands.

```
PYTHONPATH=. python3 scripts/classify-ncbi.py --disable-cutoff-ga data/m00009_ko.hmm GCF_932526225.1 data/m00009_results/classify.tsv
PYTHONPATH=. python3 scripts/classify-ncbi.py pfam-downloads/Pfam-A.hmm GCF_932526225.1 data/m00009_results/classify.tsv
```


### Cluster Proteins

Run the following script to cluster all the results from `search-genomes`

```
./scripts/cluster-dir m00009
```

### Generating MSAs

To generate MSAs and PNGs that visualize the MSAs, run the following script
for each module. This script runs the four sub-scripts below this, for each KO
number. If `cluster-ko` already ran, then the following script will also run
the two Muscle scripts for each of the clusters.

```
./scripts/generate-msas m00009
```

To generate MSAs, for a module and a KO, use the following scripts. Each
script puts a MSA in FASTA format in `data/m00009_results/m00009-alignments`

```
./scripts/muscle-ko m00009 <faa file prefix>
./scripts/hmmalign-ko m00009 <faa file prefix>
```

Note: *s (STOP codons) from BLAST search are preserved through MUSCLE by first
replacing them with Zs, then run MUSCLE, then replace Zs from the output of
MUSCLE back to *s.

SVG files (which can be opened via Chrome and other browsers) visualizing the
MSAs can generated with the following scripts.

```
./scripts/mk-msa-vis m00009 <faa file prefix> muscle
./scripts/mk-msa-vis m00009 <faa file prefix> hmmalign
```


### Other Scripts

Search in SwissProt for related proteins

```
scripts/mmseqs-swissprot-search proteins.faa results.tsv
```

Download files from NCBI

```
PYTHONPATH=. python3 scripts/ncbi-download.py GCF_932526225.1
```

