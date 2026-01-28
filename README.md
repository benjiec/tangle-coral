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
scripts/data/mmseqs-swissprot-setup
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

To load this TSV file into Tableau, remove `"` and replace them with `''`.


### Download a list of KEGG modules

```
curl https://rest.kegg.jp/list/module -o data/modules.txt
echo "Module ID\tModule Name" | cat - data/modules.txt > data/modules.tsv
rm data/modules.txt
```

The following two scripts downloads KEGG module definitions and store them as a
list of KO numbers, and, in the second script, as steps and components.

```
python3 scripts/data/fetch-kegg-module-ko.py
PYTHONPATH=. python3 scripts/data/fetch-kegg-module-def.py
```

### Download KEGG KO profile HMMs

Download the HMM profiles from `https://www.genome.jp/ftp/db/kofam/`. The
`profiles.tar.gz` file is large, so this may take awhile.

Concatenate all the .hmm files together, e.g.

```
cat profiles/*.hmm > kegg-downloads/ko.hmm
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

Also, download the `Pfam-A.clans.tsv` file into data directory.

### Prepare List of Genomes

Run the following commands to generate `data/genomes.tsv`, which includes
genome name and taxonomy information for three sets of genomes.

```
PYTHONPATH=. python3 scripts/data/fetch-genomes.py data/genomes_detect.txt
PYTHONPATH=. python3 scripts/data/fetch-genomes.py data/genomes_ref.txt
```

Also download genomic sequences, using

```
PYTHONPATH=. python3 scripts/data/ncbi-download.py data/genomes_detect.txt
PYTHONPATH=. python3 scripts/data/ncbi-download.py data/genomes_ref.txt
```


## Workflow and Scripts

The general workflow looks like the following, starting with a KEGG pathway
module and the HMM profiles for orthologs for that module

  * Detect Proteins using ortholog profile HMMs: HMM search plus refinement
  * Classify Proteins: hmmscan
  * Cluster Proteins: MMSeqs2
  * Generate MSAs: Muscle
  * Finish Protein Sequences

The following instructions use "m00009", which is the KEGG module M00009
describing the TCA cycle. Replace this number with other module IDs as
appropriate.


-### Generating Query .hmm for a KEGG Module

Use `hmmfetch` to create a smaller HMM database for the KOs of a module. Use
the following naming convention but change the module ID: `data/m00009_ko.hmm`.
This smaller HMM is used during protein detection. Full KO HMM dataset is used
during classification.


### Detect Orthologs from Genomes

The following script puts outputs in `data/m00009_results` directory

```
./scripts/detect/search-genome m00009 GCF_002042975.1
```

Or if you have a list of genome accessions in a file, do the following. Note,
these commands append to existing files, so to re-run detection on a genome,
remove `data/m00009_results/proteins.{faa/tsv}` first.

```
./scripts/detect/search-genomes m00009 data/genomes_detect.txt
```

Use the following script to compare, for a given HMM model, how NCBI annotated
proteins (i.e. in `protein.faa` and `genomic.gff`) compare against protein
found by Needle.

```
PYTHONPATH=. python3 scripts/detect/compare-gff-with-match.py \
  --best-hmm \
  kegg-downloads/ko.hmm GCF_002042975.1 data/m00009_results/protein_detected.tsv \
  --output-file <filename>
```

### Classify Proteins

The following two commands will classify detected proteins first by KEGG
ortholog, then Pfam domains.

```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --cpu 4 --disable-cutoff-ga kegg-downloads/ko.hmm m00009
PYTHONPATH=. python3 scripts/classify/classify.py --cpu 4 pfam-downloads/Pfam-A.hmm m00009
```

Classification outputs appear in `data/m00009_results/classify.tsv`. If you
are re-running classification, remove this file first.

Annotated proteins submitted to NBCI can be classified in the same way, and
added to the same output TSV, using the following two commands.

```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --disable-cutoff-ga \
  --genome-accession GCF_932526225.1 \
  kegg-downloads/ko.hmm m00009
PYTHONPATH=. python3 scripts/classify/classify.py \
  --filter-by-prev-output \
  --genome-accession GCF_932526225.1 \
  pfam-downloads/Pfam-A.hmm m00009
```

The `--filter-by-prev-output` argument first filters the curated proteins to
remove those that do not appear in the `data/m00009_results/classify.tsv` file;
only those proteins matching one or more KEGG orthologs are further classified
using Pfam.

The following helper script, `classify-ncbi`, calls the above two commands for
each accession in an accession file.

```
scripts/classify/classify-ncbi data/genomes_ref.txt
```

Use the following script to generate a `protein_ncbi.tsv` and
`protein_names.tsv` files. Both include just proteins from the NCBI reference
genomes. `protein_ncbi.tsv` is similar to `protein_detected.tsv` and enumerates
exons.  `protein_names.tsv` lists the curated names of proteins, and is used in
the Tableau workbook.

```
PYTHONPATH=. python3 scripts/classify/generate-ref-protein-tsv.py m00009 \
  data/m00009_results/protein_ncbi.tsv \
  data/m00009_results/protein_names.tsv --overwrite
```

Use the `--overwrite` option to remove old data and re-generate these files.
Otherwise data will be added to the same files, including previously generated
names, causing duplication.

To add new names from a set of new reference genomes, without removing old
entries

```
PYTHONPATH=. python3 scripts/classify/generate-ref-protein-tsv.py m00009 \
  data/m00009_results/protein_ncbi.tsv \
  data/m00009_results/protein_names.tsv --genome-accession data/genomes_new.txt
```

*IMPORTANT* The above script fails for some GFFs that are malformed, e.g.
GCF_000001735.4 (Arabidopsis) that includes weird trans-splicing in the GenBank
file which GFF does not support. In these cases, manually fixing the GFFs to
work around the errors is the best option at the moment.

Use the following script to create FASTA files for orthologs, and domains for
each ortholog, based on classification results. The FASTA files are in
`data/m00009_results/faa` directory.

```
rm data/m00009_results/faa/*
PYTHONPATH=. python3 scripts/classify/assign.py m00009
```

Note that different scoring threshold criterias are used for detected proteins
(more tolerant) vs those from reference genomes (more stringent).


### Clustering (Optional)

For each KO, run the following script to cluster assigned sequences further

```
# make sure Docker daemon is running
./scripts/cluster/cluster m00009
```

Cluster outputs are summarized in `data/m00009_results/cluster.tsv`, and
clustered FAA files are in `data/m00009_results/clusters`.

Classification and clustering results -- i.e. how detected proteins match
against KO HMM profiles and how Pfam domains map onto those proteins assigned
to a KO -- can be visualized using Tableau. A template workbook that uses the
classification output TSV and several downloaded data files (e.g.
`genomes.tsv`, `ko.tsv`, and `Pfam-A.clans.tsv`), is `data/Protein
Classification.twb`.


### Generating Multi-Sequence Alignments

To generate MSAs and PNGs that visualize the MSAs, run the following script.
The `faa_dir` argument can be either the `data/m00009_results/faa` dir, or the
`data/m00009_results/clusters` dir.

```
# make sure Docker daemon is running
./scripts/align/generate-msas m00009 <faa_dir>
```

This script creates the `data/m00009_results/alignments` dir and, for each
input FAA file, generates a MSA FAA file, a PNG visualizing the MSA, and a HMM
profile from the MSA.


### Other Scripts

Search in SwissProt for related proteins

```
scripts/mmseqs-swissprot-search proteins.faa results.tsv
```

Download files from NCBI

```
PYTHONPATH=. python3 scripts/ncbi-download.py GCF_932526225.1
```
