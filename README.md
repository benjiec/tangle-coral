# Needle

Needle is a suite of tools to curate and compare sequences of proteins by
pathway, for understudied organisms. Cryptic gene structures, unconventional
splicing models, unfinished genomes, understudied proteomes, are some of the
challenges that make applying traditional workflows (e.g. gene prediction and
classification) difficult and/or unreliable. Needle uses existing HMM profiles
to detect and classify proteins on raw genomic DNA or predicted proteomes, then
curate/cluster proteins in a phylogenetic aware manner by pathway.


## Setup

Install the following programs

  * hmmer3 package: e.g. on MacOS run `brew install hmmer`
  * MMSeqs2 docker image: `docker pull ghcr.io/soedinglab/mmseqs2`
  * Muscle aligner: `docker pull pegi3s/muscle`
  * (Optional) SwissProt DB for MMSeqs2: `scripts/data/mmseqs-swissprot-setup`
  * (Optional) NCBI docker image: `docker pull ncbi/blast`

Also create a Python virtualenv and then install required packages.

```
python3 -m venv .venv
source .venv/bin/activate
pip3 install -r requirements.txt
```

## Initial Data

There are some initial data already in `data` directory. Here are some
additional data to download.

Download

  * KEGG KO profile HMMs: `https://www.genome.jp/ftp/db/kofam/profiles.tar.gz`
    * Then concatenate all the profiles together: `cat profiles/*.hmm > kegg-downloads/ko.hmm`
    * Run `hmmpress kegg-downloads/ko.hmm`

  * KEGG ortholog table: `data/ko.tsv`

    ```
    curl https://rest.kegg.jp/list/ko -o data/ko.txt
    echo "Ortholog ID\tOrtholog Name" | cat - data/ko.txt > data/ko.tsv
    rm data/ko.txt
    # then manually remove `"` and replace them with `''` in this file
    ```

  * KEGG modules list: `data/modules.tsv` `data/module_defs.csv`

    ```
    curl https://rest.kegg.jp/list/module -o data/modules.txt
    echo "Module ID\tModule Name" | cat - data/modules.txt > data/modules.tsv
    rm data/modules.txt
    PYTHONPATH=. python3 scripts/data/fetch-kegg-module-def.py
    ```

  * Pfam HMM profiles: `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`
    * Store uncompress file at: `pfam-downloads/Pfam-A.hmm`
    * Run `hmmpress pfam-downloads/Pfam-A.hmm`
    * Also
      * `data/ko_thresholds.gz`: from KEGG FTP site `https://www.genome.jp/ftp/db/kofam/ko_list.gz`
      * `data/Pfam-A.clans.tsv`: also from above FTP site


## Prepare List of Genomes

Put NCBI genome accessions into two files

  * To run de novo HMM based protein detection: `data/genomes_detect.txt`
  * To process curated proteins submitted to NCBI: `data/genomes_ref.txt`

A genome can be in both files, but genomes in the second file must have a
`proteins.faa` file at NCBI available to download.

Run the following commands to generate `data/genomes.tsv`, which includes
genome name and taxonomy information.

```
PYTHONPATH=. python3 scripts/data/fetch-genomes.py data/genomes_detect.txt
PYTHONPATH=. python3 scripts/data/fetch-genomes.py data/genomes_ref.txt
```

Also download genomic sequences, using

```
PYTHONPATH=. python3 scripts/data/ncbi-download.py data/genomes_detect.txt
PYTHONPATH=. python3 scripts/data/ncbi-download.py data/genomes_ref.txt
```

Genome sequences are large. The above scripts store them in
`ncbi-downloads/ncbi-dataset`. You can symlink a different volume (e.g. an
external drive) to that path.


## Workflow and Scripts

The general workflow looks like the following

  * Create a HMM database for a KEGG module
  * Detect proteins using the module specific HMM database
  * Classify detected proteins using the full KEGG HMM database
  * Filter reference proteins (those from NCBI) by module specific HMM database then classify using the full KEGG HMM database
  * Assign proteins to KO numbers, and generate list of putative proteins lacking definitive assignment
  * Cluster assigned and putative proteins
  * Generate MSA for each cluster

The following instructions use "m00009", which is the KEGG module M00009
describing the TCA cycle. Replace this number with other module IDs as
appropriate.


### Create a HMM database for a KEGG module

Use `hmmfetch` to create a smaller HMM database for the KOs of a module. Use
the following naming convention but change the module ID: `data/m00009_ko.hmm`.
This smaller HMM is used during protein detection. Full KO HMM dataset is used
during classification.


### Detect proteins using the module specific HMM database

```
./scripts/detect/search-genomes m00009 data/genomes_detect.txt
```

Or if just for one genome accession

```
./scripts/detect/search-genome m00009 GCF_002042975.1
```

Note: these scripts append to existing files, so to re-run detection on a
genome, remove `data/m00009_results/proteins.{faa/tsv}` first.


### Classify detected proteins using the full KEGG HMM database

Classification outputs appear in `data/m00009_results/classify.tsv`. If you are
re-running classification, remove this file first.

The following command will classify detected proteins against full KEGG ortholog database

```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --cpu 4 --disable-cutoff-ga kegg-downloads/ko.hmm m00009
```


### Filter and classify reference proteins from NCBI

Annotated proteins submitted to NBCI should be first filtered to those relevant
to the KEGG module, then classified using the full KO HMM database. Otherwise
the classification process will run too long.

The following script does everything

```
scripts/classify/classify-ncbi m00009 data/genomes_ref.txt
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


### Assign proteins to KO numbers and curate various assets

Use the following script to create FASTA files with proteins assigned to
orthologs, as well as putative proteins, and some other TSV assets.

```
rm data/m00009_results/faa/*
PYTHONPATH=. python3 scripts/classify/assign.py m00009
```

Note that different scoring threshold criterias are used for detected proteins
(more tolerant) vs those from reference genomes (more stringent).

Afterward, use the following scripts to compute Pfam matches

```
rm data/m00009_results/candidate_pfam.tsv
PYTHONPATH=. python3 scripts/classify/classify.py \
  --cpu 4 --filter-by data/m00009_results/candidate_ko.tsv \
  pfam-downloads/Pfam-A.hmm m00009 data/m00009_results/candidate_pfam.tsv
scripts/classify/classify-pfam-ncbi m00009 data/genomes_ref.txt
```


### Cluster assigned and putative proteins

```
# make sure Docker daemon is running
./scripts/cluster/cluster m00009
```

Cluster outputs are summarized in `data/m00009_results/cluster.tsv`, and
clustered FAA files are in `data/m00009_results/clusters`.


### Generate MSA for each cluster

```
# make sure Docker daemon is running
./scripts/align/generate-msas m00009 data/m00009_results/clusters
```

### Using results

Classification and clustering results -- i.e. how detected proteins match
against KO HMM profiles and how Pfam domains map onto those proteins assigned
to a KO -- can be visualized using Tableau. Example workbooks, for m00009, are

  * `data/Protein Classification.twb`: joins several TSV files and provides
    global view into what species have what proteins for what pathway steps

  * `data/Alignment.twb`: visualizes protein alignments for clusters position
    by position, with domain and ortholog annotations.

    * This one is specific to a KO ID and requires a tabularized alignment output file

      ```
      PYTHONPATH=. python3 scripts/align/tabularize-alignment.py \
         m00009 K00164 data/m00009_results/alignments/K00164.tsv
      ```

Also `needle/duckdb.py` shows how to use duckdb to join and query the TSV
files.


### Other Scripts

Search in SwissProt for related proteins

```
scripts/mmseqs-swissprot-search proteins.faa results.tsv
```

Download files from NCBI

```
PYTHONPATH=. python3 scripts/ncbi-download.py GCF_932526225.1
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

KO to Pfam mapping was generated using KO consensus sequence. `--cut_ga` option
was used with `hmmscan`, in the following script; technically, all hits
satisfied the gathering criteria defined by Pfam. The table was manually
cleaned up into the format in `data/ko_pfam.tsv`.

```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --genome-accession _ --fasta-file kegg-downloads/ko.consensus.faa pfam-downloads/Pfam-A.hmm _ data/ko_pfam_raw.tsv
```
