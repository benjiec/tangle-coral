# Tangle: Coral

This repository keeps instructions and scripts for populating the Tangle area
"Coral". It is intended to be used with the following repositories

  * tangle-tangle: data models and scripts for downloading public information
  * tangle-needle: HMM based protein detection
  * tangle-heap: HMM and structural based classifications
  * tangle-pile: Transcriptomic data pipeline


## Setup

In addition to the above repositories, which should all be cloned and setup
according the the README files in those repositories, get the following

  * MMSeqs2 docker image: `docker pull ghcr.io/soedinglab/mmseqs2`
  * Muscle aligner: `docker pull pegi3s/muscle`

Also create a Python virtualenv and then install required packages.

```
python3 -m venv venv-coral
source venv-coral/bin/activate
pip3 install -r requirements.txt
```


## Data Setup

Setup $TANGLE_WORLD and $TANGLE_AREA environment variables.

Files from `data/` should be copied/mirroed to `$TANGLE_WORLD/areas/$TANGLE_AREA/`.

### Tangle

```
python3 scripts/area/genome-list.py | python3 scripts/world/ncbi-download.py -
python3 scripts/area/genome-list.py | python3 scripts/world/ncbi-genome-metadata.py -
scripts/world/kegg-download.sh
```

### Heap

HMM profiles for KEGG and Pfam should be downloaded to the approriate
directories according to README.md file.


### Other (move to Tangle repo soon)

Download this

```
curl https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz -o $TANGLE_WORLD/tangle/Pfam-A.clans.tsv.gz
gunzip $TANGLE_WORLD/tangle/Pfam-A.clans.tsv.gz
```

Where did this come from `pfam_go_simple.tsv`?

See `heap/helpers/ipr-pfam.py` on how to generate an UniProt to Pfam TSV, as a
DetectedTable. This file should be in
`$TANGLE_WORLD/tangle/uniprot-pfam.tsv.gz`.


## Workflow and Scripts

### Needle

Prepare genome FAA file for detection, by concatenating genome accession with
contigs - TODO need a script

Use needle repo to setup a GC batch job to detect proteins

Use needle repo `generate-ref-protein-tsv.py` script to parse GFF into detected
list - TODO need a script

### Heap

Use heap repo to filter detected FAA to remove likely unclassifiable results -
TODO need a script

Use heap repo to setup a GC batch job detect KO and detect Pfam, from both
detected proteins and reference proteins - TODO need a script

Use heap repo to assign KOs

Use heap repo to cluster proteins - need a script to generate FAAs then
cluster.

Use heap repo `tabularize-alignment.py` script to generate TSV of protein
alignments with KO and Pfam mapped on.
