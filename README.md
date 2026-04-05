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

Setup $TANGLE_WORLD and $TANGLE_AREA environment variables. The last shoudld be
set to "coral", for example.

Files from `./data/` should be copied/mirroed to `$TANGLE_WORLD/areas/$TANGLE_AREA/`.

### From `tango` repository

```
python3 scripts/area/genome-list.py | python3 scripts/world/ncbi-download.py -
python3 scripts/area/genome-list.py | python3 scripts/world/ncbi-genome-metadata.py -
scripts/world/kegg-download.sh
scripts/world/pfam-download.sh
```

### From `heap` repository

HMM profiles for KEGG and Pfam should be downloaded to the approriate
directories according to README.md file.


## Workflow and Scripts

### Needle

Pooled genomic accessions in setup script

### Heap

Heap helper: Detected protein filtered to likely protein

Merge with curated proteins and pool accessions

Classify by KO

Classify by Pfam

Assign KO, separate assignment table

Cluster assigned, in cluster TSV, with name

Cluster putative from classify TSV, with name

Visualize feature projection of putative against KO and Pfam, by cluster

Dynamically compute cluster FAA to create a muscle alignment, can use a cluster TSV and a specific cluster
