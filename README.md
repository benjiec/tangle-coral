# Tangle: Coral

This repository includes instructions and scripts for organizing coral genomic,
transcriptomic, and proteomic data. The repository uses the other "tools"
centric repositories in the tangle project group

  * `tangle`: data models, scripts for downloading public datasets
  * `needle`: HMM based protein detection
  * `heap`: HMM and structural based classifications
  * `pile`: Transcriptomic data pipeline


## Tool Setup

Clone and setup the virtual environments for each of the above repositories, in
addition to this one. Follow the README in each repository. Some repositories
require building Docker images for use locally or on Google Cloud.

For this repository, use

```
python3 -m venv venv-coral
source venv-coral/bin/activate
pip3 install -r requirements.txt
```

Organize the virtual environments like the following (assumption going forward
in this document).

```
<tangle working directory>
  └── tangle/ 
  └── needle/ 
  └── heap/ 
  └── pile/ 
  └── coral/
  └── venv-tangle/ 
  └── venv-needle/ 
  └── venv-heap/ 
  └── venv-coral/ 
``` 

Then create the following aliases for command line execution

```
alias tangle-py='venv-tangle/bin/python3'
alias needle-py='venv-needle/bin/python3'
alias heap-py='venv-heap/bin/python3'
alias coral-py='venv-coral/bin/python3'
```

Unless otherwise specified, examples in this file assume the working directory
is the current directory.


## Data Setup

See `tangle/README.md`. Setup $TANGLE_WORLD and $TANGLE_AREA environment
variables. For the coral area, the last shoudld be set to "coral".

Files from `coral/data/` should be copied/mirroed to
`$TANGLE_WORLD/areas/$TANGLE_AREA/`.


```
tangle-py tangle/scripts/area/genome-list.py | tangle-py tangle/scripts/world/ncbi-download.py -
tangle-py tangle/scripts/area/genome-list.py | tangle-py tangle/scripts/world/ncbi-genome-metadata.py -
tangle-py tangle/scripts/world/kegg-download.sh
tangle-py tangle/scripts/world/pfam-download.sh
```

HMM profiles for KEGG and Pfam should be downloaded to directories specified by
matching environment variables, according to `heap/README.md`.


## Workflow and Scripts

### Needle: protein detection

Use the following command to generate a pooled .fna file, containing contigs
from all genomes requiring protein detection.

```
rm pooled.fna
tangle-py tangle/scripts/area/genome-list.py -d | \
  tangle-py tangle/scripts/defaults.py -f ncbi_genome_fna - | \
  xargs venv-tangle/bin/python3 tangle/scripts/pool-contigs.py pooled.fna
```

Then use that pooled file to setup a Google Cloud job

```
needle-py needle/gcloud/hmm-detect/setup.py \
  --genome-accession _ --run-dir-parent runs pooled.fna
rm pooled.fna
```

### Heap: classification and clustering

Filtering

```
heap-py heap/scripts/ko-filter-target.py runs/20260402_a611f70c/test.tsv x.filtered.tsv
```

Demux, which will produce protein files per genome after filtering

```
tangle-py tangle/scripts/demux-outputs.py test.tsv test.faa test_demux.tsv detected
```

TODO

Merge with curated proteins and pool accessions - use default script to output ncbi or detected, but not return if doesn't exist

Classify by KO

Classify by Pfam

Demux again

Assign KO, separate assignment table - before or after demux? if after, for each genome?

Cluster assigned, in cluster TSV, with name

Cluster putative from classify TSV, with name

Visualize feature projection of putative against KO and Pfam, by cluster

Dynamically compute cluster FAA to create a muscle alignment, can use a cluster TSV and a specific cluster
