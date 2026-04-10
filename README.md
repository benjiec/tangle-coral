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

Then create the following aliases for command line execution (substitute full
paths)

```
alias tangle-py='venv-tangle/bin/python3'
alias needle-py='venv-needle/bin/python3'
alias heap-py='venv-heap/bin/python3'
alias coral-py='venv-coral/bin/python3'
alias pile-py='docker-compose -f pile/docker-compose.yml run --rm pile python3'
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


## RNAseq Workflows

Experiment specific scripts, to associate quants and conditions, are in
`experiments/` directory.

See [experiments/README.md](experiments/README.md) for more instructions on
processing quants, classifying transcripts, and more.


## Genomics Workflows

The goal of the genomics workflows is to end up with a picture of how many KOs
and KEGG modules are detectable in various genomes. For genomes without
predicted/curated proteins, a HMM based workflow detects and finds proteins
matching KO HMM profiles. For genomes with curated proteins, all curated
protein sequences are classified against KO HMM profiles.


### Needle: protein detection

Use the following command to generate a pooled .fna file, containing contigs
from all genomes requiring protein detection.

```
rm pooled.fna
tangle-py tangle/scripts/area/genome-list.py -d | \
  tangle-py tangle/scripts/defaults.py -f -m ncbi_genome_fna - | \
  xargs venv-tangle/bin/python3 tangle/scripts/pool-contigs.py pooled.fna
```

Then use that pooled file to setup a Google Cloud job

```
needle-py needle/gcloud/hmm-detect/setup.py \
  --genome-accession _ --run-dir-parent runs pooled.fna
rm pooled.fna
```


### Heap: classification

#### Preparing data for classification

Previous step uses HMM profiles to detect amino acid sequences and piece
sequences together into potential proteins. Many HMM profiles share the same or
very similar domains, thus resulting in multiple detected proteins with slight
differences in their sequences. Also, some detected proteins may only contain
enough amino acid sequences to partially match to a profile, and have no chance
to ever being classified as a full protein downstream.

A couple of the scripts below attempt to remove some of these junk or
duplications. Similar but not the same proteins at a locus are not being
consolidated prior to classification; while they increase classification
compute time, their presence may result in a more closely matched sequence
assigned to a profile post classification.

The following script filters them away small matches that have no chance of
being classified and assigned to a KO. Remove the `--forget-original` argument
to leave a copy of the original tsv with `.orig` suffix.

```
heap-py heap/scripts/ko-filter-target.py \
  --forget-original \
  runs/20260402_a611f70c/output_*.tsv
```

Detection was run on many genomes. The following script demultiplex the results
and creates the appropriate detection tsv and protein fastas for each genome.
Note that this script appends to each genome's TSV and fasta files, so remove
previous copies if re-running.

```
tangle-py tangle/scripts/demux-outputs.py \
  --forget-original \
  --pooled-target-fasta-suffix .faa \
  --demuxed-parent-dir `tangle-py tangle/scripts/defaults.py -m area_genomics_dir` \
  runs/20260402_a611f70c/output_*.tsv
```

Use the following script to cleanup each demuxed tsv and protein fasta file
further, to remove entries contained in other entries at the same locus. This
process takes about 2-4 hours.

```
needle-py needle/scripts/remove-contained.py \
  --forget-original \
  `tangle-py tangle/scripts/defaults.py -m area_genomics_dir`/*
```


#### Classification against KO

Use the following to generate a pooled proteins FASTA file. Likely, there are
~20M sequences.

```
tangle-py tangle/scripts/area/genome-list.py | \
  tangle-py tangle/scripts/defaults.py \
  -m area_detected_proteins \
  -m ncbi_genome_proteins \
  -f - | \
  xargs venv-tangle/bin/python3 tangle/scripts/pool-contigs.py \
  pooled_proteins.faa
```

Use the `pooled_proteins.faa` to setup a KO classification job on Google Cloud

```
heap-py heap/gcloud/hmmscan-ko/setup.py \
  --run-dir-parent runs \
  --query-database-name _ pooled_proteins.faa
```

Use the rclone option to download individual output files into an outputs
directory. Run the following to demultiplex the results.

```
tangle-py tangle/scripts/demux-outputs.py \
  --forget-original \
  runs/<run_dir>/outputs/sequence_ko_*
```

Use the following script to select proteins that match or are close to matching
a KO.

```
heap-py heap/scripts/ko-assign.py \
  --scoring-ratio-min 0.8  \
  `tangle-py tangle/scripts/defaults.py -m area_protein_ko_assigned_tsv` \
  runs/<run_dir>/sequence_ko_*.tsv
```

Then use the following script to filter each genome's protein tsv and fasta to
those that appear in the assignment file. Only do this for the detected
genomes.

```
tangle-py tangle/scripts/filter-proteins.py \
  --filter-targets-with-queries-from \
    `tangle-py tangle/scripts/defaults.py -m area_protein_ko_assigned_tsv` \
  `tangle-py tangle/scripts/defaults.py -m area_genomics_dir`/GC*
```

There are still a good number of duplicate proteins in each detected proteome,
at this point, because each HMM profile in the KO database was used to detect
proteins separately and several could have found similar proteins at the same
locus. Use the following script to further filter. Warning, this is a monster
of a script, highly recommend keeping the .orig files in case something goes
wrong.

```
tangle-py tangle/scripts/filter-clustered-proteins.py \
  --ko-classification-tsv \
    `tangle-py tangle/scripts/defaults.py -m area_protein_ko_assigned_tsv` \
  `tangle-py tangle/scripts/defaults.py -m area_genomics_dir`/GC*
```


#### Generating a manifest

Use the following two commands to create an manifest of both the detected and
NCBI protein sequences.

```
tangle-py tangle/scripts/area/genome-list.py | \
  tangle-py tangle/scripts/defaults.py \
  -m area_detected_proteins -f - | \
  xargs venv-tangle/bin/python3 tangle/scripts/manifest.py \
    --sequence-source hmm-detected \
    --sequence-type protein \
    `tangle-py tangle/scripts/defaults.py -m area_sequence_manifest_tsv`

tangle-py tangle/scripts/area/genome-list.py | \
  tangle-py tangle/scripts/defaults.py \
  -m ncbi_genome_proteins -f - | \
  xargs venv-tangle/bin/python3 tangle/scripts/manifest.py \
    --append \
    --sequence-source ncbi \
    --sequence-type protein \
    `tangle-py tangle/scripts/defaults.py -m area_sequence_manifest_tsv`
```


#### Classification against Pfam

Use the following to generate a pooled proteins FASTA file, now much smaller
than the input used for KO classification.

```
tangle-py tangle/scripts/area/genome-list.py | \
  tangle-py tangle/scripts/defaults.py \
  -m area_detected_proteins \
  -m ncbi_genome_proteins \
  -f - | \
  xargs venv-tangle/bin/python3 tangle/scripts/pool-contigs.py \
  pooled_proteins.faa
```

Submit a HMM scan job on Google Cloud

```
heap-py heap/gcloud/hmmscan-pfam/setup.py \
  --run-dir-parent runs \
  --query-database-name _ pooled_proteins.faa
```

Use the rclone option to download individual output files into an outputs
directory, then use the following to demultiplex the outputs and concatenate.

```
tangle-py tangle/scripts/demux-outputs.py \
  --forget-original \
  runs/<run_dir>/sequence_pfam_*.tsv

cat <run_dir>/sequence_pfam_*.tsv > sequence_pfam_full.tsv
{ head -1 sequence_pfam_full.tsv; grep -v query_database sequence_pfam_full.tsv; } > sequence_pfam.tsv; rm sequence_pfam_full.tsv
mv sequence_pfam.tsv `tangle-py tangle/scripts/defaults.py -m area_protein_pfam_tsv`
```


### MMSeqs: Clustering

Cluster assigned, in cluster TSV, with name

Cluster putative from classify TSV, with name

Visualize feature projection of putative against KO and Pfam, by cluster
