# Tangle: Coral

This repository includes instructions and scripts for organizing coral genomic,
transcriptomic, and proteomic data. The repository uses the other "tools"
centric repositories in the tangle project group

  * `tangle`: data models, scripts for downloading public datasets
  * `needle`: HMM based protein detection
  * `heap`: HMM and structural based classifications
  * `pile`: Transcriptomic data pipeline


## BigQuery Tables

The following tables are loaded into BigQuery.

  * Public data
    * `kegg_modules`: names of modules, join with `module_id`
    * `kegg_orthologs`: names of KOs, join with `ortholog_id`
    * `kegg_module_definitions`: step by step specification of modules, join with `ortholog_id` and `module_id`
    * `pfam_domains`: names of Pfam families, join with `pfam_accession` and use `pfam_description`
    * `pfam_go`: mapping of Pfam families to GO terms, join with `pfam_accession` and use `go_description`

  * Genomic data
    * `genomes`: list of genome accessions and taxonomy, join with `Genome Accession`
    * `genomic_sequences`: protein manifest, join with `sequence_accession` *and* `sequence_database`
      * For genomic data (as oppose to RNAseq) `sequence_database` refers to a genome accession
    * `global_detected`: detected features with different semantics based on content
      * warning: there may be sequences in this table that do not appear in the `genomic_sequences` table
      * `target_database` is "KO": maps protein sequences (`query_accession` and `query_database`) to KEGG orthologs (`target_accession`)
      * `target_database` is "Pfam-A": maps protein sequences (`query_accession` and `query_database`) to Pfam families (`LEFT(target_accession,7)`)
      * other intended uses
        * exons mapping: `query_database` and `target_database` are genome accessions, `query_database` and `query_accession` identify a contig,
          and `target_database` and `target_accession` identify a protein sequence.
        * structural detection: `target_database` is `afdb_swissprot`, `query_database` and `query_accession` identify a protein sequence

   * Experiment data
     * `experiment_sequences`: list of sequences, join with `sequence_id`, filter by `sequence_database` as experiment
     * `experiment_detected`: detected features, much like `global_detected`
       * Join `sequence_id` with `query_accession`, `experiment_id` or `sequence_database` with `query_database`
     * `experiment_transcript_counts`: RNAseq transcript counts, join with `sequence_id` and `experiment_id`
     * `experiment_deseq2_tall`: tall table of DESeq2 analysis, join with `sequence_id` and `experiment_id`
       * `analysis_type` column describes base and test groups


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

Download `uniprot_sprot.dat.gz` from UniProt FTP, then run

```
tangle-py tangle/scripts/world/parse-uniprot-dat.py uniprot_sprot.dat.gz
```

to generate uniprot TSVs, there are two.


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
differences in their sequences.

Detection was run on many genomes. The following script demultiplex the results
and creates the appropriate detection tsv and protein fastas for each genome.

IMPORTANT: this script appends to each genome's TSV and fasta files, so REMOVE
PREVIOUS DATA if re-running.

```
# CANNOT run concurrently on different inputs because aggregates data by genome
tangle-py tangle/scripts/demux-outputs.py \
  --forget-original \
  --set-batch \
  --pooled-target-fasta-suffix .faa \
  --demuxed-parent-dir `tangle-py tangle/scripts/defaults.py -m area_genomics_dir` \
  runs/20260402_a611f70c/output_*.tsv
```

Use the following script to cleanup each demuxed tsv and protein fasta file
further, to remove entries contained in other entries at the same locus. This
process takes about 2-4 hours.

```
# can run concurrently on different inputs
needle-py needle/scripts/remove-contained.py \
  --forget-original \
  `tangle-py tangle/scripts/defaults.py -m area_genomics_dir`/*
```

Similar but not the same proteins at a locus are not being consolidated prior
to classification; while they increase classification compute time, their
presence may result in a more closely matched sequence assigned to a profile
post classification.


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
# can run concurrently on different inputs
tangle-py tangle/scripts/demux-outputs.py \
  --set-batch \
  --forget-original \
  runs/<run_dir>/outputs/sequence_ko_*
```

If you ran the `hmmscan-ko` job on Google Cloud, the results are already
filtered to remove those far below the KO HMM thresholds. Then you can just do

```
cat runs/<run_dir>/outputs/sequence_ko_*.tsv > sequence_ko_full.tsv
{ head -1 sequence_ko_full.tsv; grep -v query_database sequence_ko_full.tsv; } > sequence_ko.tsv; rm sequence_ko_full.tsv
mv sequence_ko.tsv `tangle-py tangle/scripts/defaults.py -m area_protein_ko_assigned_tsv`
```

However, if results were not filtered, run the following script

```
# CANNOT run concurrently on different inputs because aggregates outcome into one file
heap-py heap/scripts/ko-assign.py \
  --scoring-ratio-min 0.8  \
  `tangle-py tangle/scripts/defaults.py -m area_protein_ko_assigned_tsv` \
  runs/<run_dir>/outputs/sequence_ko_*.tsv
```

Then use the following script to filter each genome's protein tsv and fasta to
those that appear in the assignment file. Only do this for the detected
genomes.

```
# can run concurrently on different inputs
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
# can run concurrently on different inputs
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
# can run concurrently on different inputs
tangle-py tangle/scripts/demux-outputs.py \
  --set-batch \
  --forget-original \
  runs/<run_dir>/sequence_pfam_*.tsv

cat <run_dir>/sequence_pfam_*.tsv > sequence_pfam_full.tsv
{ head -1 sequence_pfam_full.tsv; grep -v query_database sequence_pfam_full.tsv; } > sequence_pfam.tsv; rm sequence_pfam_full.tsv
mv sequence_pfam.tsv `tangle-py tangle/scripts/defaults.py -m area_protein_pfam_tsv`
```


## Loading to BigQuery

Create dataset like this, with "tangle" as the dataset name

```
bq mk --location=[REGION] tangle_coral
```

### Generic assets

Load assets using

```
tangle-py tangle/scripts/bq-load-assets.py --dataset-name tangle_coral
```

Load UniProt to Pfam and GO mappings

```
tangle-py tangle/scripts/bq-schema.py \
  --check $TANGLE_WORLD/tangle/uniprot_pfam.tsv.gz \
  tangle.detected

tangle-py tangle/scripts/bq-schema.py \
  --check $TANGLE_WORLD/tangle/uniprot_go.tsv.gz \
  --table-name UniPGoTable \
  tangle.uniprot

tangle-py tangle/scripts/bq-schema.py \
  --check $TANGLE_WORLD/tangle/uniprot.tsv.gz \
  --table-name UniPTable \
  tangle.uniprot

# make sure the above runs successfully

tangle-py tangle/scripts/bq-schema.py \
  tangle.detected > tangle_detected.schema.json

tangle-py tangle/scripts/bq-schema.py \
  --table-name UniPGoTable \
  tangle.uniprot > tangle_uniprot_go.schema.json

tangle-py tangle/scripts/bq-schema.py \
  --table-name UniPTable \
  tangle.uniprot > tangle_uniprot.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.global_detected \
  $TANGLE_WORLD/tangle/uniprot_pfam.tsv.gz \
  ./tangle_detected.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.uniprot_go \
  $TANGLE_WORLD/tangle/uniprot_go.tsv.gz \
  ./tangle_uniprot_go.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.uniprot \
  $TANGLE_WORLD/tangle/uniprot.tsv.gz \
  ./tangle_uniprot.schema.json

rm ./tangle_detected.schema.json
rm ./tangle_uniprot_go.schema.json
rm ./tangle_uniprot.schema.json
```

### Genomics: KO and Pfam mappings

Load sequence KO assignments

```
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_protein_ko_assigned_tsv` \
  tangle.detected
# make sure the above runs successfully

tangle-py tangle/scripts/bq-schema.py \
  tangle.detected > tangle_detected.schema.json
bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.global_detected \
  `tangle-py tangle/scripts/defaults.py -m area_protein_ko_assigned_tsv` \
  ./tangle_detected.schema.json
rm ./tangle_detected.schema.json
```

Load sequence Pfam assignments

```
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_protein_pfam_tsv` \
  tangle.detected
# make sure the above runs successfully

tangle-py tangle/scripts/bq-schema.py \
  tangle.detected > tangle_detected.schema.json
bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.global_detected \
  `tangle-py tangle/scripts/defaults.py -m area_protein_pfam_tsv` \
  ./tangle_detected.schema.json
rm ./tangle_detected.schema.json
```

### Genomics: detected protein to contig mappings

**SKIP for now, as we don't have similar mappings for proteins from NCBI** The
latter requires parsing GFFs, and there are some malformed GFFs in NCBI that
makes parsing less automated.

Load genomic contig to protein mappings, from all the genomes

```
cat `tangle-py tangle/scripts/defaults.py -m area_genomics_dir`/GC*/proteins.tsv > proteins.tsv.all
{ head -1 proteins.tsv.all; grep -v query_database proteins.tsv.all; } > proteins.tsv; rm proteins.tsv.all

tangle-py tangle/scripts/bq-schema.py \
  --check proteins.tsv \
  tangle.detected
# make sure the above runs successfully

tangle-py tangle/scripts/bq-schema.py \
  tangle.detected > tangle_detected.schema.json
bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.global_detected \
  proteins.tsv \
  ./tangle_detected.schema.json
rm ./tangle_detected.schema.json
```

### Sequence and genome manifests

Use the following to load the manifest of all protein identifiers

```
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_sequence_manifest_tsv` \
  tangle.manifest
# make sure the above runs successfully

tangle-py tangle/scripts/bq-schema.py \
  tangle.manifest > tangle_manifest.schema.json
bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.genomic_sequences \
  `tangle-py tangle/scripts/defaults.py -m area_sequence_manifest_tsv` \
  ./tangle_manifest.schema.json
rm ./tangle_manifest.schema.json
```

And the following to load list of genomes

```
bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  --autodetect \
  tangle_coral.genomes \
  `tangle-py tangle/scripts/defaults.py -m area_genome_taxon_tsv`
```

### Tables for RNAseq experiments

Each experiment directory in area experiment dir (i.e. `tangle-py
tangle/scripts/defaults.py -m area_experiments_dir`) should have the following
files

  * `sequence_list.tsv`
  * `sequence_ko.tsv`
  * `sequence_pfam.tsv`
  * `sequence_data.tsv`
  * `des2_tall.tsv`

If the `sequence_ko.tsv` and `sequence_pfam.tsv` files are in the older format, use the following

```
coral-py coral/experiments/helpers/update-detected.py \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_ko.tsv EXP_PM41342399 KO
coral-py coral/experiments/helpers/update-detected.py \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_pfam.tsv EXP_PM41342399 Pfam-A --pfam-mode
```

Use the following script to standardize the column orders for
`sequence_data.tsv` and `des2_tall.tsv`

```
coral-py coral/experiments/helpers/update-exp.py \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399` EXP_PM41342399
```

Then, use the following to load into BigQuery, repeating for each experiment

```
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_list.tsv \
  tangle.manifest
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_ko.tsv \
  tangle.detected
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_pfam.tsv \
  tangle.detected
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_data.tsv \
  --table-name TranscriptCountsTable \
  tangle.exp
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/des2_tall.tsv \
  --table-name DESeq2Table \
  tangle.exp

# make sure the above runs successfully

tangle-py tangle/scripts/bq-schema.py tangle.manifest > tangle_manifest.schema.json
tangle-py tangle/scripts/bq-schema.py tangle.detected > tangle_detected.schema.json
tangle-py tangle/scripts/bq-schema.py --table-name TranscriptCountsTable tangle.exp > exp_transcript_counts.schema.json
tangle-py tangle/scripts/bq-schema.py --table-name DESeq2Table tangle.exp > exp_deseq2_tall.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_sequences \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_list.tsv \
  ./tangle_manifest.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_detected \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_ko.tsv \
  ./tangle_detected.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_detected \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_pfam.tsv \
  ./tangle_detected.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_transcript_counts \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/sequence_data.tsv \
  ./exp_transcript_counts.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_deseq2_tall \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM41342399`/des2_tall.tsv \
  ./exp_deseq2_tall.schema.json

rm ./tangle_manifest.schema.json
rm ./tangle_detected.schema.json
rm ./exp_transcript_counts.schema.json
rm ./exp_deseq2_tall.schema.json
```

To just load a file into the `experiment_detected` table (rows can be removed
with the unique `batch` value later) do the following, replacing
"sequence_fs.tsv" with your TSV file

```
tangle-py tangle/scripts/bq-schema.py \
  --check sequence_fs.tsv \
  tangle.detected

# make sure the above runs completely

tangle-py tangle/scripts/bq-schema.py tangle.detected > tangle_detected.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_detected \
  sequence_fs.tsv \
  ./tangle_detected.schema.json
```
