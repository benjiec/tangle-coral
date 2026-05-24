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
    * `uniprot`: mapping of UniProt accession to function description, join with `uniprot_accession`
    * `uniprot_go`: mapping of UniProt accession to GO terms, join with `uniprot_accession`
    * `orthodb_uniprot_groups`: mapping of UniProt accessions to OrthoDB groups, join with `uniprot_accession`
    * `kegg_ko_thresholds`: KO HMM thresholds and consensus lenth, join ortholog ID with `model`

  * Genomic data
    * `genomes`: list of genome accessions and taxonomy, join with `Genome Accession`
    * `genomic_sequences`: protein manifest, join with `sequence_accession` *and* `sequence_database`
      * For genomic data (as oppose to RNAseq) `sequence_database` refers to a genome accession
    * `global_detected`: detected features with different semantics based on content
      * warning: there may be sequences in this table that do not appear in the `genomic_sequences` table
      * `target_database` is "KO": maps protein sequences (`query_accession` and `query_database`) to KEGG orthologs (`target_accession`)
      * `target_database` is "Pfam-A": maps protein sequences (`query_accession` and `query_database`) to Pfam families (`LEFT(target_accession,7)`)
      * `target_database` is "afdb_swissprot": maps protein sequences (`query_accession` and `query_database`) to SwissProt (`target_accession`)
      * other intended uses
        * exons mapping: `query_database` and `target_database` are genome accessions, `query_database` and `query_accession` identify a contig,
          and `target_database` and `target_accession` identify a protein sequence.

   * Experiment data
     * `experiment_transcripts`: list of transcripts, join with `sequence_id`, filter by `sequence_database` as experiment
     * `experiment_transcript_proteins`: join table of transcripts to proteins, see below on how to use
     * `experiment_detected`: protein domains or ortholog matches, much like `global_detected`
       * Join protein accessions with `query_accession`, `experiment_id` or `sequence_database` with `query_database`
     * `experiment_transcript_genes`: mapping transcript IDs to gene IDs that are then counted, see below
     * `experiment_gene_counts`: RNAseq transcript counts, join with `gene_id` and `experiment_id`
     * `experiment_deseq2_tall`: tall table of DESeq2 analysis, join with `gene_id` and `experiment_id`
       * `analysis_type` column describes base and test groups


### Experiment Transcript to Proteins Table

The `experiment_transcript_proteins` table uses the `tangle.DetectedTable`
format to specify how transcript IDs and protein IDs are related. Some of the
important columns for joining/analyzing data are

  * `query_database`, `target_database`: join with `experiment_id` from other tables
  * `query_accession` when `query_type` is "transcript": transcript ID or accession
  * `query_accession` when `query_type` is "cds": cds ID, not used outside of this table
  * `target_accession` when `target_type` is "cds": cds ID
  * `target_accession` when `target_type` is "protein": protein ID or accession, joinable with `experiment_detected` accessions

Transcript counts and differential expression statistics are based on
transcript IDs, in `experiment_transcripts`, `experiment_gene_counts`, and
`experiment_deseq2_tall` tables. `experiment_transcripts` table can join with
`experiment_transcript_proteins` twice, first to translate transcript ID to
cds ID, then from cds ID to a protein ID or accession. The outcome of this
join can then be further joined with `experiment_detected` to find Pfam domains
or KEGG ortholog matches.

It may be important to understand where a Pfam or KO match appears on the
transcript. The following terms may be more human readable, and can be computed
from the above tables and joins.

  * `transcript_start`, `transcript_end`
    * First `experiment_transcript_proteins` join, `query_start` and
      `query_end`. If a cds is on reverse strand, then the start value is
      larger than the end value.
  * `cds_start_on_transcript`, `cds_end_on_transcript`
    * First `experiment_transcript_proteins` join, `target_start` and
      `target_end`.
  * `protein_start_on_cds`, `protein_end_on_cds`
    * Second `experiment_transcript_proteins` join, `query_start` and
      `query_end`.
  * `protein_start`, `protein_end`
    * Second `experiment_transcript_proteins` join, `target_start` and
      `target_end`.
  * `protein_start_on_transcript`, `protein_end_on_transcript`: RNA bp location
    on transcript where protein starts and ends
    * See below for formulas
  * `protein_length`
    * This is `abs(protein_start-protein_end)+1`
  * `match_start_on_protein`, `match_end_on_protein`
    * `experiment_detected`, `query_start` and `query_end`
  * `matched_start`, `matched_end`
    * `experiment_detected`, `target_start` and `target_end`
  * `protein_match_length`
    * This is `abs(match_start_on_protein-match_end_on_protein)+1`
  * `matched_length`
    * This is `abs(matched_start-matched_end)+1`
  * `protein_match_percent`
    * This is `100*protein_match_length/protein_length`
  * `ko_consensus_length`
    * Left join `experiment_detected` with `kegg_ko_thresholds` on `target_database = "KO"` and `target_accession = kegg_ko_thresholds.model`
    * This is the `mlen` column in the `kegg_ko_thresholds` table
  * `matched_percent`
    * This is `100*matched_length/ko_consensus_length`

Formula for `protein_start_on_transcript`

```
IF [Transcript Start] < [Transcript End]
THEN
  [Protein Start on CDS]+([CDS Start on Transcript]-[Transcript Start])
ELSE
  [Protein End on CDS]+([CDS Start on Transcript]-[Transcript End])
END
```

Formula for `protein_end_on_transcript`

```
IF [Transcript Start] < [Transcript End]
THEN
  [Protein End on CDS]+([CDS Start on Transcript]-[Transcript Start])
ELSE
  [Protein Start on CDS]+([CDS Start on Transcript]-[Transcript End])
END
```

### Experiment Transcript to Genes Table

RNASeq experiments count and analyze differential expression by genes,
potentially aggregating multiple transcripts of gene isoforms together. While
each experiment's transcriptome is represented by a list of transcripts, the
counts and deseq2 tables refer to gene IDs. Use the following join

  * `sequence_id` and `sequence_database` columns from `experiment_transcripts`
  * Join `transcript_id` and `experiment_id` columns of `experiment_transcript_genes`
  * Join `gene_id` and `experiment_id` columns of `experiment_gene_counts` or `experiemnt_deseq2_tall`


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
tangle/scripts/world/kegg-download.py
tangle/scripts/world/pfam-download.sh
tangle/scripts/world/odb-download.sh
```

Parse the OrthoDB files into a TSV file that joins UniProt IDs to OrthoDB groups.

```
tangle-py tangle/scripts/world/parse-odb.py
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

Use the following commands to a) generate a pooled .fna file, containing
contigs from all genomes requiring protein detection, and b) setup a Google
Cloud job.

```
rm pooled.fna
tangle-py tangle/scripts/area/genome-list.py -d | \
  tangle-py tangle/scripts/defaults.py -f -m ncbi_genome_fna - | \
  xargs venv-tangle/bin/python3 tangle/scripts/pool-contigs.py pooled.fna

needle-py needle/gcloud/hmm-detect/setup.py \
  --genome-accession _ --run-dir-parent runs pooled.fna
rm pooled.fna
```

After run completes, use the following to demux the outputs by genome accession
and create .faa and .tsv files for each detected genome.

IMPORTANT: this script appends to each genome's TSV and fasta files, so REMOVE
PREVIOUS DATA if re-running.

```
# CANNOT run concurrently on different inputs because aggregates data by genome
tangle-py tangle/scripts/demux-outputs.py \
  --forget-original \
  --set-batch \
  --pooled-target-fasta-suffix .faa \
  --demuxed-parent-dir `tangle-py tangle/scripts/defaults.py -m area_genomics_dir` \
  runs/<run_dir>/outputs/output_*.tsv
```


### Heap: classification

#### Classification by KO

Use the following commands to a) generate a pooled proteins FASTA file, and b)
generate a Google Cloud job.

```
tangle-py tangle/scripts/area/genome-list.py | \
  tangle-py tangle/scripts/defaults.py \
  -m area_detected_proteins \
  -m ncbi_genome_proteins \
  -f - | \
  xargs venv-tangle/bin/python3 tangle/scripts/pool-contigs.py \
  pooled_proteins.faa

heap-py heap/gcloud/hmmscan-ko/setup.py \
  --run-dir-parent runs \
  --query-database-name _ pooled_proteins.faa
```

Run the following to demultiplex the results.

```
# can run concurrently on different inputs
tangle-py tangle/scripts/demux-outputs.py \
  --set-batch \
  --forget-original \
  runs/<run_dir>/outputs/sequence_ko_*

cat runs/<run_dir>/outputs/sequence_ko_*.tsv > sequence_ko_full.tsv
{ head -1 sequence_ko_full.tsv; grep -v query_database sequence_ko_full.tsv; } > sequence_ko.tsv; rm sequence_ko_full.tsv
mv sequence_ko.tsv `tangle-py tangle/scripts/defaults.py -m area_protein_ko_assigned_tsv`
```

Use the following script to filter each genome's protein tsv and fasta to those
that appear in the assignment file. Only do this for the detected genomes.

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

Use the following commands to a) generate a pooled proteins FASTA file, and b)
setup a Google Cloud job.

```
tangle-py tangle/scripts/area/genome-list.py | \
  tangle-py tangle/scripts/defaults.py \
  -m area_detected_proteins \
  -m ncbi_genome_proteins \
  -f - | \
  xargs venv-tangle/bin/python3 tangle/scripts/pool-contigs.py \
  pooled_proteins.faa

heap-py heap/gcloud/hmmscan-pfam/setup.py \
  --run-dir-parent runs \
  --query-database-name _ pooled_proteins.faa
```

Use the following to demultiplex the outputs and concatenate.

```
# can run concurrently on different inputs
tangle-py tangle/scripts/demux-outputs.py \
  --set-batch \
  --forget-original \
  runs/<run_dir>/outputs/sequence_pfam_*.tsv

cat runs/<run_dir>/outputs/sequence_pfam_*.tsv > sequence_pfam_full.tsv
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

Load KEGG KO thresholds

```
tangle-py tangle/scripts/bq-schema.py \
  --check ko_thresholds.tsv \
  tangle.kegg

# make sure the above runs successfully

tangle-py tangle/scripts/bq-schema.py \
  tangle.kegg > tangle_kegg_ko_thresholds.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.kegg_ko_thresholds \
  ko_thresholds.tsv \
  ./tangle_kegg_ko_thresholds.schema.json

rm ./tangle_kegg_ko_thresholds.schema.json
```

Load OrthoDB groups and join table with UniProt accessions

```
tangle-py tangle/scripts/bq-schema.py \
  --check $TANGLE_WORLD/tangle/odb_uniprot_groups.tsv.gz \
  tangle.orthodb

# make sure the above runs successfully

tangle-py tangle/scripts/bq-schema.py \
  tangle.orthodb > tangle_odb_uniprot_groups.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.orthodb_uniprot_groups \
  $TANGLE_WORLD/tangle/odb_uniprot_groups.tsv.gz \
  ./tangle_odb_uniprot_groups.schema.json

rm ./tangle_odb_uniprot_groups.schema.json
```


### Tables for genomics data

The following tables can be removed, if reloading data

  * `genomic_sequences`
  * `genomes` - if there are new genomes
  * `global_detected` - because SwissProt to Pfam data is in here as well, remove old data by `batch` value using Google Cloud console


#### KO and Pfam assignments

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

#### Detected protein to contig mappings

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

#### Sequence and genome manifests

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

  * `transcript_list.tsv`
  * `transcript_proteins.tsv`
  * `sequence_ko.tsv`
  * `sequence_pfam.tsv`
  * `transcript_genes.tsv`
  * `gene_counts.tsv`
  * `des2_tall.tsv`

Use the following to load into BigQuery, repeating for each experiment

```
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/transcript_list.tsv \
  tangle.manifest
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/transcript_proteins.tsv \
  tangle.detected
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_ko.tsv \
  tangle.detected
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_pfam.tsv \
  tangle.detected
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/gene_counts.tsv \
  --table-name GeneCountsTable \
  tangle.exp
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/transcript_genes.tsv \
  --table-name TranscriptGenesTable \
  tangle.exp
tangle-py tangle/scripts/bq-schema.py \
  --check `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/des2_tall.tsv \
  --table-name DESeq2Table \
  tangle.exp

# make sure the above runs successfully

tangle-py tangle/scripts/bq-schema.py tangle.manifest > tangle_manifest.schema.json
tangle-py tangle/scripts/bq-schema.py tangle.detected > tangle_detected.schema.json
tangle-py tangle/scripts/bq-schema.py --table-name GeneCountsTable tangle.exp > exp_gene_counts.schema.json
tangle-py tangle/scripts/bq-schema.py --table-name TranscriptGenesTable tangle.exp > exp_transcript_genes.schema.json
tangle-py tangle/scripts/bq-schema.py --table-name DESeq2Table tangle.exp > exp_deseq2_tall.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_transcripts \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/transcript_list.tsv \
  ./tangle_manifest.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_detected \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/transcript_proteins.tsv \
  ./tangle_detected.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_detected \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_ko.tsv \
  ./tangle_detected.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_detected \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/sequence_pfam.tsv \
  ./tangle_detected.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_gene_counts \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/gene_counts.tsv \
  ./exp_gene_counts.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_transcript_genes \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/transcript_genes.tsv \
  ./exp_transcript_genes.schema.json

bq load \
  --source_format=CSV \
  --field_delimiter='\t' \
  --skip_leading_rows=1 \
  tangle_coral.experiment_deseq2_tall \
  `tangle-py tangle/scripts/defaults.py -m area_experiment PM34593802`/des2_tall.tsv \
  ./exp_deseq2_tall.schema.json

rm ./tangle_manifest.schema.json
rm ./tangle_detected.schema.json
rm ./exp_gene_counts.schema.json
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
