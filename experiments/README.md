# Analysis Tasks

## Download data from NCBI SRA

To install tools

```
pip3 install pysradb
brew install awscli
```

Use the following to get a list of SRR accessions to download

```
pysradb metadata PRJNA591730 > metadata.tsv
```

Then put together a script to download them using

```
aws s3 cp s3://sra-pub-run-odp/sra/SRR36764133/SRR36764133 . --no-sign-request
```

To extract fastq files, use NCBI SRA toolkit (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)

```
fasterq-dump SRR9331965
```

## Running Salmon to map reads

Install Salmon using Conda into an environment, like "salmon_env"

Then run the following before you can start using salmon

```
conda init
source ~/.bash_profile
conda activate salmon_env
```

### Building decoy-aware inde

```
grep "^>" c_goreaui.fna | cut -d " " -f 1 | sed -e 's/>//g' > c_goreaui.decoys.txt
cat c_goreaui.transcripts.fna c_goreaui.fna > c_goreaui.gentrome.fna   # transcript MUST COME FIRST
salmon index -t c_goreaui.gentrome.fna -d c_goreaui.decoys.txt -i c_goreaui.salmon_index

grep "^>" aten.fna | cut -d " " -f 1 | sed -e 's/>//g' > aten.decoys.txt
cat aten.transcripts.fna aten.fna > aten.gentrome.fna
salmon index -t aten.gentrome.fna -d aten.decoys.txt -i aten.salmon_index
```

### Quant using Salmon

Change the "-p 2" to however many CPUs/cores you would like to use

```
salmon quant -i aten.salmon_index \
             -l A \
             -1 sra/SRR9331965_1.fastq \
             -2 sra/SRR9331965_2.fastq \
             --validateMappings \
             -o aten-quants/SRR9331965 -p 2

salmon quant -i c_goreaui.salmon_index \
             -l A \
             -1 sra/SRR9331965_1.fastq \
             -2 sra/SRR9331965_2.fastq \
             --validateMappings \
             -o c_goreaui-quants/SRR9331965 -p 2
```

### Process quants

Each directory may have its own script to process quants into `sequence_data.tsv`.

```
python3 experiments/doi:10.1126_sciadv.aba2498/process_salmon_quants.py \
  experiments/doi:10.1126_sciadv.aba2498 data/exp_results/doi:10.1126_sciadv.aba2498
python3 scripts/data/update-genome-accession.py \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv \
  GCA_014633955.1 --when doi:10.1126/sciadv.aba2498-a_tenuis
python3 scripts/data/update-genome-accession.py \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv \
  GCA_947184155.2 --when doi:10.1126/sciadv.aba2498-c_goreaui
```

## Classify transcripts to KOs

```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --disable-cutoff --genome-accession _ \
  --fasta-file experiments/10.1126_sciadv.aba2498/c_goreaui.faa \
  kegg-downloads/ko.hmm _ experiments/10.1126_sciadv.aba2498/c_goreaui_ko.tsv

PYTHONPATH=. python3 scripts/classify/classify.py \
  --disable-cutoff --genome-accession _ \
  --fasta-file experiments/10.1126_sciadv.aba2498/aten.faa \
  kegg-downloads/ko.hmm _ experiments/10.1126_sciadv.aba2498/aten_ko.tsv

cat experiments/10.1126_sciadv.aba2498/c_goreaui_ko.tsv experiments/10.1126_sciadv.aba2498/aten_ko.tsv > experiments/10.1126_sciadv.aba2498/sequence_ko.tsv
```

Then filter the classified file to only keep classifications above the KO
threshold. For example,

```
python3 scripts/analysis/assign-ko.py experiments/10.1126_sciadv.aba2498/sequence_ko.tsv
mv experiments/10.1126_sciadv.aba2498/sequence_ko.tsv_filtered experiments/10.1126_sciadv.aba2498/sequence_ko.tsv
```

You can also directly classify/identify Pfam domains from the proteome

```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --genome-accession _ \
  --fasta-file data/exp_results/doi:10.1126_sciadv.aba2498/proteins.faa \
  pfam-downloads/Pfam-A.hmm _ data/exp_results/doi:10.1126_sciadv.aba2498/sequence_pfam.tsv \
   --cpu 4
```

Then, copy the .ko.tsv and .faa files to
data/exp_results/10.1126_sciadv.aba2498 directory. Usually, following these
conventions

  * `sequence_data.tsv`: data file
  * `sequence_ko.tsv`: sequence to KO mapping
  * `sequence_pfam.tsv`: sequence to Pfam mapping
  * `sequence_list.tsv`: from des2-merge script below
  * `des2_tall.tsv`: des2 data, from des2-merge script below


## Run DESeq2

DESeq2 analysis on RNAseq or proteomics results can be done using two scripts.
First, the `des2-simple.py` script does simple comparison of every two types of
cohorts at a given timepoint, or every two timepoints at the same cohort. For
example,

```
python3 scripts/analysis/des2-simple.py \
  --timepoint 1 --min-count 5 \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv data/exp_results/doi:10.1126_sciadv.aba2498

python3 scripts/analysis/des2-simple.py \
  --cohort SS8 --min-count 5 \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv data/exp_results/doi:10.1126_sciadv.aba2498
```

Then, the `des2-factor.py` script uses the `cohort_metadata.tsv` file to select
split by category, and then within each category, does comparison of every two
factor values. For example,

```
python3 scripts/analysis/des2-factor.py \
  --factor NutrientLevel --category Temperature \
  data/exp_results/doi:10.1093_pcp_pcac175/sequence_data.tsv \
  data/exp_results/doi:10.1093_pcp_pcac175/cohort_metadata.tsv \
  data/exp_results/doi:10.1093_pcp_pcac175

python3 scripts/analysis/des2-factor.py \
  --factor Temperature --category NutrientLevel \
  data/exp_results/doi:10.1093_pcp_pcac175/sequence_data.tsv \
  data/exp_results/doi:10.1093_pcp_pcac175/cohort_metadata.tsv \
  data/exp_results/doi:10.1093_pcp_pcac175
```

Use the following to merge multiple DES2 TSV files into a single tall TSV file,
first argument is output directory. IMPORTANT: for RNAseq data, make sure the
mapped file from Salmon, or other tool, has been converted to refer to protein
sequence IDs in the proteins.faa file.

```
PYTHONPATH=. python3 scripts/analysis/des2-merge.py \
  data/exp_results/doi:10.1126_sciadv.aba2498 \
  data/exp_results/doi:10.1126_sciadv.aba2498/proteins.faa \
  data/exp_results/doi:10.1126_sciadv.aba2498/deseq2_*.tsv
```
