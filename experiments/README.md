# Analysis Tasks


## Download data from NCBI SRA

To install tools

```
pip3 install pysradb
```

Use the following to get a list of SRR accessions to download

```
pysradb metadata PRJNA591730 > metadata.tsv
```

To extract fastq files, use NCBI SRA toolkit (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)

```
fasterq-dump SRR9331965
```


## Conda environments

Install Salmon using Conda into an environment, like "salmon_env". Also install samtools in that env.

```
conda create -n trinity_env -c bioconda -c conda-forge trinity
conda activate trinity_env
conda install -c bioconda bbmap salmon bowtie2 samtools transdecoder
```

Then run the following before you can start using salmon

```
conda init
source ~/.bash_profile
conda activate trinity_env
```


## Trinity de novo assembly of transcriptome

Reduce size of reads, starting with an existing reference genome.

```
bowtie2-build c_goreaui.fna c_goreaui.bowtie.index

bowtie2 -x GCA_014633955.1_index \
        -1 SRR9331959_1.fastq -2 SRR9331959_2.fastq \
        -p 8 \
        --no-unal \
        --un-conc SRR9331959_host_removed.fastq \
        -S SRR9331959_to_host.sam
bowtie2 --local -x c_goreaui.bowtie.index \
        -1 SRR9331959_host_removed.1.fastq -2 SRR9331959_host_removed.2.fastq \
        -p 8 \
        --no-unal \
        --al-conc SRR9331959_algae_captured.fastq \
        -S SRR9331959_algae_capture.sam

bbnorm.sh \
        in1=./SRR9331959_algae_captured.1.fastq \
        in2=./SRR9331959_algae_captured.2.fastq \
        out1=SRR9331959_norm.1.fq \
        out2=SRR9331959_norm.2.fq \
        target=50
```

De novo assembly

```
Trinity --seqType fq \         
        --left SRR9331959_norm.1.fq \
        --right SRR9331959_norm.2.fq \
        --max_memory 20G \
        --CPU 8 \
        --output trinity_SRR9331959_algae --no_normalize_reads

TrinityStats.pl trinity_SRR9331959_algae.Trinity.fasta
```

Cluster the transcripts, so if isoforms are grouped into different genes, they
can be clustered together. Note that this step is geared towards differential
functional analysis, not genetics. We are favoring aggregating, rather than
separating, alleles, for example.

```
scripts/cluster/mmseqs-cluster-trinity-transcripts \
  trinity_SRR9331959_algae.Trinity.fasta
```

Use TransDecoder to compute ORFs

```
TransDecoder.LongOrfs -t trinity_SRR9331959_algae.Trinity.fasta_rep_seq.fna.gz
TransDecoder.Predict -t trinity_SRR9331959_algae.Trinity.fasta_rep_seq.fna.gz -T 8
```


## Aligning against reference rather than assembled transcriptome

Building decoy-aware index, using genomic sequence as decoys, as reads may
contain UTRs or pre-mRNAs not in the curated transcriptome (often just START to
STOP).

```
grep "^>" c_goreaui.fna | cut -d " " -f 1 | sed -e 's/>//g' > c_goreaui.decoys.txt
cat c_goreaui.transcripts.fna c_goreaui.fna > c_goreaui.gentrome.fna   # transcript MUST COME FIRST
salmon index -t c_goreaui.gentrome.fna -d c_goreaui.decoys.txt -i c_goreaui.salmon_index

grep "^>" aten.fna | cut -d " " -f 1 | sed -e 's/>//g' > aten.decoys.txt
cat aten.transcripts.fna aten.fna > aten.gentrome.fna
salmon index -t aten.gentrome.fna -d aten.decoys.txt -i aten.salmon_index
```

## Quant using Salmon

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

## Process quants

Each directory must have its own script to process quants into `sequence_data.tsv`.

You can put multiple genomes in a single dataset file, i.e.
`sequence_data.tsv`. This is preferred if multiple species are present in a
sample, as DESeq2 can normalize across samples, based on more data per sample.
However, make sure, in the mapping workflow, that sequence IDs are unique
between species.


```
python3 experiments/doi:10.1126_sciadv.aba2498/process_salmon_quants.py \
  experiments/doi:10.1126_sciadv.aba2498 data/exp_results/doi:10.1126_sciadv.aba2498
```

Helpers to update the genome accession to something more formal

```
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

First, identify an internal control gene, EIF5B has been shown to work well for
algae and likely coral as well.

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

Use the `--include-timepoint 0` flag with the `--cohort` flag to limit to
specific timepoints, when comparing timepoint data within a cohort.

There is also a `des2-specific.py` script to compare any two cohort and/or
timepoints

```
python3 scripts/analysis/des2-specific.py --cohort1 34C --timepoint1 0 --cohort2 27C --timepoint2 192 --min-count 5 \
  data/exp_results/doi:10.1093_ismejo_wraf268/sequence_data.tsv data/exp_results/doi:10.1093_ismejo_wraf268
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


## Investigate transcript fidelity

Use bowtie2 to align reads against a reference genome, already indexed. to
obtain alignments with soft clipping. This allows investigation of specific
transcripts, how reads line up against the transcript, to see if there may be
multiple copies of portions of the transcripts may exist elsewhere in
transcriptome, or highly expressed domains either incorrectly assembled or
mapped at this transcript, or if the transcript is a full legitimate
transcript.

```
bowtie2-build c_goreaui.transcriptome.fna c_goreaui.transcriptome.bowtie.index

bowtie2 --local --very-sensitive-local -p 8 \
  -x c_goreaui.transcriptome.bowtie.index \
  -1 SRR9331959_1.fastq -2 SRR9331959_2.fastq \         
  -S SRR9331959_c_goreaui.sam

samtools view -S -b SRR9331959_c_goreaui.sam > SRR9331959_unsorted.bam
samtools sort SRR9331959_unsorted.bam -o SRR9331959_sorted.bam
samtools index SRR9331959_sorted.bam

samtools view SRR9331959_sorted.bam "lcl|CAMXCT020004035.1_cds_CAL1160776.1_35975" > SRR9331959_CAL1160776.1.sam
python3 scripts/analysis/pileup.py SRR9331959_CAL1160776.1.sam SRR9331959_CAL1160776.1.png
```

You can also search for a specific sequence using bowtie2, even though other
tools may be better for this

```
bowtie2 --local -f -p 8 \
        -x c_goreaui.transcriptome.bowtie.index \
        -U query.fasta \
        -S results.sam
```
