Running Salmon to map reads
---------------------------

Install Salmon using Conda into an environment, like "salmon_env"

Then run the following before you can start using salmon

```
conda init
source ~/.bash_profile
conda activate salmon_env
```

0. Obtaining fastq files

Download SRA files from NCBI at: https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=549921

To extract fastq files, use NCBI SRA toolkit (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)

```
fasterq-dump SRR9331965
```

1. Building decoy-aware inde

```
grep "^>" c_goreaui.fna | cut -d " " -f 1 | sed -e 's/>//g' > c_goreaui.decoys.txt
cat c_goreaui.transcripts.fna c_goreaui.fna > c_goreaui.gentrome.fna   # transcript MUST COME FIRST
salmon index -t c_goreaui.gentrome.fna -d c_goreaui.decoys.txt -i c_goreaui.salmon_index

grep "^>" aten.fna | cut -d " " -f 1 | sed -e 's/>//g' > aten.decoys.txt
cat aten.transcripts.fna aten.fna > aten.gentrome.fna
salmon index -t aten.gentrome.fna -d aten.decoys.txt -i aten.salmon_index
```

2. Quant using salmon

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

3. Process quants

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

4. Classify transcripts to KOs

```
PYTHONPATH=. python3 scripts/classify/classify.py \
  --disable-cutoff --genome-accession _ \
  --fasta-file experiments/10.1126_sciadv.aba2498/c_goreaui.faa \
  kegg-downloads/ko.hmm _ experiments/10.1126_sciadv.aba2498/c_goreaui.ko.tsv

PYTHONPATH=. python3 scripts/classify/classify.py \
  --disable-cutoff --genome-accession _ \
  --fasta-file experiments/10.1126_sciadv.aba2498/aten.faa \
  kegg-downloads/ko.hmm _ experiments/10.1126_sciadv.aba2498/aten.ko.tsv
```

Then, copy the .ko.tsv and .faa files to
data/exp_results/10.1126_sciadv.aba2498 directory.

5. Run DESeq2, see repo README.md
