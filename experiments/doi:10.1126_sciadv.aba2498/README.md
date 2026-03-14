## Read Mapping

Raw seq reads were downloaded from NCBI using BioProject code referenced in
paper. 

Decoy and transcript fasta files for C. goreaui and A. tenius were constructed.
Salmon was used to map reads onto these transcripts.

DESeq2 performs count normalization by default.

```
python3 scripts/analysis/des2-simple.py --timepoint 1 --min-count 5 \
  --genome-accession GCA_014633955.1 \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv data/exp_results/doi:10.1126_sciadv.aba2498

python3 scripts/analysis/des2-simple.py --timepoint 1 --min-count 5 \
  --genome-accession GCA_947184155.2 \
  data/exp_results/doi:10.1126_sciadv.aba2498/sequence_data.tsv data/exp_results/doi:10.1126_sciadv.aba2498

PYTHONPATH=. python3 scripts/analysis/des2-merge.py \
  data/exp_results/doi:10.1126_sciadv.aba2498 \
  data/exp_results/doi:10.1126_sciadv.aba2498/proteins.faa \
  data/exp_results/doi:10.1126_sciadv.aba2498/deseq2_*.tsv
```


## Looking at where reads map to

```
salmon quant -i c_goreaui.salmon_index \
             -l A \
             -1 SRR9331959_1.fastq \
             -2 SRR9331959_2.fastq \
             --validateMappings \
             -o c_goreaui-quants/SRR9331959 -p 2 --writeMappings=SRR9331959_c_goreaui.bam

samtools view -S -b SRR9331959_c_goreaui.bam > SRR9331959_unsorted.bam
samtools sort SRR9331959_unsorted.bam -o SRR9331959_sorted.bam
samtools index SRR9331959_sorted.bam

samtools view SRR9331959_sorted.bam "lcl|CAMXCT020004035.1_cds_CAL1160776.1_35975" > SRR9331959_CAL1160776.1.reads.txt

bowtie2 --local --very-sensitive-local -p 8 \
  -x c_goreaui.bowtie.index \
  -1 SRR9331959_1.fastq -2 SRR9331959_2.fastq \         
  -S SRR9331959_bowtie_local.sam
```
