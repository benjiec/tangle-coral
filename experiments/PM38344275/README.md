This paper provides three separate transcriptomes, assembled using Trinity. The
NGS reads were from other papers.

Transcriptomes https://pubmed.ncbi.nlm.nih.gov/38344275/
  * Cladocopium proliferum SCF055, reads from PRJNA723630, https://pubmed.ncbi.nlm.nih.gov/35383179/
  * Durusdinium trenchii CCMP2556, reads from PRJNA508937, https://pubmed.ncbi.nlm.nih.gov/31693775/
  * Prorocentrum cordatum CCMP1329, reads from PRJEB54915, https://pubmed.ncbi.nlm.nih.gov/37996937/

The downloaded transcripts, from PM38344275 data archive, as

```
cladocopium_proliferum.fna
durusdinium_trenchii.fna
prorocentrum_cordatum.fna
```

Use the following to uniquefy transcript names

```
coral-py coral/experiments/PM34593802/unique-acc.py \
  cladocopium_proliferum.fna
coral-py coral/experiments/PM34593802/unique-acc.py \
  durusdinium_trenchii.fna
coral-py coral/experiments/PM34593802/unique-acc.py \
  prorocentrum_cordatum.fna
```

Then move these .fna files into separate pile transcriptome directories, under
pile project PM38344275.

See Pile README. Run TransDecoder to generate protein fasta files. Save the
.gff3 files in case we need them later.

Concatentate the protein FASTAs together into one `proteins.faa`.
