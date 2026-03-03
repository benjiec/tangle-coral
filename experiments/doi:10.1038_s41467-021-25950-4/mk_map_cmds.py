manifest = """
SRR6255820: orbicella, control, A
SRR6255822: orbicella, control, B
SRR6255825: orbicella, control, C
SRR6255844: orbicella, treatment, A
SRR6255862: orbicella, treatment, B
SRR6255864: orbicella, treatment, C
SRR6255880: pseudodiploria, control, A
SRR6255882: pseudodiploria, control, B
SRR6255888: pseudodiploria, control, C
SRR6255887: pseudodiploria, treatment, A
SRR6255886: pseudodiploria, treatment, B
SRR6255890: pseudodiploria, treatment, C
SRR6255889: siderastrea, control, A
SRR6256312: siderastrea, control, B
SRR6256323: siderastrea, control, C
SRR6256322: siderastrea, treatment, A
SRR6256324: siderastrea, treatment, B
SRR6256325: siderastrea, treatment, C
"""

for line in manifest.strip().split("\n"):
    sra_rec = line.split(": ")[0]
    holobiont = line.split(": ")[1].split(",")[0]
    print(f"salmon quant -i {holobiont}_host.salmon_index -l A --validateMappings -p 2 -o host_quants/{sra_rec} -1 reads/{sra_rec}_1.fastq -2 reads/{sra_rec}_2.fastq")
    print(f"salmon quant -i {holobiont}_symb.salmon_index -l A --validateMappings -p 2 -o symb_quants/{sra_rec} -1 reads/{sra_rec}_1.fastq -2 reads/{sra_rec}_2.fastq")
