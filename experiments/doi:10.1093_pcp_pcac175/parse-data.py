import csv

data_fn = "alldata.csv"

with open(data_fn, "r") as f:
    reader = csv.DictReader(f, delimiter=",")
    rows = [row for row in reader]

data = []

for row in rows:
    acc = row["accession"]
    for key, value in row.items():
        value = None if value == "NA" else value
        if key.startswith("LC"):
            data.append(dict(sequence_id=acc, genome_accession="GCA_001939145.1", count=float(value) if value is not None else "", timepoint=1, sample=key, cohort="LowNutrient25C"))
        elif key.startswith("LT"):
            data.append(dict(sequence_id=acc, genome_accession="GCA_001939145.1", count=float(value) if value is not None else "", timepoint=1, sample=key, cohort="LowNutrient34C"))
        elif key.startswith("HNC"):
            data.append(dict(sequence_id=acc, genome_accession="GCA_001939145.1", count=float(value) if value is not None else "", timepoint=1, sample=key, cohort="HighNutrient25C"))
        elif key.startswith("HNT"):
            data.append(dict(sequence_id=acc, genome_accession="GCA_001939145.1", count=float(value) if value is not None else "", timepoint=1, sample=key, cohort="HighNutrient34C"))
        elif key.startswith("HLC"):
            data.append(dict(sequence_id=acc, genome_accession="GCA_001939145.1", count=float(value) if value is not None else "", timepoint=1, sample=key, cohort="Imbalanced25C"))
        elif key.startswith("HLT"):
            data.append(dict(sequence_id=acc, genome_accession="GCA_001939145.1", count=float(value) if value is not None else "", timepoint=1, sample=key, cohort="Imbalanced34C"))

metadata = []
metadata.append(dict(cohort="LowNutrient25C", term="Temperature", value=25, unit="C"))
metadata.append(dict(cohort="LowNutrient25C", term="N", value=0))
metadata.append(dict(cohort="LowNutrient25C", term="P", value=0))
metadata.append(dict(cohort="HighNutrient25C", term="Temperature", value=25, unit="C"))
metadata.append(dict(cohort="HighNutrient25C", term="N", value=1))
metadata.append(dict(cohort="HighNutrient25C", term="P", value=1))
metadata.append(dict(cohort="Imbalanced25C", term="Temperature", value=25, unit="C"))
metadata.append(dict(cohort="Imbalanced25C", term="N", value=1))
metadata.append(dict(cohort="Imbalanced25C", term="P", value=0))
metadata.append(dict(cohort="LowNutrient34C", term="Temperature", value=34, unit="C"))
metadata.append(dict(cohort="LowNutrient34C", term="N", value=0))
metadata.append(dict(cohort="LowNutrient34C", term="P", value=0))
metadata.append(dict(cohort="HighNutrient34C", term="Temperature", value=34, unit="C"))
metadata.append(dict(cohort="HighNutrient34C", term="N", value=1))
metadata.append(dict(cohort="HighNutrient34C", term="P", value=1))
metadata.append(dict(cohort="Imbalanced34C", term="Temperature", value=34, unit="C"))
metadata.append(dict(cohort="Imbalanced34C", term="N", value=1))
metadata.append(dict(cohort="Imbalanced34C", term="P", value=0))

with open("sequence_data.tsv", "w") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(data[0].keys()))
    writer.writeheader()
    for d in data:
        writer.writerow(d)

with open("cohort_metadata.tsv", "w") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(metadata[0].keys()))
    writer.writeheader()
    for d in metadata:
        writer.writerow(d)
