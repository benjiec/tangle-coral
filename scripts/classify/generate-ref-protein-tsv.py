import os
import csv
import argparse
import contextlib
from needle.gff import parse_gff_to_hits
from scripts.defaults import DefaultPath
from needle.match import ProteinsTSV
from needle.classify import ClassifyTSV

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generates protein TSV for classified NCBI RefSeq proteins.")
    ap.add_argument("module_id")
    ap.add_argument("output_protein_tsv", help="Output Protein TSV file")
    ap.add_argument("output_name_tsv", help="Output Protein TSV file")
    ap.add_argument("--overwrite", help="Overwrite old data", action="store_true", default=False)
    args = ap.parse_args()

    classify_tsv = f"data/{args.module_id}_results/classify.tsv"
    classify_rows = ClassifyTSV.from_tsv_to_rows(classify_tsv)

    # get list of genome accessions that we have GFF files for
    genome_accessions = set(
        [row[ClassifyTSV.HDR_GENOME_ACCESSION] for row in classify_rows \
         if os.path.exists(DefaultPath.ncbi_genome_gff(row[ClassifyTSV.HDR_GENOME_ACCESSION]))]
    )

    classified_protein_accessions = set(
        [row[ClassifyTSV.HDR_PROTEIN_ACCESSION] for row in classify_rows \
         if row[ClassifyTSV.HDR_GENOME_ACCESSION] in genome_accessions]
    )

    print(f"{len(genome_accessions)} genome accessions")
    name_tsv_fieldnames = ["protein_accession", "name"]

    if args.overwrite:
        print("overwriting previous data")
        with contextlib.suppress(FileNotFoundError):
            os.remove(args.output_protein_tsv)
            os.remove(args.output_name_tsv)

        with open(args.output_name_tsv, "w") as f:
            writer = csv.DictWriter(f, fieldnames=name_tsv_fieldnames, delimiter='\t')
            writer.writeheader()
    else:
        print("appending to previous data")

    for gac in genome_accessions:
        print(gac)
        gff_path = DefaultPath.ncbi_genome_gff(gac)
        proteins = parse_gff_to_hits(gff_path)
        proteins = [p for p in proteins if p.query_accession in classified_protein_accessions]
        ProteinsTSV.append_to_tsv(args.output_protein_tsv, proteins, gac)

        with open(args.output_name_tsv, "a") as f:
            writer = csv.DictWriter(f, fieldnames=name_tsv_fieldnames, delimiter='\t')
            for p in proteins:
                if p._product_name is not None:
                    writer.writerow(dict(
                        protein_accession=p.query_accession,
                        name=p._product_name
                    ))
