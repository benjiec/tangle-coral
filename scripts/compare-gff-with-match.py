#!/usr/bin/env python3
import os
import re
import sys
import csv
import tempfile
import argparse
import itertools
import subprocess
from typing import List, Dict, Set

from needle.detect import DOM_EVALUE_LIMIT
from needle.match import extract_subsequence_strand_sensitive, read_fasta_as_dict
from needle.hmm import hmmscan_file
from needle.gff import parse_gff_to_hits
from defaults import DefaultPath


if __name__ == "__main__":
    # gff_proteins: query = protein accession, target = contig accession
    # hmm scan rows: query = protein, target = hmm name
    # needle rows: query = hmm name, target = contig accession

    ap = argparse.ArgumentParser(description="Join hmmscan domtblout with GFF-derived proteins.")
    ap.add_argument("hmm_file", help="Query HMM file")
    ap.add_argument("genome_accession")
    ap.add_argument("needle_match_file")
    ap.add_argument("--output-file", default=None)
    args = ap.parse_args()

    genome_accession = args.genome_accession
    faa_file = DefaultPath.ncbi_genome_faa(genome_accession)
    gff_file = DefaultPath.ncbi_genome_gff(genome_accession)
    fna_file = DefaultPath.ncbi_genome_fna(genome_accession)
    genomic_fasta = read_fasta_as_dict(fna_file)

    gff_proteins = parse_gff_to_hits(gff_file, protein_id_attr="protein_id", genomic_sequences=genomic_fasta)
    print("parsed GFF protein hits", len(gff_proteins))

    all_needle_rows = []
    with open(args.needle_match_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
           row["query_start"] = int(row["query_start"])
           row["query_end"] = int(row["query_end"])
           row["target_start"] = int(row["target_start"])
           row["target_end"] = int(row["target_end"])
           all_needle_rows.append(row)
    print("total Needle rows", len(all_needle_rows))

    hmmscan_rows = hmmscan_file(args.hmm_file, faa_file, cutoff=False)
    hmmscan_rows = [row for row in hmmscan_rows if row["dom_evalue"] < DOM_EVALUE_LIMIT]
    print("hmmscan", len(hmmscan_rows), "matches")

    keyf = lambda hmm_row: (hmm_row["target_name"], hmm_row["query_name"])
    hmmscan_rows = sorted(hmmscan_rows, key=keyf)

    out_rows = []
    for (hmm_name, hmm_query_protein_name), group in itertools.groupby(hmmscan_rows, key=keyf):

        gff_protein = [p for p in gff_proteins if p.query_accession == hmm_query_protein_name]
        assert len(gff_protein) == 1
        gff_protein = gff_protein[0]

        protein_hmm_rows = list(group)
        hmm_len = protein_hmm_rows[0]["target_length"]
        hmm_eval = protein_hmm_rows[0]["seq_evalue"]
        protein_len = protein_hmm_rows[0]["query_length"]

        hmm_aa_matched = [0 for i in range(0, hmm_len)]
        protein_aa_matched = [0 for i in range(0, protein_len)]
        for row in protein_hmm_rows:
            for hmm_i in range(row["hmm_from"]-1, row["hmm_to"]):
                hmm_aa_matched[hmm_i] = 1
            if row["ali_from"] < row["ali_to"]:
                for protein_i in range(row["ali_from"]-1, row["ali_to"]):
                    protein_aa_matched[protein_i] = 1
            else:
                for protein_i in range(row["ali_to"]-1, row["ali_from"]):
                    protein_aa_matched[protein_i] = 0

        hmm_aa_matched = sum(hmm_aa_matched)
        protein_aa_matched = sum(protein_aa_matched)
        protein_perc = round(protein_aa_matched * 100 / protein_len, 1)
        hmm_perc = round(hmm_aa_matched * 100 / hmm_len, 1)

        out = []
        out.append(genome_accession)
        out.append(hmm_query_protein_name)
        out.append(gff_protein.target_accession)
        out.append(gff_protein.target_start)
        out.append(gff_protein.target_end)
        out.append(protein_len)
        out.append(hmm_name)
        out.append(hmm_len)
        out.append(hmm_eval)
        out.append(hmm_aa_matched)
        out.append(hmm_perc)
        out.append(protein_aa_matched)
        out.append(protein_perc)

        needle_rows = [row for row in all_needle_rows \
                       if row["protein_hit_id"].startswith(hmm_name) and row["target_accession"] == gff_protein.target_accession]

        print(hmm_query_protein_name, "hmm", hmm_name, hmm_eval, len(protein_hmm_rows), "needle", len(needle_rows))
        for row in protein_hmm_rows:
            print("    hmm", row["target_name"], row["hmm_from"], row["hmm_to"], "matches", gff_protein.query_accession, row["ali_from"], row["ali_to"],
                  "of", row["query_length"], row["dom_evalue"])
        for m in gff_protein.matches:
            print("    gff", m.target_start, m.target_end, "=> protein", m.query_start, m.query_end)

        for row in needle_rows:
            print(" needle", row["protein_hit_id"], row["target_start"], row["target_end"], "=> hmm model", row["query_start"], row["query_end"])

        if len(needle_rows) == 0:
            status = "NOT FOUND"
            out.append(status)

        else:
            needle_target_left = min([min(row["target_start"], row["target_end"]) for row in needle_rows])
            needle_target_right = max([max(row["target_start"], row["target_end"]) for row in needle_rows])
            gff_target_left = min(gff_protein.target_start, gff_protein.target_end)
            gff_target_right = max(gff_protein.target_start, gff_protein.target_end)

            if needle_target_left > gff_target_right or needle_target_right < gff_target_left:
                status = "NOT FOUND (wrong locus)"
                out.append(status)

            else:
                needle_aa_matched_hmm = [0 for i in range(0, hmm_len)]
                for row in needle_rows:
                    for i in range(row["query_start"]-1, row["query_end"]):
                        needle_aa_matched_hmm[i] = 1
                needle_aa_matched = sum(needle_aa_matched_hmm)
                needle_aa_perc = round(needle_aa_matched * 100 / hmm_len, 1)

                status = "FOUND"
                substatus = []
                if float(abs(needle_aa_matched - hmm_aa_matched) / hmm_aa_matched) > 0.1:
                    substatus.append("len diff")
                if len(substatus):
                    status += " ("+", ".join(substatus)+")"

                out.append(status)
                out.append(needle_aa_matched)
                out.append(needle_aa_perc)
                out.append(hmm_perc)

        out_rows.append(out)

    if args.output_file:
        output_f = open(args.output_file, "w")
        output_f.write("genome\tprotein accession\ttarget accession\ttarget start\ttarget end\tprotein len\t"+\
                       "hmm name\thmm len\thmm evalue\thmm matched aa\thmm perc\tprotein matched aa\tprotein perc\tstatus\tneedle aa to hmm\tneedle matched perc\thmm perc\n")
        for out in out_rows:
            output_f.write("\t".join([str(x) for x in out])+"\n")
        output_f.close()
