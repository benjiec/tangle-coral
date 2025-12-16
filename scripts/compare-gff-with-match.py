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

from needle.match import extract_subsequence_strand_sensitive, read_fasta_as_dict
from needle.match import compute_three_frame_translations
from needle.hmm import hmmscan_file, HMMCollection
from needle.gff import parse_gff_to_hits
from needle.io import export_protein_hits
from defaults import DefaultPath


def protein_enzyme_equality_score(enzyme_name, protein_name):
    score = 0
    enzyme_words = re.split(r"\W+", enzyme_name)
    for i, word in enumerate(enzyme_words):
        if word.endswith("ase") or word == "metabolism":
            if i > 0:
                if len(enzyme_words[i-1]) > 3 and enzyme_words[i-1] in protein_name:
                    score += 1
                    if enzyme_words[i-1]+" "+word in protein_name:
                        score *= 2
                elif len(enzyme_words[i-2]) > 3 and enzyme_words[i-2] in protein_name:
                    score += 1
                    if enzyme_words[i-2]+" "+word in protein_name:
                        score *= 2
    return score


def load_accessions(fasta_path):
    d = read_fasta_as_dict(fasta_path)
    return list(d.keys())


def read_ko_names(path: str) -> Dict[str, str]:
    name_by_accession = {}
    with open(path, "r") as f:
        for raw_line in f:
            if not raw_line:
                continue
            line = raw_line.rstrip("\n")
            if not line:
                continue
            header_content = line.strip()
            accession = header_content.split("\t")[0]
            name_by_accession[accession] = " ".join(header_content.split("\t")[1:])
    return name_by_accession


def read_fasta_sequence_names(path: str) -> Dict[str, str]:
    name_by_accession = {}
    with open(path, "r") as f:
        for raw_line in f:
            if not raw_line:
                continue
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                header_content = line[1:].strip()
                accession = header_content.split(None, 1)[0]
                name_by_accession[accession] = header_content
    return name_by_accession


def run_cmd(cmd: List[str]) -> None:
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {' '.join(cmd)}", file=sys.stderr)
        raise


def print_comparison(protein_name, status, hmm_rows, gff_hit, needle_rows):
    print("")
    print(protein_name)
    print(status)

    protein_aa_len = 0
    protein_len = 0
    for row in hmm_rows:
        protein_aa_len += (abs(row["ali_to"] - row["ali_from"]) + 1)
        protein_len = row["query_length"]

    print("    found by hmmscan, on", gff_hit.target_accession, gff_hit.target_start, gff_hit.target_end, "of", gff_hit.query_accession)
    for row in hmm_rows:
        print("    hmm model", row["target_name"], row["target_accession"], row["hmm_from"], row["hmm_to"], "matches", gff_hit.query_accession, row["ali_from"], row["ali_to"],
              "of", row["query_length"], row["dom_evalue"], row["dom_score"])

    if not needle_rows:
        print("    DID NOT FIND NEEDLE RESULT")
        for m in gff_hit.matches:
            print("    gff match", m.target_start, m.target_end, m.query_start, m.query_end)
    else:
        for row in needle_rows:
            print("    needle dna match", row["protein_hit_id"], row["target_start"], row["target_end"], "covering hmm model", row["query_start"], row["query_end"]) 
        for m in gff_hit.matches:
            print("    gff match", m.target_start, m.target_end, m.query_start, m.query_end)


def compare(hmm_file, query_accession, genome_accession, gff_proteins, protein_seq_names, genomic_fasta, needle_rows, output_f, not_found_f):

    # Run hmmscan producing domtblout
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = ["hmmpress", "-f", hmm_file]
        run_cmd(cmd)
        hmmscan_rows = hmmscan_file(hmm_file, faa_file)

    # gff_proteins: query = protein accession, target = contig accession
    # hmm rows: query = protein, target = hmm name
    # needle rows: query = hmm name, target = contig accession

    print(hmm_file)
    print("received", len(hmmscan_rows), "matches from hmmscan")
    queries_from_hmmscan = list(set([m["query_name"] for m in hmmscan_rows]))

    # Filter to proteins that showed up in hmmscan output (by protein_id == query name)
    filtered_hits = [pm for pm in gff_proteins if pm.query_accession in queries_from_hmmscan]
    print("gff proteins filtered by hmmscan results to", len(filtered_hits))
    target_accessions = list(set([pm.target_accession for pm in filtered_hits]))

    # Load genomic fasta and set match target sequence
    for pm in filtered_hits:
        for m in pm.matches:
            m.target_sequence = extract_subsequence_strand_sensitive(genomic_fasta[m.target_accession], m.target_start, m.target_end)

    needle_rows = [row for row in needle_rows if row["protein_hit_id"].startswith(query_accession) and row["target_accession"] in target_accessions]
    print("needle results matching", query_accession, "on same target accessions", len(needle_rows))

    def sorter(gff_hit):
        hmm_rows = [row for row in hmmscan_rows if row["query_name"] == gff_hit.query_accession]
        return max([row["hmm_to"] for row in hmm_rows]) - \
               min([row["hmm_from"] for row in hmm_rows])
    filtered_hits = sorted(filtered_hits, key=sorter, reverse=True)

    not_found_hits = []

    for gff_hit in filtered_hits:
        hmm_rows = [row for row in hmmscan_rows if row["query_name"] == gff_hit.query_accession]
        hmm_len = hmm_rows[0]["target_length"]
        hmm_eval = hmm_rows[0]["seq_evalue"]
        protein_len = hmm_rows[0]["query_length"]

        hmm_aa_matched = [0 for i in range(0, hmm_len)]
        protein_aa_matched = [0 for i in range(0, protein_len)]
        for row in hmm_rows:
            for hmm_i in range(row["hmm_from"]-1, row["hmm_to"]):
                hmm_aa_matched[hmm_i] = 1
            if row["ali_from"] < row["ali_to"]:
                for protein_i in range(row["ali_from"]-1, row["ali_to"]):
                    protein_aa_matched[protein_i] = 1
            else:
                for protein_i in range(row["ali_to"]-1, row["ali_from"]):
                    protein_aa_matched[protein_i] = 0

        # compute the exon length where hmm matched
        matched_exons = []
        for m in gff_hit.matches:
            fragment_matched_by_hmm = False
            for i in range(m.query_start-1, m.query_end):
                if i < len(protein_aa_matched) and protein_aa_matched[i] > 0:
                    fragment_matched_by_hmm = True
            if fragment_matched_by_hmm:
                matched_exons.append(abs(m.target_end-m.target_start)+1)
        matched_exon_len = round(sum(matched_exons) / len(matched_exons), 1)

        hmm_aa_matched = sum(hmm_aa_matched)
        protein_aa_matched = sum(protein_aa_matched)
        protein_perc = round(protein_aa_matched * 100 / protein_len, 1)
        hmm_perc = round(hmm_aa_matched * 100 / hmm_len, 1)

        gff_dna_len = 0
        for m in gff_hit.matches:
            gff_dna_len += (abs(m.target_start - m.target_end) + 1)
        gff_target_left = min([min(m.target_start, m.target_end) for m in gff_hit.matches])
        gff_target_right = max([max(m.target_start, m.target_end) for m in gff_hit.matches])

        protein_name = protein_seq_names[gff_hit.query_accession]
        protein_acc = gff_hit.query_accession

        needle_rows_on_target = [row for row in needle_rows if row["target_accession"] == gff_hit.target_accession]

        if len(needle_rows_on_target) == 0:
            status = "NOT FOUND"
            print_comparison(protein_name, status, hmm_rows, gff_hit, None)
            if output_f:
                output_f.write(f"{query_accession}\t{genome_accession}\t{protein_name}\t{protein_acc}\t{status}\t"+\
                               f"{hmm_len}\t{hmm_perc}\t{protein_len}\t{protein_perc}\t{matched_exon_len}\t{hmm_eval}\t\t0\t{hmm_aa_matched}\n")
            not_found_hits.append(gff_hit)

        else:
            keyf = lambda row: row["protein_hit_id"]
            needle_rows_on_target = sorted(needle_rows_on_target, key=keyf)

            for needle_protein_hit_id, needle_rows_for_hit in itertools.groupby(needle_rows_on_target, keyf):
                needle_rows_for_hit = list(needle_rows_for_hit)

                needle_target_left = min([min(row["target_start"], row["target_end"]) for row in needle_rows_for_hit])
                needle_target_right = max([max(row["target_start"], row["target_end"]) for row in needle_rows_for_hit])

                if needle_target_left > gff_target_right or \
                   needle_target_right < gff_target_left:
                    status = "NOT SAME LOCUS"
                    # print_comparison(protein_name, status, hmm_rows, gff_hit, needle_rows_for_hit)

		    # not actually reporting this, because we are also not
		    # reporting other needle matches that were not discovered
		    # by hmmscan

                else:
                    needle_aa_matched_hmm = [0 for i in range(0, hmm_len)]
                    needle_dna_len = 0
                    for row in needle_rows_for_hit:
                        for i in range(row["query_start"]-1, row["query_end"]):
                            needle_aa_matched_hmm[i] = 1
                        needle_dna_len += (abs(row["target_start"] - row["target_end"]) + 1)
                    needle_aa_matched = sum(needle_aa_matched_hmm)

                    status = "FOUND"
                    substatus = []
                    if float(abs(needle_aa_matched - hmm_aa_matched) / hmm_aa_matched) > 0.1:
                        substatus.append("LEN DIFF HMMSCAN")
                    if len(substatus):
                        status += " ("+", ".join(substatus)+")"

                    print_comparison(protein_name, status, hmm_rows, gff_hit, needle_rows_for_hit)
                    if output_f:
                        output_f.write(f"{query_accession}\t{genome_accession}\t{protein_name}\t{protein_acc}\t{status}\t"+\
                                       f"{hmm_len}\t{hmm_perc}\t{protein_len}\t{protein_perc}\t{matched_exon_len}\t{hmm_eval}\t"+\
                                       f"{needle_protein_hit_id}\t{needle_aa_matched}\t{hmm_aa_matched}\n")

    if output_f:
        output_f.flush()

    if not_found_f:
        for gff_hit in not_found_hits:
            for m in gff_hit.matches:
                x3f = compute_three_frame_translations(m.target_sequence, 1, len(m.target_sequence))
                for frame_start, _2, aa in x3f:
                  not_found_f.write(f">{gff_hit.protein_hit_id}_{m.target_start}_{frame_start}\n{aa}\n")
        not_found_f.flush()


if __name__ == "__main__":

    ap = argparse.ArgumentParser(description="Join hmmscan domtblout with GFF-derived proteins.")
    ap.add_argument("query_fasta", help="Query fasta used for searching")
    ap.add_argument("genome_accession")
    ap.add_argument("needle_match_file")
    ap.add_argument("hmm_collection", help="HMM collection, either ko or pfam", choices=['ko', 'pfam'])
    ap.add_argument("--output-file", default=None)
    ap.add_argument("--not-found-file", default=None)
    args = ap.parse_args()

    genome_accession = args.genome_accession
    ncbi_download_dir = DefaultPath.ncbi_genome_dir(genome_accession)

    faa_file = DefaultPath.ncbi_genome_faa(genome_accession)
    gff_file = DefaultPath.ncbi_genome_gff(genome_accession)
    fna_file = DefaultPath.ncbi_genome_fna(genome_accession)
    assert fna_file
    protein_seq_names = read_fasta_sequence_names(faa_file)
    genomic_fasta = read_fasta_as_dict(fna_file)

    # Parse GFF -> ProteinHit list
    gff_proteins = parse_gff_to_hits(gff_file, protein_id_attr="protein_id", genomic_sequences=genomic_fasta)
    print("parsed protein hits from GFF, total", len(gff_proteins))

    needle_rows = []
    with open(args.needle_match_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
           row["query_start"] = int(row["query_start"])
           row["query_end"] = int(row["query_end"])
           row["target_start"] = int(row["target_start"])
           row["target_end"] = int(row["target_end"])
           needle_rows.append(row)
    print("total needle rows", len(needle_rows))

    output_f = None
    if args.output_file:
        output_f = open(args.output_file, "w")
        output_f.write("query\tgenome\tprotein name\tprotein accession\tstatus\t"+\
                       "hmm len\thmm match perc\tprotein len\tprotein match perc\texon avg len\thmmscan evalue\t"+\
                       "needle protein id\tneedle aa matched\thmm aa matched\n")

    not_found_f = None
    if args.not_found_file:
        not_found_f = open(args.not_found_file, "w")

    accessions = load_accessions(args.query_fasta)
    if args.hmm_collection == "ko":
        hmm_collection = HMMCollection(DefaultPath.ko_hmm(), accessions)
    else:
        hmm_collection = HMMCollection(DefaultPath.pfam_hmm(), accessions)

    try:
        for acc in accessions:
            compare(hmm_collection.get(acc), acc, genome_accession, gff_proteins, protein_seq_names, genomic_fasta, needle_rows, output_f, not_found_f)
    finally:
        hmm_collection.clean()

    if output_f:
        output_f.close()
    if not_found_f:
        not_found_f.close()
