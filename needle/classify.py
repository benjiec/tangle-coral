# Classification output format: protein id, genome accession, hmm-db, hmm
# accession, hmm coordinates, protein coordinates, dom evalue, dom score,
# threshold, pass threshold, protein query hit rank

import os
import csv
import itertools
from pathlib import Path
from .seq import read_fasta_as_dict, write_fasta_from_dict, extract_subsequence
from .hmm import hmmscan_file


class ClassifyTSV(object):

    HDR_PROTEIN_ACCESSION = "protein_accession"
    HDR_GENOME_ACCESSION = "genome_accession"
    HDR_HMM_DB = "hmm_db"
    HDR_HMM_ACCESSION = "hmm_accession"
    HDR_HMM_START = "hmm_start"
    HDR_HMM_END = "hmm_end"
    HDR_PROTEIN_START = "protein_start"
    HDR_PROTEIN_END = "protein_end"
    HDR_DOM_EVALUE_COND = "dom_evalue_cond"
    HDR_DOM_EVALUE = "dom_evalue"
    HDR_DOM_SCORE = "dom_score"
    HDR_SCORE_THRESHOLD = "score_threshold"
    HDR_DOM_RANK = "dom_rank_for_protein"

    HEADERS = [
        HDR_PROTEIN_ACCESSION,   # 0
        HDR_GENOME_ACCESSION,    # 1
        HDR_HMM_DB,              # 2
        HDR_HMM_ACCESSION,       # 3
        HDR_HMM_START,           # 4
        HDR_HMM_END,             # 5
        HDR_PROTEIN_START,       # 6
        HDR_PROTEIN_END,         # 7
        HDR_DOM_EVALUE_COND,     # 8
        HDR_DOM_EVALUE,          # 9
        HDR_DOM_SCORE,           # 10
        HDR_SCORE_THRESHOLD,     # 11
        HDR_DOM_RANK,            # 12
    ]

    @staticmethod
    def from_tsv_to_rows(tsv_path):
        rows = []
        with open(tsv_path, "r") as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                row = {k: row[k] for k in ClassifyTSV.HEADERS}
                row[ClassifyTSV.HDR_HMM_START] = int(row[ClassifyTSV.HDR_HMM_START])
                row[ClassifyTSV.HDR_HMM_END] = int(row[ClassifyTSV.HDR_HMM_END])
                row[ClassifyTSV.HDR_PROTEIN_START] = int(row[ClassifyTSV.HDR_PROTEIN_START])
                row[ClassifyTSV.HDR_PROTEIN_END] = int(row[ClassifyTSV.HDR_PROTEIN_END])
                row[ClassifyTSV.HDR_DOM_EVALUE_COND] = float(row[ClassifyTSV.HDR_DOM_EVALUE_COND])
                row[ClassifyTSV.HDR_DOM_EVALUE] = float(row[ClassifyTSV.HDR_DOM_EVALUE])
                row[ClassifyTSV.HDR_DOM_SCORE] = float(row[ClassifyTSV.HDR_DOM_SCORE])
                row[ClassifyTSV.HDR_SCORE_THRESHOLD] = float(row[ClassifyTSV.HDR_SCORE_THRESHOLD]) if row[ClassifyTSV.HDR_SCORE_THRESHOLD] else None
                row[ClassifyTSV.HDR_DOM_RANK] = int(row[ClassifyTSV.HDR_DOM_RANK])
                rows.append(row)

        return rows

    @staticmethod
    def to_tsv_from_hmmscan_rows(hmm_db_name, tsv_path, hmm_rows, protein_genome_accession_dict, score_threshold_dict):

        # hmmscan (not hmmsearch): target is hmm profile, query is protein ID

        tacf = lambda row: row["target_accession"].strip() if row["target_accession"].strip() and row["target_accession"].strip() != "-" else row["target_name"].strip()
        qacf = lambda row: row["query_accession"].strip() if row["query_accession"].strip() and row["query_accession"].strip() != "-" else row["query_name"].strip()
        sorted_score_for_protein = {}
        hmm_rows = sorted(hmm_rows, key=qacf)
        for protein_accession, group in itertools.groupby(hmm_rows, key=qacf):
            scores = sorted([row["dom_score"] for row in list(group)], reverse=True)
            sorted_score_for_protein[protein_accession] = scores

        if not os.path.exists(tsv_path) or os.path.getsize(tsv_path) == 0:
            with open(tsv_path, "w") as f:
                writer = csv.DictWriter(f, fieldnames=ClassifyTSV.HEADERS, delimiter='\t')
                writer.writeheader()

        with open(tsv_path, "a") as f:
            writer = csv.DictWriter(f, fieldnames=ClassifyTSV.HEADERS, delimiter='\t')

            for row in hmm_rows:
                protein_accession = qacf(row)
                genome_accession = protein_genome_accession_dict[protein_accession]
                hmm_accession = tacf(row)
                score_threshold = "" if score_threshold_dict is None or hmm_accession not in score_threshold_dict else score_threshold_dict[hmm_accession]

                data = {
                    ClassifyTSV.HDR_PROTEIN_ACCESSION:  protein_accession,
                    ClassifyTSV.HDR_GENOME_ACCESSION:  genome_accession,
                    ClassifyTSV.HDR_HMM_DB:  hmm_db_name,
                    ClassifyTSV.HDR_HMM_ACCESSION:  hmm_accession,
                    ClassifyTSV.HDR_HMM_START:  row["hmm_from"],
                    ClassifyTSV.HDR_HMM_END:  row["hmm_to"],
                    ClassifyTSV.HDR_PROTEIN_START:  row["ali_from"],
                    ClassifyTSV.HDR_PROTEIN_END:  row["ali_to"],
                    ClassifyTSV.HDR_DOM_EVALUE_COND:  row["dom_evalue_cond"],
                    ClassifyTSV.HDR_DOM_EVALUE:  row["dom_evalue"],
                    ClassifyTSV.HDR_DOM_SCORE:  row["dom_score"],
                    ClassifyTSV.HDR_SCORE_THRESHOLD: score_threshold,
                    ClassifyTSV.HDR_DOM_RANK: sorted_score_for_protein[protein_accession].index(row["dom_score"])+1
                }
                writer.writerow(data)


def classify(hmm_file, proteins_faa, cutoff_ga, output_tsv_path, protein_genome_accession_dict, score_threshold_dict, hmm_db_name = None):

    hmm_rows = hmmscan_file(hmm_file, proteins_faa, cutoff=cutoff_ga)
    if hmm_db_name is None:
        hmm_db_name = Path(hmm_file).stem
    ClassifyTSV.to_tsv_from_hmmscan_rows(hmm_db_name, output_tsv_path, hmm_rows, protein_genome_accession_dict, score_threshold_dict)


def group_by_assignment(classify_rows, ortholog_hmm_db_name, score_to_threshold_ratio):

    # assignment criteria
    # hmm_db == <ortholog_hmm_db_name> AND
    # dom_rank_for_protein == 1 AND
    # dom_score / score_threshold >= <score_to_threshold_ratio>

    assignment_rows = [row for row in classify_rows if row[ClassifyTSV.HDR_HMM_DB] == ortholog_hmm_db_name and \
                                                       row[ClassifyTSV.HDR_DOM_RANK] == 1 and \
                                                       row[ClassifyTSV.HDR_DOM_SCORE] / row[ClassifyTSV.HDR_SCORE_THRESHOLD] >= score_to_threshold_ratio]
                 
    assignments = {row[ClassifyTSV.HDR_PROTEIN_ACCESSION]: row[ClassifyTSV.HDR_HMM_ACCESSION] for row in assignment_rows}
    assigned_rows = [row for row in classify_rows if (row[ClassifyTSV.HDR_HMM_DB] != ortholog_hmm_db_name and \
                                                      row[ClassifyTSV.HDR_PROTEIN_ACCESSION] in assignments) or \
                                                     (row[ClassifyTSV.HDR_HMM_DB] == ortholog_hmm_db_name and \
                                                      row[ClassifyTSV.HDR_PROTEIN_ACCESSION] in assignments and \
                                                      row[ClassifyTSV.HDR_HMM_ACCESSION] == assignments[row[ClassifyTSV.HDR_PROTEIN_ACCESSION]])]

    keyf = lambda row: assignments[row[ClassifyTSV.HDR_PROTEIN_ACCESSION]]
    assigned_rows = sorted(assigned_rows, key=keyf)
    return itertools.groupby(assigned_rows, keyf)


def assign_ko(classify_rows, ortholog_hmm_db_name, proteins_faa, output_dir, domain_hmm_db_name = None, score_to_threshold_ratio = 0.75):

    assignments = group_by_assignment(classify_rows, ortholog_hmm_db_name, score_to_threshold_ratio)

    if type(proteins_faa) == type(""):
        proteins_faa = [proteins_faa]

    proteins_seq_dict = {}
    for faa_file in proteins_faa:
        proteins_seq_dict |= read_fasta_as_dict(faa_file)

    for ko_id, ko_rows in assignments:
        ko_rows = list(ko_rows)
        ko_fasta_prefix = os.path.join(output_dir, ko_id)
        protein_ids = set([row[ClassifyTSV.HDR_PROTEIN_ACCESSION] for row in ko_rows])

        # protein fasta - just protein sequences assigned to KO
        ko_proteins = {k:v for k,v in proteins_seq_dict.items() if k in protein_ids}
        write_fasta_from_dict(ko_proteins, ko_fasta_prefix+".faa")

        if domain_hmm_db_name is not None:
            domain_rows = [row for row in ko_rows if row[ClassifyTSV.HDR_HMM_DB] == domain_hmm_db_name]
            keyf = lambda row: row[ClassifyTSV.HDR_HMM_ACCESSION]
            domain_rows = sorted(domain_rows, key=keyf)

            for domain_accession, group in itertools.groupby(domain_rows, keyf):
                domain_seqs = {}
                for domain_row in group:
                    protein_id = domain_row[ClassifyTSV.HDR_PROTEIN_ACCESSION]
                    protein_start = domain_row[ClassifyTSV.HDR_PROTEIN_START]
                    protein_end = domain_row[ClassifyTSV.HDR_PROTEIN_END]
                    seq = extract_subsequence(proteins_seq_dict[protein_id], protein_start, protein_end)
                    domain_seq_accession = f"{protein_id}_{protein_start}_{protein_end}_{domain_accession}"
                    domain_seqs[domain_seq_accession] = seq
                write_fasta_from_dict(domain_seqs, ko_fasta_prefix+"_"+domain_accession+".faa")
