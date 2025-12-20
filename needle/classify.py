# Classification output format: protein id, genome accession, hmm-db, hmm
# accession, hmm coordinates, protein coordinates, dom evalue, dom score,
# threshold, pass threshold, protein query hit rank

import os
import csv
import itertools
from .match import ProteinsTSV
from .hmm import hmmscan_file


class ClassifyTSV(object):

    HEADERS = [
        "protein_accession",     # 0
        "genome_accession",      # 1
        "hmm_accession",         # 2
        "hmm_start",             # 3
        "hmm_end",               # 4
        "protein_start",         # 5
        "protein_end",           # 6
        "dom_evalue_cond",       # 7
        "dom_evalue",            # 8
        "dom_score",             # 9
        "score_threshold",       # 10
        "dom_rank_for_protein"   # 11
    ]

    @staticmethod
    def to_tsv_from_hmmscan_rows(tsv_path, hmm_rows, protein_defs, score_threshold_dict):

        # hmmscan (not hmmsearch): target is hmm profile, query is protein ID

        tacf = lambda row: row["target_accession"].strip() if row["target_accession"].strip() and row["target_accession"].strip() != "-" else row["target_name"].strip()
        qacf = lambda row: row["query_accession"].strip() if row["query_accession"].strip() and row["query_accession"].strip() != "-" else row["query_name"].strip()
        sorted_eval_for_protein = {}
        hmm_rows = sorted(hmm_rows, key=qacf)
        for protein_accession, group in itertools.groupby(hmm_rows, key=qacf):
            evalues = sorted([row["dom_evalue"] for row in list(group)])
            sorted_eval_for_protein[protein_accession] = evalues

        with open(tsv_path, "w") as f:
            writer = csv.DictWriter(f, fieldnames=ClassifyTSV.HEADERS, delimiter='\t')
            writer.writeheader()

            for row in hmm_rows:
                protein_accession = qacf(row)
                genome_accession = protein_defs[protein_accession]["genome_accession"]
                hmm_accession = tacf(row)
                score_threshold = "" if score_threshold_dict is None or hmm_accession not in score_threshold_dict else score_threshold_dict[hmm_accession]

                data = {
                    ClassifyTSV.HEADERS[0]:  protein_accession,
                    ClassifyTSV.HEADERS[1]:  genome_accession,
                    ClassifyTSV.HEADERS[2]:  hmm_accession,
                    ClassifyTSV.HEADERS[3]:  row["hmm_from"],
                    ClassifyTSV.HEADERS[4]:  row["hmm_to"],
                    ClassifyTSV.HEADERS[5]:  row["ali_from"],
                    ClassifyTSV.HEADERS[6]:  row["ali_to"],
                    ClassifyTSV.HEADERS[7]:  row["dom_evalue_cond"],
                    ClassifyTSV.HEADERS[8]:  row["dom_evalue"],
                    ClassifyTSV.HEADERS[9]:  row["dom_score"],
                    ClassifyTSV.HEADERS[10]: score_threshold,
                    ClassifyTSV.HEADERS[11]: sorted_eval_for_protein[protein_accession].index(row["dom_evalue"])+1
                }
                writer.writerow(data)


def classify(hmm_file, proteins_faa, proteins_tsv, cutoff_ga, output_tsv_path, score_threshold_dict):

    protein_match_rows = ProteinsTSV.from_tsv_to_rows(proteins_tsv)
    protein_defs = {
      row["protein_hit_id"]: dict(genome_accession=row["genome_accession"], target_accession=row["target_accession"])
      for row in protein_match_rows
    }

    hmm_rows = hmmscan_file(hmm_file, proteins_faa, cutoff=cutoff_ga)
    ClassifyTSV.to_tsv_from_hmmscan_rows(output_tsv_path, hmm_rows, protein_defs, score_threshold_dict)
