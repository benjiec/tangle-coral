# Classification output format: protein id, genome accession, hmm-db, hmm
# accession, hmm coordinates, protein coordinates, dom evalue, dom score,
# threshold, pass threshold, protein query hit rank

import os
import csv
import itertools
from pathlib import Path
from .hmm import hmmscan_file


class ClassifyTSV(object):

    HEADERS = [
        "protein_accession",     # 0
        "genome_accession",      # 1
        "hmm_db",                # 2
        "hmm_accession",         # 3
        "hmm_start",             # 4
        "hmm_end",               # 5
        "protein_start",         # 6
        "protein_end",           # 7
        "dom_evalue_cond",       # 8
        "dom_evalue",            # 9
        "dom_score",             # 10
        "score_threshold",       # 11
        "dom_rank_for_protein"   # 12
    ]

    @staticmethod
    def to_tsv_from_hmmscan_rows(hmm_db_name, tsv_path, hmm_rows, protein_genome_accession_dict, score_threshold_dict):

        # hmmscan (not hmmsearch): target is hmm profile, query is protein ID

        tacf = lambda row: row["target_accession"].strip() if row["target_accession"].strip() and row["target_accession"].strip() != "-" else row["target_name"].strip()
        qacf = lambda row: row["query_accession"].strip() if row["query_accession"].strip() and row["query_accession"].strip() != "-" else row["query_name"].strip()
        sorted_eval_for_protein = {}
        hmm_rows = sorted(hmm_rows, key=qacf)
        for protein_accession, group in itertools.groupby(hmm_rows, key=qacf):
            evalues = sorted([row["dom_evalue"] for row in list(group)])
            sorted_eval_for_protein[protein_accession] = evalues

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
                    ClassifyTSV.HEADERS[0]:  protein_accession,
                    ClassifyTSV.HEADERS[1]:  genome_accession,
                    ClassifyTSV.HEADERS[2]:  hmm_db_name,
                    ClassifyTSV.HEADERS[3]:  hmm_accession,
                    ClassifyTSV.HEADERS[4]:  row["hmm_from"],
                    ClassifyTSV.HEADERS[5]:  row["hmm_to"],
                    ClassifyTSV.HEADERS[6]:  row["ali_from"],
                    ClassifyTSV.HEADERS[7]:  row["ali_to"],
                    ClassifyTSV.HEADERS[8]:  row["dom_evalue_cond"],
                    ClassifyTSV.HEADERS[9]:  row["dom_evalue"],
                    ClassifyTSV.HEADERS[10]:  row["dom_score"],
                    ClassifyTSV.HEADERS[11]: score_threshold,
                    ClassifyTSV.HEADERS[12]: sorted_eval_for_protein[protein_accession].index(row["dom_evalue"])+1
                }
                writer.writerow(data)


def classify(hmm_file, proteins_faa, cutoff_ga, output_tsv_path, protein_genome_accession_dict, score_threshold_dict, hmm_db_name = None):

    hmm_rows = hmmscan_file(hmm_file, proteins_faa, cutoff=cutoff_ga)
    if hmm_db_name is None:
        hmm_db_name = Path(hmm_file).stem
    ClassifyTSV.to_tsv_from_hmmscan_rows(hmm_db_name, output_tsv_path, hmm_rows, protein_genome_accession_dict, score_threshold_dict)
