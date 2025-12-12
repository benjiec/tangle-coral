import re
import os
import shutil
import subprocess
import tempfile


def run_command(cmd: str):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def parse_hmmsearch_domtbl(domtbl_path):
    expected_header = "# target name  accession  tlen  query name  accession  qlen  E-value  score  bias  #  of  c-Evalue  i-Evalue  score  bias  from  to  from  to"
    idx_target = 0
    idx_target_acc = 1
    idx_t_len = 2
    idx_query = 3
    idx_q_len = 5
    idx_eval = 6
    idx_score = 7
    idx_h_from = 15
    idx_h_to = 16
    idx_a_from = 17
    idx_a_to = 18

    # first, sanity check these indices
    expected_header_parts = re.split(r'\s\s+', expected_header)
    assert expected_header_parts[idx_eval] == "E-value"
    assert expected_header_parts[idx_score] == "score"
    assert expected_header_parts[idx_target] == "# target name"
    assert expected_header_parts[idx_target_acc] == "accession"
    assert expected_header_parts[idx_t_len] == "tlen"
    assert expected_header_parts[idx_query] == "query name"
    assert expected_header_parts[idx_q_len] == "qlen"
    assert expected_header_parts[idx_h_from] == "from"
    assert expected_header_parts[idx_h_to] == "to"
    assert expected_header_parts[idx_a_from] == "from"
    assert expected_header_parts[idx_a_to] == "to"

    has_headers = False
    expected_header = " ".join(expected_header.split())

    matches = []
    with open(domtbl_path, "r") as domf:
        for line in domf:
            if " ".join(line.split()).startswith(expected_header):
                has_headers = True
            if not line or line.startswith("#") or has_headers is False:
                continue
            parts = line.strip().split()
            match = dict(
                target_name = parts[idx_target],
                target_accession = parts[idx_target_acc],
                query_name = parts[idx_query],
                evalue = float(parts[idx_eval]),
                score = float(parts[idx_score]),
                query_length = int(parts[idx_q_len]),
                hmm_from = int(parts[idx_h_from]),
                hmm_to = int(parts[idx_h_to]),
                target_length = int(parts[idx_t_len]),
                ali_from = int(parts[idx_a_from]),
                ali_to = int(parts[idx_a_to])
            )
            matches.append(match)

    assert has_headers
    return matches


def hmmsearch(hmm_file_name, sequences):
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "cands.faa")
        domtbl_path = os.path.join(tmpdir, "out.domtbl")
        with open(fasta_path, "w") as f:
            for i, cand in enumerate(sequences):
                f.write(f">cand_{i}\n{cand}\n")
        cmd = ["hmmsearch", "--cut_ga", "--domtblout", domtbl_path, hmm_file_name, fasta_path]
        run_command(cmd)
        return parse_hmmsearch_domtbl(domtbl_path)


def hmmscan_file(hmm_file_name, fasta_path):
    with tempfile.TemporaryDirectory() as tmpdir:
        domtbl_path = os.path.join(tmpdir, "out.domtbl")
        cmd = ["hmmscan", "--cut_ga", "--domtblout", domtbl_path, hmm_file_name, fasta_path]
        run_command(cmd)
        return parse_hmmsearch_domtbl(domtbl_path)


def hmmfetch(hmm_file_name, acc, output_file_name):
    cmd = ["hmmfetch", "-o", output_file_name, hmm_file_name, acc]
    run_command(cmd)


class HMMCollection(object):

    def __init__(self, big_hmm_file, accession_ids):
        self.__by_accession = {}

        accession_ids = list(set(accession_ids))
        with tempfile.TemporaryDirectory(delete=False) as temp_dir:
            self.__temp_dir = temp_dir
            for acc in accession_ids:
                with tempfile.NamedTemporaryFile(dir=temp_dir, delete=False, suffix=".hmm") as tmpf:
                    tmpf.close()
                    print("fetching hmm into temp file", tmpf.name)
                    hmmfetch(big_hmm_file, acc, tmpf.name)
                    self.__by_accession[acc] = tmpf.name

    def get(self, accession):
        return self.__by_accession[accession]


    def clean(self):
        for fn in self.__by_accession.values():
            print("removing", fn)
            os.remove(fn)
        print("removing", self.__temp_dir)
        shutil.rmtree(self.__temp_dir)
        self.__by_accession = None
        self.__temp_dir = None
