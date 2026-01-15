import os
import csv
import tempfile
import unittest

from needle.classify import ClassifyTSV, classify, group_by_assignment, assign_ko
from needle.match import ProteinHit, Match, ProteinsTSV
from needle.seq import write_fasta_from_dict, read_fasta_as_dict
import needle.classify as classify_mod


class TestClassifyTSV(unittest.TestCase):

    def setUp(self):
        self.orig_hmmscan_file = classify_mod.hmmscan_file
        self.hmmscan_called_args = []
        self.fake_hmm_rows = []

        def _fake_hmmscan(*args, **kwargs):
            self.hmmscan_called_args.append(list(args)+[kwargs])
            return self.fake_hmm_rows

        classify_mod.hmmscan_file = _fake_hmmscan

        m1 = Match(query_accession="q1", target_accession="t1",
                   query_start=2, query_end=5, target_start=10001, target_end=10012,
                   e_value=1e-5, identity=0.9, target_sequence="A"*12)
        p1 = ProteinHit(matches=[m1], query_start=2, query_end=5, target_start=10001, target_end=10012)
        m2 = Match(query_accession="q2", target_accession="t2",
                   query_start=2, query_end=5, target_start=10001, target_end=10012,
                   e_value=1e-5, identity=0.9, target_sequence="A"*12)
        p2 = ProteinHit(matches=[m2], query_start=2, query_end=5, target_start=10001, target_end=10012)
        self.proteins = [p1, p2]
        self.protein_genome_accession_dict = { p1.protein_hit_id: "g1", p2.protein_hit_id: "g2" }
        self.set_default_fake_hmm_rows()

    def tearDown(self):
        classify_mod.hmmscan_file = self.orig_hmmscan_file

    def set_default_fake_hmm_rows(self):
        assert len(self.proteins) == 2
        self.fake_hmm_rows = [
            dict(target_name = "h1",
                 target_accession = "h1",
                 target_length = 51,
                 query_name = self.proteins[0].protein_hit_id,
                 query_accession = self.proteins[0].protein_hit_id,
                 query_length = 61,
                 seq_evalue = 0.001,
                 seq_score = 30,
                 dom_evalue = 0.002,
                 dom_evalue_cond = 0.003,
                 dom_score = 21,
                 hmm_from = 3,
                 hmm_to = 6,
                 ali_from = 2,
                 ali_to = 5 
            ),
            dict(target_name = "h2",
                 target_accession = "h2",
                 target_length = 52,
                 query_name = self.proteins[0].protein_hit_id,
                 query_accession = self.proteins[0].protein_hit_id,
                 query_length = 62,
                 seq_evalue = 0.0001,
                 seq_score = 29,
                 dom_evalue = 0.0002,
                 dom_evalue_cond = 0.0003,
                 dom_score = 28,
                 hmm_from = 3,
                 hmm_to = 6,
                 ali_from = 2,
                 ali_to = 5 
            ),
            dict(target_name = "h3",
                 target_accession = "h3",
                 target_length = 53,
                 query_name = self.proteins[1].protein_hit_id,
                 query_accession = self.proteins[1].protein_hit_id,
                 query_length = 63,
                 seq_evalue = 0.00001,
                 seq_score = 23,
                 dom_evalue = 0.00002,
                 dom_evalue_cond = 0.00003,
                 dom_score = 24,
                 hmm_from = 3,
                 hmm_to = 6,
                 ali_from = 2,
                 ali_to = 5 
            )
        ]

    def test_classify_calls_hmmscan_correctly(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv", mode="w") as tmpf:
            tmpf.close()
            classify("fake_hmm_file_name", "fake_proteins_faa_name", True, tmpf.name, self.protein_genome_accession_dict, None)
            self.assertEqual(self.hmmscan_called_args[0], ["fake_hmm_file_name", "fake_proteins_faa_name", {"cutoff": True}])
            classify("fake_hmm_file_name", "fake_proteins_faa_name", False, tmpf.name, self.protein_genome_accession_dict, None)
            self.assertEqual(self.hmmscan_called_args[1], ["fake_hmm_file_name", "fake_proteins_faa_name", {"cutoff": False}])
            os.remove(tmpf.name)

    def read_classify_tsv_outputs(self, fn):
        with open(fn, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = [r for r in reader]
            return rows

    def test_classify_tsv_includes_genome_accession_from_def_file(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv", mode="w") as tmpf:
            tmpf.close()

            classify("fake_hmm_file_name", "fake_proteins_faa_name", True, tmpf.name, self.protein_genome_accession_dict, None)
            rows = self.read_classify_tsv_outputs(tmpf.name)
            self.assertEqual(len(rows), 3)

            self.assertEqual(rows[0]["protein_accession"], self.proteins[0].protein_hit_id)
            self.assertEqual(rows[0]["genome_accession"], "g1")
            self.assertEqual(rows[1]["protein_accession"], self.proteins[0].protein_hit_id)
            self.assertEqual(rows[1]["genome_accession"], "g1")
            self.assertEqual(rows[2]["protein_accession"], self.proteins[1].protein_hit_id)
            self.assertEqual(rows[2]["genome_accession"], "g2")

            os.remove(tmpf.name)

    def test_classify_tsv_appends_to_output_file(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv", mode="w") as tmpf:
            tmpf.close()

            classify("fake_hmm_file_name", "fake_proteins_faa_name", True, tmpf.name, self.protein_genome_accession_dict, None)
            rows = self.read_classify_tsv_outputs(tmpf.name)
            self.assertEqual(len(rows), 3)

            classify("fake_hmm_file_name", "fake_proteins_faa_name", True, tmpf.name, self.protein_genome_accession_dict, None)
            rows = self.read_classify_tsv_outputs(tmpf.name)
            self.assertEqual(len(rows), 6)

            os.remove(tmpf.name)

    def test_includes_evalue_from_hmmscan_and_sets_target_and_query_correctly(self):
        # hmmscan returns both accession and name

        with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv", mode="w") as tmpf:
            tmpf.close()

            classify("fake_hmm_file_name", "fake_proteins_faa_name", True, tmpf.name, self.protein_genome_accession_dict, None)
            rows = self.read_classify_tsv_outputs(tmpf.name)
            self.assertEqual(len(rows), 3)

            self.assertEqual(rows[0]["dom_evalue"], str(self.fake_hmm_rows[0]["dom_evalue"]))
            self.assertEqual(rows[0]["dom_evalue_cond"], str(self.fake_hmm_rows[0]["dom_evalue_cond"]))
            self.assertEqual(rows[0]["dom_score"], str(self.fake_hmm_rows[0]["dom_score"]))
            self.assertEqual(rows[0]["protein_accession"], self.fake_hmm_rows[0]["query_accession"])
            self.assertEqual(rows[0]["protein_start"], str(self.fake_hmm_rows[0]["ali_from"]))
            self.assertEqual(rows[0]["protein_end"], str(self.fake_hmm_rows[0]["ali_to"]))
            self.assertEqual(rows[0]["hmm_accession"], self.fake_hmm_rows[0]["target_accession"])
            self.assertEqual(rows[0]["hmm_start"], str(self.fake_hmm_rows[0]["hmm_from"]))
            self.assertEqual(rows[0]["hmm_end"], str(self.fake_hmm_rows[0]["hmm_to"]))

            self.assertEqual(rows[1]["dom_evalue"], str(self.fake_hmm_rows[1]["dom_evalue"]))
            self.assertEqual(rows[1]["dom_evalue_cond"], str(self.fake_hmm_rows[1]["dom_evalue_cond"]))
            self.assertEqual(rows[1]["dom_score"], str(self.fake_hmm_rows[1]["dom_score"]))
            self.assertEqual(rows[1]["protein_accession"], self.fake_hmm_rows[1]["query_accession"])
            self.assertEqual(rows[1]["protein_start"], str(self.fake_hmm_rows[1]["ali_from"]))
            self.assertEqual(rows[1]["protein_end"], str(self.fake_hmm_rows[1]["ali_to"]))
            self.assertEqual(rows[1]["hmm_accession"], self.fake_hmm_rows[1]["target_accession"])
            self.assertEqual(rows[1]["hmm_start"], str(self.fake_hmm_rows[1]["hmm_from"]))
            self.assertEqual(rows[1]["hmm_end"], str(self.fake_hmm_rows[1]["hmm_to"]))

            self.assertEqual(rows[2]["dom_evalue"], str(self.fake_hmm_rows[2]["dom_evalue"]))
            self.assertEqual(rows[2]["dom_evalue_cond"], str(self.fake_hmm_rows[2]["dom_evalue_cond"]))
            self.assertEqual(rows[2]["dom_score"], str(self.fake_hmm_rows[2]["dom_score"]))
            self.assertEqual(rows[2]["protein_accession"], self.fake_hmm_rows[2]["query_accession"])
            self.assertEqual(rows[2]["protein_start"], str(self.fake_hmm_rows[2]["ali_from"]))
            self.assertEqual(rows[2]["protein_end"], str(self.fake_hmm_rows[2]["ali_to"]))
            self.assertEqual(rows[2]["hmm_accession"], self.fake_hmm_rows[2]["target_accession"])
            self.assertEqual(rows[2]["hmm_start"], str(self.fake_hmm_rows[2]["hmm_from"]))
            self.assertEqual(rows[2]["hmm_end"], str(self.fake_hmm_rows[2]["hmm_to"]))

            os.remove(tmpf.name)

        # without accession
        self.fake_hmm_rows[0]["query_accession"] = "-"
        self.fake_hmm_rows[0]["target_accession"] = "-"
        self.fake_hmm_rows[1]["query_accession"] = " "
        self.fake_hmm_rows[1]["target_accession"] = " "
        self.fake_hmm_rows[2]["query_accession"] = ""
        self.fake_hmm_rows[2]["target_accession"] = ""

        with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv", mode="w") as tmpf:
            tmpf.close()

            classify("fake_hmm_file_name", "fake_proteins_faa_name", True, tmpf.name, self.protein_genome_accession_dict, None)
            rows = self.read_classify_tsv_outputs(tmpf.name)
            self.assertEqual(len(rows), 3)

            self.assertEqual(rows[0]["protein_accession"], self.fake_hmm_rows[0]["query_name"])
            self.assertEqual(rows[0]["hmm_accession"], self.fake_hmm_rows[0]["target_name"])
            self.assertEqual(rows[1]["protein_accession"], self.fake_hmm_rows[1]["query_name"])
            self.assertEqual(rows[1]["hmm_accession"], self.fake_hmm_rows[1]["target_name"])
            self.assertEqual(rows[2]["protein_accession"], self.fake_hmm_rows[2]["query_name"])
            self.assertEqual(rows[2]["hmm_accession"], self.fake_hmm_rows[2]["target_name"])

            os.remove(tmpf.name)

    def test_includes_scoring_threshold(self):
        scoring_threshold = {"h1": 3, "h3": 4}

        with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv", mode="w") as tmpf:
            tmpf.close()

            classify("fake_hmm_file_name", "fake_proteins_faa_name", True, tmpf.name, self.protein_genome_accession_dict, scoring_threshold)
            rows = self.read_classify_tsv_outputs(tmpf.name)
            self.assertEqual(len(rows), 3)

            self.assertEqual(rows[0]["score_threshold"], "3")
            self.assertEqual(rows[1]["score_threshold"], "")
            self.assertEqual(rows[2]["score_threshold"], "4")

            os.remove(tmpf.name)

    def test_ranks_hits_for_protein(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv", mode="w") as tmpf:
            tmpf.close()

            classify("fake_hmm_file_name", "fake_proteins_faa_name", True, tmpf.name, self.protein_genome_accession_dict, None)
            rows = self.read_classify_tsv_outputs(tmpf.name)
            self.assertEqual(len(rows), 3)

            self.assertEqual(rows[0]["dom_rank_for_protein"], "2") # score 21 in first row vs 28 in second row
            self.assertEqual(rows[1]["dom_rank_for_protein"], "1") # score 21 in first row vs 28 in second row
            self.assertEqual(rows[2]["dom_rank_for_protein"], "1") # different protein

            os.remove(tmpf.name)


class TestAssignment(unittest.TestCase):

    def _classify_row(self, protein_id, hmm_db, hmm_acc, dom_score, score_threshold, dom_rank, protein_start = None, protein_end = None):
        return {
            ClassifyTSV.HDR_PROTEIN_ACCESSION: protein_id,
            ClassifyTSV.HDR_HMM_DB: hmm_db,
            ClassifyTSV.HDR_HMM_ACCESSION: hmm_acc,
            ClassifyTSV.HDR_HMM_START: 1,
            ClassifyTSV.HDR_HMM_END: 20,
            ClassifyTSV.HDR_PROTEIN_START: protein_start if protein_start else 1,
            ClassifyTSV.HDR_PROTEIN_END: protein_end if protein_end else 20,
            ClassifyTSV.HDR_DOM_SCORE: dom_score,
            ClassifyTSV.HDR_SCORE_THRESHOLD: score_threshold,
            ClassifyTSV.HDR_DOM_RANK: dom_rank
        }

    def test_group_by_assignment_groups_by_score_threshold(self):
        rows = [
            self._classify_row("p1", "ko", "k1", 10, 100, 2),    # excluded, rank != 1
            self._classify_row("p1", "ko", "k2", 11, 100, 1),    # excluded, score too low
            self._classify_row("p1", "ko", "k3", 11, 100, 1),    # excluded, score too low
            self._classify_row("p1", "pf", "pf1", 100, 100, 1),  # excluded, hmm db does not match
            self._classify_row("p2", "ko", "k1", 75, 100, 1),    # assigns k1 to p2
            self._classify_row("p2", "ko", "k2", 90, 100, 2),    # higher score, but rank != 1
            self._classify_row("p2", "pf0", "pf1", 100, 100, 1), # returned
            self._classify_row("p2", "pf0", "pf2", 10, 100, 2),  # returned
            self._classify_row("p3", "ko", "k2", 90, 100, 1),    # assigns k2 to p3
            self._classify_row("p4", "ko", "k1", 100, 100, 1),   # assigns k1 to p4
            self._classify_row("p4", "pf", "pf1", 100, 100, 1),  # returned
        ]

        assignments = group_by_assignment(rows, "ko", 0.75)
        assignments = {hmm_acc: list(group) for hmm_acc, group in assignments}

        # nothing assigned to k3, returns assigned KO#s
        self.assertEqual(sorted(list(assignments.keys())), ["k1", "k2"])

        self.assertEqual([(row["protein_accession"], row["hmm_accession"]) for row in assignments["k1"]],
                         [("p2", "k1"),
                          ("p2", "pf1"),
                          ("p2", "pf2"),
                          ("p4", "k1"),
                          ("p4", "pf1")])
        # even though p3 is assigned, no additional rows are returned
        self.assertEqual([(row["protein_accession"], row["hmm_accession"]) for row in assignments["k2"]],
                         [("p3", "k2")])

    def test_creates_protein_fasta_by_assignment(self):
        rows = [
            self._classify_row("p1", "ko", "k1", 10, 100, 2),    # excluded, rank != 1
            self._classify_row("p1", "ko", "k2", 11, 100, 1),    # excluded, score too low
            self._classify_row("p1", "ko", "k3", 11, 100, 1),    # excluded, score too low
            self._classify_row("p1", "pf", "pf1", 100, 100, 1),  # excluded, hmm db does not match
            self._classify_row("p2", "ko", "k1", 95, 100, 1),    # assigns k1 to p2
            self._classify_row("p2", "ko", "k2", 98, 100, 2),    # higher score, but rank != 1
            self._classify_row("p2", "pf0", "pf1", 100, 100, 1), # returned
            self._classify_row("p2", "pf0", "pf2", 10, 100, 2),  # returned
            self._classify_row("p3", "ko", "k2", 95, 100, 1),    # assigns k2 to p3
            self._classify_row("p4", "ko", "k1", 100, 100, 1),   # assigns k1 to p4
            self._classify_row("p4", "pf", "pf1", 100, 100, 1, 2, 18),  # returned
        ]

        proteins = {
            "p1": "A"*20,
            "p2": "K"*20,
            "p3": "L"*20,
            "p4": "LMK"*20,
        }

        # creates protein fasta
        with tempfile.TemporaryDirectory() as tmpdir:
            proteins_faa = tmpdir+"/proteins.faa"
            write_fasta_from_dict(proteins, proteins_faa)
            assign_ko(rows, "ko", proteins_faa, tmpdir)
            self.assertEqual(os.path.exists(tmpdir+"/k1.faa"), True)
            self.assertEqual(os.path.exists(tmpdir+"/k2.faa"), True)
            self.assertEqual(os.path.exists(tmpdir+"/k3.faa"), False)

            k1_seqs = read_fasta_as_dict(tmpdir+"/k1.faa")
            self.assertEqual(list(k1_seqs.keys()), ["p2", "p4"]) 
            self.assertEqual(k1_seqs["p2"], "K"*20)
            self.assertEqual(k1_seqs["p4"], "LMK"*20)

            k1_seqs = read_fasta_as_dict(tmpdir+"/k2.faa")
            self.assertEqual(list(k1_seqs.keys()), ["p3"])
            self.assertEqual(k1_seqs["p3"], "L"*20)

        # creates domain fasta with subsequence, and only from specified hmm db
        with tempfile.TemporaryDirectory() as tmpdir:
            proteins_faa = tmpdir+"/proteins.faa"
            write_fasta_from_dict(proteins, proteins_faa)
            assign_ko(rows, "ko", proteins_faa, tmpdir, "pf")
            self.assertEqual(os.path.exists(tmpdir+"/k1_pf1.faa"), True)
            self.assertEqual(os.path.exists(tmpdir+"/k1_pf2.faa"), False)
            self.assertEqual(os.path.exists(tmpdir+"/k2_pf1.faa"), False)
            self.assertEqual(os.path.exists(tmpdir+"/k2_pf2.faa"), False)

            k1_seqs = read_fasta_as_dict(tmpdir+"/k1_pf1.faa")
            self.assertEqual(list(k1_seqs.keys()), ["p4_2_18_pf1"])
            self.assertEqual(k1_seqs["p4_2_18_pf1"], "MKLMKLMKLMKLMKLMK")

    def test_appends_to_existing_fasta(self):
        rows = [
            self._classify_row("p1", "ko", "k1", 10, 100, 2),    # excluded, rank != 1
            self._classify_row("p1", "ko", "k2", 11, 100, 1),    # excluded, score too low
            self._classify_row("p1", "ko", "k3", 11, 100, 1),    # excluded, score too low
            self._classify_row("p1", "pf", "pf1", 100, 100, 1),  # excluded, hmm db does not match
            self._classify_row("p2", "ko", "k1", 95, 100, 1),    # assigns k1 to p2
            self._classify_row("p2", "ko", "k2", 98, 100, 2),    # higher score, but rank != 1
            self._classify_row("p2", "pf0", "pf1", 100, 100, 1), # returned
            self._classify_row("p2", "pf0", "pf2", 10, 100, 2),  # returned
            self._classify_row("p3", "ko", "k2", 95, 100, 1),    # assigns k2 to p3
            self._classify_row("p4", "ko", "k1", 100, 100, 1),   # assigns k1 to p4
            self._classify_row("p4", "pf", "pf1", 100, 100, 1, 2, 18),  # returned
        ]

        proteins = {
            "p1": "A"*20,
            "p2": "K"*20,
            "p3": "L"*20,
            "p4": "LMK"*20,
        }

        # creates protein fasta, only p2 and p4 assigned to k1
        with tempfile.TemporaryDirectory() as tmpdir:
            proteins_faa = tmpdir+"/proteins.faa"
            write_fasta_from_dict(proteins, proteins_faa)
            assign_ko(rows, "ko", proteins_faa, tmpdir)
            self.assertEqual(os.path.exists(tmpdir+"/k1.faa"), True)
            k1_seqs = read_fasta_as_dict(tmpdir+"/k1.faa")
            self.assertEqual(list(k1_seqs.keys()), ["p2", "p4"]) 

        # does not override existing fasta content
        with tempfile.TemporaryDirectory() as tmpdir:
            proteins_faa = tmpdir+"/proteins.faa"
            write_fasta_from_dict(proteins, proteins_faa)
            some_existing_entries = { "x1": "A"*20 }
            write_fasta_from_dict(some_existing_entries, tmpdir+"/k1.faa")
            assign_ko(rows, "ko", proteins_faa, tmpdir)
            self.assertEqual(os.path.exists(tmpdir+"/k1.faa"), True)
            k1_seqs = read_fasta_as_dict(tmpdir+"/k1.faa")
            self.assertEqual(list(k1_seqs.keys()), ["x1", "p2", "p4"]) 
