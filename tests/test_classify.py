import os
import csv
import tempfile
import unittest

from needle.classify import ClassifyTSV, classify
from needle.match import ProteinHit, Match, ProteinsTSV
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
                 seq_score = 20,
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

            self.assertEqual(rows[0]["dom_rank_for_protein"], "2")
            self.assertEqual(rows[1]["dom_rank_for_protein"], "1")
            self.assertEqual(rows[2]["dom_rank_for_protein"], "1")

            os.remove(tmpf.name)
