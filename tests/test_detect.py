import os
import tempfile
import unittest
from needle.detect import Results, hmm_search_genome_sequence
import needle.detect as detect_mod


class TestParseDetectResults(unittest.TestCase):

    def test_parse_ncbi_header_and_extract_sequences_reverse(self):
        """Verify reverse-direction hit: target sequence normalized to 5'->3' and reverse-complemented for translation."""

        # Create synthetic FASTA and TSV with NCBI-style headers; reverse-direction on target (sstart > send)
        with tempfile.TemporaryDirectory() as tmpdir:
            target_fasta_path = os.path.join(tmpdir, "target.faa")
            results_tsv_path = os.path.join(tmpdir, "results.tsv")

            # >S1 DNA: AAACCCGGGTTT
            with open(target_fasta_path, "w") as f:
                f.write(">S1 synthetic subject\n")
                f.write("AAACCCGGGTTT\n")

            header = "\t".join(Results.PRODUCER_HEADER)
            # qseqid sseqid evalue pident qstart qend sstart send sseq
            row = ["Q1", "S1", "1e-5", "99.9", "2", "5", "10", "7", "ABCD"]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(row) + "\n")

            res = Results(results_tsv_path, target_fasta_path)
            matches = res.matches()
            self.assertEqual(len(matches), 1)
            m = matches[0]
            self.assertEqual(m.query_accession, "Q1")
            self.assertEqual(m.target_accession, "S1")
            self.assertAlmostEqual(m.e_value, 1e-5, places=12)
            self.assertAlmostEqual(m.identity, 99.9, places=3)
            self.assertEqual(m.query_start, 2)
            self.assertEqual(m.query_end, 5)
            self.assertEqual(m.target_start, 10)
            self.assertEqual(m.target_end, 7)

            # Target 7..10 on AAACCCGGGTTT -> "GGGT"; reverse-complement -> "ACCC"
            self.assertEqual(m.target_sequence, "ACCC")
            # Ensure translated target (revcomp path) is a valid amino-acid sequence length <= query segment
            _ = m.target_sequence_translated()
            # Matched sequence preserved
            self.assertEqual(m.matched_sequence, "ABCD")

    def test_parse_ncbi_header_and_extract_sequences_forward(self):
        """Verify forward-direction hit: target sequence 5'->3' matches the query AA after translation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            target_fasta_path = os.path.join(tmpdir, "t.fna")
            results_tsv_path = os.path.join(tmpdir, "r.tsv")

            # Target DNA forward: ATG GAA TTT -> MEF at positions 5..13
            with open(target_fasta_path, "w") as f:
                f.write(">Tfwd\n")
                f.write("NNNNATGGAATTTNNNN\n")  # 1..4 N; 5..13 coding; 14..17 N

            header = "\t".join(Results.PRODUCER_HEADER)
            # Forward direction: sstart < send
            row = ["Qfwd", "Tfwd", "0", "100.0", "1", "3", "5", "13", "MEF"]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(row) + "\n")

            res = Results(results_tsv_path, target_fasta_path)
            ms = res.matches()
            self.assertEqual(len(ms), 1)
            m = ms[0]
            # Target DNA 5..13 5'->3'
            self.assertEqual(m.target_sequence, "ATGGAATTT")
            # Translated target equals matched_sequence equals query segment
            self.assertEqual(m.target_sequence_translated(), "MEF")
            self.assertEqual(m.matched_sequence, "MEF")


class TestHMMSearchGenomeSequence(unittest.TestCase):

    def setUp(self):
        self.orig = detect_mod.hmmsearch_sequence_dict
        self.hmmsearch_calls_translated = []

        def _fake_hmmsearch(hmm_file, translated): 
            self.hmmsearch_calls_translated.append(translated)
            return [
                dict(target_name = k,
                     target_accession = "t1",
                     query_name = "q",
                     query_accession = "q1",
                     seq_evalue = 0.001,
                     seq_score = 20,
                     dom_evalue = 0.001,
                     dom_score = 20,
                     query_length = 60,
                     hmm_from = 3,
                     hmm_to = 6,
                     target_length = 7,
                     ali_from = 2,
                     ali_to = 5 
                ) for k,v in translated.items()
            ]

        detect_mod.hmmsearch_sequence_dict = _fake_hmmsearch

    def tearDown(self):
        detect_mod.hmmsearch_sequence_dict = self.orig

    def test_hmm_search_genome_sequence_windows_correctly_and_returns_right_coordinates(self):
        target_accession = "t1"
        target_sequence = "A"*20+"T"*20+"G"*20
        win = 20
        win_overlap = 3
        contig_start = None
        contig_end = None

        detected = hmm_search_genome_sequence(None, target_accession, target_sequence, win, win_overlap, contig_start, contig_end)

        # calls hmmsearch_sequence_dict once for every window, there are 3 windows, each with 6 translations
        self.assertEqual(len(self.hmmsearch_calls_translated), 3)
        self.assertEqual(len(self.hmmsearch_calls_translated[0]), 6)
        self.assertEqual(len(self.hmmsearch_calls_translated[1]), 6)
        self.assertEqual(len(self.hmmsearch_calls_translated[2]), 6)

        # translations are correct
        # 3 windows: 1..20+3, 21..40+3, 41..60
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_fwd_0"], "KKKKKKN")
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_fwd_1"], "KKKKKKI")
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_fwd_2"], "KKKKKKF")
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_rev_0"], "KFFFFFF")
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_rev_1"], "NFFFFFF")
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_rev_2"], "IFFFFFF")
        # check other two windows are correct at first frame
        self.assertEqual(self.hmmsearch_calls_translated[1]["t1_fwd_0"], "FFFFFFL")
        self.assertEqual(self.hmmsearch_calls_translated[2]["t1_fwd_0"], "GGGGGG")

        # converted returned coordinates to dna coordinates by each frame, and uses Results.PRODUCER_HEADER format
        # 3 windows: 1..20+3, 21..40+3, 41..60
        self.assertEqual(len(detected), 6*3)
        # window 1..20+3, frame 0 - see setUp(): matched query start and end is 3..6, ali is 2..5
        # ali coordinate 2..5 is aa coordinate, and should be translated to dna coordinate based on frame+window
        self.assertEqual(detected[0][Results.PRODUCER_HEADER.index(Results.H_QSEQID)], "q1")
        self.assertEqual(detected[0][Results.PRODUCER_HEADER.index(Results.H_SSEQID)], "t1")
        self.assertEqual(detected[0][Results.PRODUCER_HEADER.index(Results.H_EVALUE)], 0.001) 
        self.assertEqual(detected[0][Results.PRODUCER_HEADER.index(Results.H_QSTART)], 3)
        self.assertEqual(detected[0][Results.PRODUCER_HEADER.index(Results.H_QEND)], 6)
        self.assertEqual(detected[0][Results.PRODUCER_HEADER.index(Results.H_SSTART)], 4)
        self.assertEqual(detected[0][Results.PRODUCER_HEADER.index(Results.H_SEND)], 15)
        # two more forward frames
        self.assertEqual(detected[1][Results.PRODUCER_HEADER.index(Results.H_QSEQID)], "q1")
        self.assertEqual(detected[1][Results.PRODUCER_HEADER.index(Results.H_SSEQID)], "t1")
        self.assertEqual(detected[1][Results.PRODUCER_HEADER.index(Results.H_EVALUE)], 0.001) 
        self.assertEqual(detected[1][Results.PRODUCER_HEADER.index(Results.H_QSTART)], 3)
        self.assertEqual(detected[1][Results.PRODUCER_HEADER.index(Results.H_QEND)], 6)
        self.assertEqual(detected[1][Results.PRODUCER_HEADER.index(Results.H_SSTART)], 5)
        self.assertEqual(detected[1][Results.PRODUCER_HEADER.index(Results.H_SEND)], 16)
        self.assertEqual(detected[2][Results.PRODUCER_HEADER.index(Results.H_SSTART)], 6)
        self.assertEqual(detected[2][Results.PRODUCER_HEADER.index(Results.H_SEND)], 17)
        # three reverse frames, starting at 23
        self.assertEqual(detected[3][Results.PRODUCER_HEADER.index(Results.H_SSTART)], 20)
        self.assertEqual(detected[3][Results.PRODUCER_HEADER.index(Results.H_SEND)], 9)
        self.assertEqual(detected[4][Results.PRODUCER_HEADER.index(Results.H_SSTART)], 19)
        self.assertEqual(detected[4][Results.PRODUCER_HEADER.index(Results.H_SEND)], 8)
        self.assertEqual(detected[5][Results.PRODUCER_HEADER.index(Results.H_SSTART)], 18)
        self.assertEqual(detected[5][Results.PRODUCER_HEADER.index(Results.H_SEND)], 7)
        # next window, fwd
        self.assertEqual(detected[6][Results.PRODUCER_HEADER.index(Results.H_SSTART)], 24)
        self.assertEqual(detected[6][Results.PRODUCER_HEADER.index(Results.H_SEND)], 35)
        self.assertEqual(detected[7][Results.PRODUCER_HEADER.index(Results.H_SSTART)], 25)
        self.assertEqual(detected[7][Results.PRODUCER_HEADER.index(Results.H_SEND)], 36)
        self.assertEqual(detected[8][Results.PRODUCER_HEADER.index(Results.H_SSTART)], 26)
        self.assertEqual(detected[8][Results.PRODUCER_HEADER.index(Results.H_SEND)], 37)

    def test_hmm_search_genome_sequence_windows_honest_contig_coord_window(self):
        target_accession = "t1"
        target_sequence = "A"*20+"T"*20+"G"*20
        win = 20
        win_overlap = 3
        contig_left_0b = 1
        contig_right_excl_0b = 23

        detected = hmm_search_genome_sequence(None, target_accession, target_sequence, win, win_overlap, contig_left_0b, contig_right_excl_0b)

        # full contig is 2..23 (or 1..22 in 0b space)
        # 3 windows: 2..min(21+3,23), 22..23
        # calls hmmsearch_sequence_dict once for every window, there are 2 windows, but second window is not big enough so no translations
        self.assertEqual(len(self.hmmsearch_calls_translated), 1)
        self.assertEqual(len(self.hmmsearch_calls_translated[0]), 6)

        # translations are correct
        # 2..min(21+3,23) or 2..23, last codon is 20..22
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_fwd_0"], "KKKKKKI")
        # 3..23, last codon is 21..23
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_fwd_1"], "KKKKKKF")
        # 4..23, last codon is 19..21
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_fwd_2"], "KKKKKN")
        # on rev, first codon is 23..21
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_rev_0"], "KFFFFFF")
        # on rev, first codon is 22..20
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_rev_1"], "NFFFFFF")
        # on rev, first codon is 21..19
        self.assertEqual(self.hmmsearch_calls_translated[0]["t1_rev_2"], "IFFFFF")
