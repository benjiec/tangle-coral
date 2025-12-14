import unittest

from needle.seq import extract_subsequence, extract_subsequence_strand_sensitive, compute_three_frame_translations


class TestExtractSubsequence(unittest.TestCase):

    def test_forward_in_bounds(self):
        seq = "ABCDEFGH"  # 1..8
        # 2..5 => BCDE
        self.assertEqual(extract_subsequence(seq, 2, 5), "BCDE")
        # 1..8 => full
        self.assertEqual(extract_subsequence(seq, 1, 8), "ABCDEFGH")

    def test_reversed_coords_same_slice(self):
        seq = "ABCDEFGH"
        # reversed coordinates (5..2) should return same slice as (2..5)
        self.assertEqual(extract_subsequence(seq, 5, 2), "BCDE")

    def test_end_beyond_length_is_clamped(self):
        seq = "ABCDEFGH"  # len=8
        # 6..12 => 6..8 => FGH
        self.assertEqual(extract_subsequence(seq, 6, 12), "FGH")
        # reversed with beyond-length
        self.assertEqual(extract_subsequence(seq, 12, 6), "FGH")

    def test_start_beyond_length_returns_none(self):
        seq = "ABCDEFGH"  # len=8
        # left (min) > len(seq) -> None
        self.assertIsNone(extract_subsequence(seq, 9, 10))
        self.assertIsNone(extract_subsequence(seq, 10, 9))

    def test_invalid_coords_non_positive(self):
        seq = "ABCDEFGH"
        self.assertIsNone(extract_subsequence(seq, 0, 5))
        self.assertIsNone(extract_subsequence(seq, 5, 0))
        self.assertIsNone(extract_subsequence(seq, -1, 2))

    def test_none_full_sequence_returns_none(self):
        self.assertIsNone(extract_subsequence(None, 1, 3))

    def test_empty_sequence_returns_none(self):
        # len("") == 0, left=min(1,1)=1 > 0 -> None
        self.assertIsNone(extract_subsequence("", 1, 1))

    def test_reversed_coords_strand_sensitie(self):
        seq = "TAGTCAAA"
        # reversed coordinates (5..2) should return same slice as (2..5)
        self.assertEqual(extract_subsequence_strand_sensitive(seq, 5, 2), "GACT")


class TestThreeFrameTranslations(unittest.TestCase):

    def test_three_frame_translation_computes_per_frame_sequence_and_dna_coordinates(self):
        genomic = "ATGCGATGACTTCGTTATGCTT"

        # fwd
        self.assertEqual(
          compute_three_frame_translations(genomic, 2, 17),
          [
            (2, 16, "CDDFV"),
            (3, 17, "AMTSL"),
            (4, 15, "R*LR"),
          ]
        )

        # rev
        # TGCGATGACTTCGTT -> AACGAAGTCATCGCA
        self.assertEqual(
          compute_three_frame_translations(genomic, 16, 2),
          [
            (16, 2, "NEVIA"),
            (15, 4, "TKSS"),
            (14, 3, "RSHR")
          ]
        )
