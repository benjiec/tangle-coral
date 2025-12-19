import os
import tempfile
import unittest
import io

from needle.io import (
    write_nucmatch_rows,
    write_fasta_record,
    export_protein_hits,
)
import needle.io as io_mod

from needle.match import ProteinHit, Match


class TestIO(unittest.TestCase):

    def test_write_nucmatch_rows(self):
        a = Match("QX","TX",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        b = Match("QX","TX",4,6,10,18,0.0,90.0,False); b.target_sequence="GAAGTGGGG"
        pm = ProteinHit([a,b],1,6,1,18)
        pid = pm.protein_hit_id
        buf = io.StringIO()
        write_nucmatch_rows(buf, pm, "G", 1.1, 2.2)
        rows = [r for r in buf.getvalue().splitlines() if r]
        self.assertEqual(len(rows), 2)
        cols0 = rows[0].split("\t")
        self.assertEqual(cols0[0], pid)
        self.assertEqual(cols0[1], "G")
        self.assertEqual(cols0[2], "TX")
        self.assertEqual(cols0[3], "1")
        self.assertEqual(cols0[4], "9")
        self.assertEqual(cols0[5], "QX")
        self.assertEqual(cols0[6], "1")
        self.assertEqual(cols0[7], "3")
        self.assertEqual(cols0[8], "1.1")
        self.assertEqual(cols0[9], "2.2")
        cols1 = rows[1].split("\t")
        self.assertEqual(cols1[0], pid)
        self.assertEqual(cols0[1], "G")
        self.assertEqual(cols0[2], "TX")
        self.assertEqual(cols1[3], "10")
        self.assertEqual(cols1[4], "18")
        self.assertEqual(cols0[5], "QX")
        self.assertEqual(cols1[6], "4")
        self.assertEqual(cols1[7], "6")
        self.assertEqual(cols0[8], "1.1")
        self.assertEqual(cols0[9], "2.2")

    def test_write_fasta_record(self):
        a = Match("QX","TX",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        pm = ProteinHit([a],1,3,1,9)
        pid = pm.protein_hit_id
        buf = io.StringIO()
        write_fasta_record(buf, pm)
        s = buf.getvalue().splitlines()
        self.assertEqual(s[0], f">{pid}")
        self.assertTrue(s[1].startswith("MEF"))

    def test_export_protein_hits_filters_and_writes(self):
        # pm1: single block => eligible
        a = Match("Q1","T1",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        pm1 = ProteinHit([a],1,3,1,9, hmm_file="ignored.hmm")
        # pm2: overlapping blocks => not single sequence
        b1 = Match("Q2","T2",1,3,1,9,0.0,100.0,False); b1.target_sequence="ATGGAATTT"
        b2 = Match("Q2","T2",3,5,10,18,0.0,100.0,False); b2.target_sequence="GAAGTGGGG"
        pm2 = ProteinHit([b1,b2],1,5,1,18, hmm_file="ignored.hmm")
        with tempfile.TemporaryDirectory() as d:
            p1 = os.path.join(d, "prot.faa")
            p2 = os.path.join(d, "nuc.tsv")
            p3 = os.path.join(d, "prot_fastas")
            orig = io_mod.hmmsearch
            try:
                io_mod.hmmsearch = lambda hmm, seqs, gap_removal: [dict(seq_score=55.0, seq_evalue=1e-6)]
                export_protein_hits("GENOMEZ", [pm1, pm2], p1, p2, p3)
            finally:
                io_mod.hmmsearch = orig
            with open(p1) as f:
                lines = [l.strip() for l in f if l.strip()]
            with open(p2) as f:
                lines2 = [l.strip() for l in f if l.strip()]
            faa_path = os.path.join(p3, "Q1.faa")
            with open(faa_path) as f:
                lines3 = [l.strip() for l in f if l.strip()]
            # fasta: 2 lines (header+seq)
            self.assertEqual(len(lines), 2)
            # nuc rows: only pm1 yields rows => 1 line + header
            self.assertEqual(len(lines2), 2)
            # fasta: 2 lines (header+seq)
            self.assertEqual(len(lines3), 2)

    def test_export_protein_hits_appends_not_overwrites(self):
        a = Match("QX","TX",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        pm = ProteinHit([a],1,3,1,9, hmm_file="ignored.hmm")
        with tempfile.TemporaryDirectory() as d:
            p1 = os.path.join(d, "prot.faa")
            p2 = os.path.join(d, "nuc.tsv")
            p3 = os.path.join(d, "fastas")
            orig = io_mod.hmmsearch
            try:
                io_mod.hmmsearch = lambda hmm, seqs, gap_removal: [dict(seq_score=10.0, seq_evalue=2e-5)]
                export_protein_hits("GENOME1", [pm], p1, p2, p3)
                # capture initial counts
                with open(p1) as f:
                    lines_prot_1 = [l for l in f if l.strip()]
                with open(p2) as f:
                    lines_nuc_1 = [l for l in f if l.strip()]
                faa_path = os.path.join(p3, "QX.faa")
                with open(faa_path) as f:
                    lines_faa_1 = [l for l in f if l.strip()]
                # second export appends
                export_protein_hits("GENOME1", [pm], p1, p2, p3)
                with open(p1) as f:
                    lines_prot_2 = [l for l in f if l.strip()]
                with open(p2) as f:
                    lines_nuc_2 = [l for l in f if l.strip()]
                with open(faa_path) as f:
                    lines_faa_2 = [l for l in f if l.strip()]
            finally:
                io_mod.hmmsearch = orig
            # proteins.tsv: header + 1 row initially; after append expect +2 lines
            self.assertEqual(len(lines_prot_1), 2)
            self.assertEqual(len(lines_prot_2), 4)
            # nuc.tsv: header + 1 row initially; after append expect +1 row
            self.assertEqual(len(lines_nuc_1), 2)
            self.assertEqual(len(lines_nuc_2), 3)
            # faa file: 2 lines initially; after append expect +2 lines
            self.assertEqual(len(lines_faa_1), 2)
            self.assertEqual(len(lines_faa_2), 4)


if __name__ == "__main__":
    unittest.main()
