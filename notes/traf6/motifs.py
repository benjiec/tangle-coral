import re
import sys
import argparse
from tangle import unique_batch
from tangle.sequence import read_fasta_as_dict
from tangle.detected import DetectedTable

MOTIFS = (
  (r"P.E..[FYWDE]", "traf6_binding_strict"),
  (r"P.[QE]..[FYWHDEILVKR]", "traf6_binding_permissive"),
  (r"P.[QE].[ST]", "traf2_binding_strict"),
  (r"P.[QEN].[STEDAG]", "traf2_binding_variant"),
  (r"P.{1,2}[QEN].{1,2}[STEDAG]", "traf2_binding_permissive"),
)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args(argv)

    sequences = read_fasta_as_dict(args.input)
    batch = unique_batch()

    motifs = [(re.compile(pat), motif) for pat, motif in MOTIFS]
    rows = []

    for acc,seq in sequences.items():
        for pat,motif in motifs:
            matches = [(m.start(), m.end(), m.group()) for m in pat.finditer(seq)]
            for start, end, matched in matches:
                row = dict(
                  detection_type="sequence",
                  detection_method="other",
                  batch=batch,
                  query_accession=acc,
                  query_database="_",
                  target_accession=motif,
                  target_database="Motif",
                  target_model=matched,
                  query_start=start,
                  query_end=end,
                  target_start=1,
                  target_end=len(matched)
                )
                rows.append(row)

    DetectedTable.write_tsv(args.output, rows)


if __name__ == "__main__":
    sys.exit(main())
