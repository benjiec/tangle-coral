from sieve.rules import HMMAlignment, Leader, Pfam, Rules, Sequence

βTrCP_pattern = r"[DE][ST]G[LIVMFYWA].{1,2}[ST]"

inhibitor_rule = Rules(
    # anything that matches Ankyrin repeat (3 copies)
    Pfam.matches("PF12796")

    # a few known negatives
    & ~Pfam.matches("PF00554")
    & ~Pfam.matches("PF08424")
    & ~Pfam.matches("PF01876")
)

canonical_ikb_rule = Rules(
    # anything that matches Ankyrin repeat (3 copies)
    Pfam.matches("PF12796")
    & Sequence.length_between(300, 1200)

    # a few known negatives
    & ~Pfam.matches("PF00554")
    & ~Pfam.matches("PF08424")
    & ~Pfam.matches("PF01876")

    # required IKK phosphorylation and then degradation by βTrCP
    & Sequence.matches_regex(βTrCP_pattern).relativeToPfam("PF12796", -100, -20)
)

nucleus_inhinitor_rule = Rules(
    # anything that matches Ankyrin repeat (3 copies)
    Pfam.matches("PF12796")

    # a few known negatives
    & ~Pfam.matches("PF00554")
    & ~Pfam.matches("PF08424")
    & ~Pfam.matches("PF01876")

    # NLS
    & (Sequence.matches_regex(
         r"[KR]{3,5}|P[KR]{3,4}"
       ).relativeToPfam("PF12796", -90, -20) |
       Sequence.matches_regex(
         r"[KR]{2}.{10,14}[KR]{3,5}"
       ).relativeToPfam("PF12796", -90, -20))
)

# also BCL3 does not have DSG..S
