from sieve.rules import HMMAlignment, Leader, Pfam, Rules, Sequence

ikb_rule = Rules(
    # anything that matches Ankyrin repeat (3 copies)
    Pfam.matches("PF12796")
    & Sequence.length_between(300, 1200)

    # a few known negatives
    & ~Pfam.matches("PF00554")
    & ~Pfam.matches("PF08424")
    & ~Pfam.matches("PF01876")

    # required IKK phosphorylation and then degradation by βTrCP
    & Sequence.matches_regex(r"[DE][ST]G[LIVMFYWA].{1,2}[ST]").relativeToPfam("PF12796", -100, -20)
)

bcl3_rule = Rules(
    # anything that matches Ankyrin repeat (3 copies)
    Pfam.matches("PF12796")
    & Sequence.length_between(300, 1200)

    # a few known negatives
    & ~Pfam.matches("PF00554")
    & ~Pfam.matches("PF08424")
    & ~Pfam.matches("PF01876")

    # does not have a DSG..S site
    & ~Sequence.matches_regex(r"[DE][ST]G[LIVMFYWA].{1,2}[ST]").relativeToPfam("PF12796", -100, -20)

    # NLS
    & (Sequence.matches_regex(
         r"[KR]{3,5}|P[KR]{3,4}"
       ).relativeToPfam("PF12796", -90, -20) |
       Sequence.matches_regex(
         r"[KR]{2}.{10,14}[KR]{3,5}"
       ).relativeToPfam("PF12796", -90, -20))
)
