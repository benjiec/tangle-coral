from sieve.rules import HMMAlignment, Leader, Pfam, Rules, Sequence

βTrCP_pattern = r"[DE][ST]G[LIVMFYWA].{1,2}[ST]"

ank_rule = Rules(
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
    & Sequence.matches_regex(βTrCP_pattern).relativeToPfam("PF12796", -180, -10)
)

# BCL3 requires localization=nucleus by DeepLoc
# also BCL3 does not have DSG..S
