from sieve.rules import HMMAlignment, Leader, Pfam, Rules, Sequence

# lysine (K) within 40 AAs upstream of a DSG(\phi).S degron motif
βTrCP_pattern = r"K.{1,40}?[DEST][ST]G[LIVMFYWA].{1,2}[ST]"


# all encompassing rule to filter sequences so we can send them to DeepLoc
no_deeploc_rule = Rules(

    Pfam.matches("PF00023").times(3, 11)
    & Sequence.length_between(300, 1000)
    & Pfam.matches_only("PF00023", "PF12796", "PF13606", "PF13637", "PF13857")

    & Pfam.matches("PF00023").betweenAA(30, 800, all_matches=True)
)


canonical_ikb_rule = Rules(
    Leader().localize_at("Cytoplasm")

    & Pfam.matches("PF00023").times(3, 11)
    & Sequence.length_between(300, 900)
    & Pfam.matches_only("PF00023", "PF12796", "PF13606", "PF13637", "PF13857")
   
    # Canonical IκBs require an intact N-terminal Signal-Receiving Domain (SRD)
    # ahead of Ankyrin Repeat 1. Minimum Leader Length: The SRD needs at least
    # 30 to 70 amino acids to house the target IKK phosphorylation serines and
    # lysine acceptor residues for ubiquitin attachment.
    #
    # Also want to limit the last match to be between 80 AAs from end of
    # protein, to allow for PEST.
    #
    & Pfam.matches("PF00023").betweenAA(30, 800, all_matches=True)

    # required IKK phosphorylation and then degradation by βTrCP
    & Sequence.matches_regex(βTrCP_pattern).relativeToPfam("PF00023", -200, -10)
)

canonical_bcl3_rule = Rules(
    Leader().localize_at("Nucleus")

    & Pfam.matches("PF00023").times(3, 11)
    & Sequence.length_between(300, 1000)
    & Pfam.matches_only("PF00023", "PF12796", "PF13606", "PF13637", "PF13857")

    # BCL-3 features an N-terminal leader of ~118 amino acids preceding Ankyrin
    # Repeat 1 (which begins at residue ~119 in human BCL-3). This region
    # contains its basic Nuclear Localization Signal (NLS) and
    # proline/serine-rich regulation sites. Conservating setting start of ANK
    # at 80.
    #
    # Also want to limit the last match to be between 80 to 200 AAs from end of
    # protein, to allow for a TAD.
    #
    & Pfam.matches("PF00023").betweenAA(80, 800, all_matches=True)

    # Does not have a DSG..S site
    & ~Sequence.matches_regex(βTrCP_pattern).relativeToPfam("PF00023", -200, -10)
)
