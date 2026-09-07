from sieve.rules import HMMAlignment, Leader, Pfam, Rules, Sequence

# lysine (K) within 40 AAs upstream of a DSG(\phi).S degron motif
βTrCP_pattern = r"K.{1,40}?[DEST][ST]G[LIVMFYWA].{1,2}[ST]"

ANK_PF = "PF00023"
All_ANK_PFs = ("PF00023", "PF12796", "PF13606", "PF13637", "PF13857")
min_protein_length = 300
IkB_hypothetical_length = 900
Bcl3_hypothetical_length = 1000


# all encompassing rule to filter sequences so we can send them to DeepLoc
no_deeploc_rule = Rules(

    Pfam.matches(ANK_PF).times(4, 9)
    & Sequence.length_between(min_protein_length, max(IkB_hypothetical_length, Bcl3_hypothetical_length))
    & Pfam.matches_only(*All_ANK_PFs)

    & Pfam.matches(ANK_PF).betweenAA(60, max(IkB_hypothetical_length, Bcl3_hypothetical_length)-100, all_matches=True)
)


canonical_ikb_rule = Rules(
    Leader().localize_at("Cytoplasm")

    & Pfam.matches(ANK_PF).times(4, 9)
    & Sequence.length_between(min_protein_length, IkB_hypothetical_length)
    & Pfam.matches_only(*All_ANK_PFs)
   
    # Canonical IκBs require an intact N-terminal Signal-Receiving Domain (SRD)
    # ahead of Ankyrin Repeat 1. Minimum Leader Length: NES+SRD. The SRD needs
    # at least 30 to 70 amino acids to house the target IKK phosphorylation
    # serines and lysine acceptor residues for ubiquitin attachment. Also
    # leaving room at the end for PEST sequence.

    & Pfam.matches(ANK_PF).betweenAA(60, IkB_hypothetical_length-50, all_matches=True)

    # required IKK phosphorylation and then degradation by βTrCP
    & Sequence.matches_regex(βTrCP_pattern).relativeToPfam("PF00023", -200, -10)
)

canonical_bcl3_rule = Rules(
    Leader().localize_at("Nucleus")

    & Pfam.matches(ANK_PF).times(4, 9)
    & Sequence.length_between(min_protein_length, Bcl3_hypothetical_length)
    & Pfam.matches_only(*All_ANK_PFs)

    # BCL-3 features an N-terminal leader of ~118 amino acids preceding Ankyrin
    # Repeat 1 (which begins at residue ~119 in human BCL-3). This region
    # contains its basic Nuclear Localization Signal (NLS) and
    # proline/serine-rich regulation sites. Conservating setting start of ANK
    # at 80. Also leaving room at the end to allow for a TAD.

    & Pfam.matches(ANK_PF).betweenAA(80, Bcl3_hypothetical_length-150, all_matches=True)

    # Does not have a DSG..S site
    & ~Sequence.matches_regex(βTrCP_pattern).relativeToPfam("PF00023", -200, -10)
)
