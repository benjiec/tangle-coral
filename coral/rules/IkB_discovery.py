from sieve.rules import HMMAlignment, Leader, Pfam, Rules, Sequence

# lysine (K) within 40 AAs upstream of a DSG(\phi).S degron motif
βTrCP_pattern = r"K.{1,40}?[DEST][ST]G[LIVMFYWA].{1,2}[ST]"

ANK_PF = "PF00023"
All_ANK_PFs = ("PF00023", "PF12796", "PF13606", "PF13637", "PF13857")
min_protein_length = 300
IkB_hypothetical_length = 900
Bcl3_hypothetical_length = 1000

permissive_ard_rule = Rules(
    Sequence.length_between(min_protein_length, 2000)
    & Pfam.matches(ANK_PF).times(4,9).betweenAA(60, 1800, all_matches=True).spansAA(130,280)
)

ikb_like_leader_rule = Rules(
    # export bound RHD-containing dimers back out
    Leader().is_NES()
)

kb_rule = Rules(
    # regulatory mechanism
    & TFMotifs.has(
        "GM.5.0.Rel",
    ).betweenBED(-1500, 500)
)

ikb_like_rule = Rules(
    # export bound RHD-containing dimers back out
    Leader().is_NES()

    # regulatory mechanism
    & TFMotifs.has(
        "GM.5.0.Rel",
    ).betweenBED(-1500, 500)
)

canonical_ikb_rule = Rules(
    Leader().is_NES()
    & Pfam.matches(ANK_PF).times(4, 9)
    & Sequence.length_between(min_protein_length, IkB_hypothetical_length)
    & Pfam.matches_only(*All_ANK_PFs)
   
    # Canonical IκBs require an intact N-terminal Signal-Receiving Domain (SRD)
    # ahead of Ankyrin Repeat 1. Minimum Leader Length: NES+SRD. The SRD needs
    # at least 30 to 70 amino acids to house the target IKK phosphorylation
    # serines and lysine acceptor residues for ubiquitin attachment. Also
    # leaving room at the end for PEST sequence, up to 50 AAs.

    & Pfam.matches(ANK_PF).betweenAA(60, IkB_hypothetical_length-50, all_matches=True)

    # required IKK phosphorylation and then degradation by βTrCP
    & Sequence.matches_regex(βTrCP_pattern).relativeToPfam("PF00023", -200, -10)
)

canonical_bcl3_rule = Rules(
    Leader().is_NLS()
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
