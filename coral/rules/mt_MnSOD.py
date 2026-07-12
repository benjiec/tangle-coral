from sieve.rules import HMMAlignment, KO, Leader, Pfam, Rules, TFMotifs
from sieve.result_filters import Field, FieldRegex, LeaderCall

rule_tangle_curated = Rules(
    # required for mitochondrial function rather than secreted
    Leader().upstreamOfPfam("PF00081").betweenAA(-45,0).is_mTP()

    # required domains
    & Pfam.matches("PF00081")
    & Pfam.matches("PF02777")
    & KO.matches("K04564")

    # required catalytic residues to form metal binding pocket
    & HMMAlignment("k04564.hmm").is_at("H", 61)
    & HMMAlignment("k04564.hmm").is_at("H", 106)
    & HMMAlignment("k04564.hmm").is_at("D", 187)
    & HMMAlignment("k04564.hmm").is_at("H", 191)

    # proton gate-keeper
    & HMMAlignment("k04564.hmm").matches_regex("[QDK]", 196)

    # regulatory mechanism
    & TFMotifs.has_within(
        20,
        "GM.5.0.Rel",
        "GM.5.0.bZIP",
        min_score_threshold=8,
    )
)

rule_fasta = Rules(
    # required for mitochondrial function rather than secreted
    Leader().upstreamOfPfam("PF00081").betweenAA(-45,0).is_mTP()

    # required domains
    & Pfam.matches("PF00081")
    & Pfam.matches("PF02777")
    & KO.matches("K04564")

    # required catalytic residues to form metal binding pocket
    & HMMAlignment("k04564.hmm").is_at("H", 61)
    & HMMAlignment("k04564.hmm").is_at("H", 106)
    & HMMAlignment("k04564.hmm").is_at("D", 187)
    & HMMAlignment("k04564.hmm").is_at("H", 191)

    # proton gate-keeper
    & HMMAlignment("k04564.hmm").matches_regex("[QDK]", 196)
)

is_positive = (
      FieldRegex(r"HMMAlignment.+").all().eq("true")
    & FieldRegex(r"Pfam.matches.+").all().eq("true")
    & LeaderCall("mTP").ge(10)
)
