from sieve.rules import HMMAlignment, KO, Leader, Pfam, Rules, TFMotifs
from sieve.result_filters import Field, FieldRegex, LeaderCall

rule_tangle_curated = Rules(
    Leader().upstreamOfPfam("PF00081").betweenAA(-70,-15)

    # requires SP+cTP architecture
    & HMMAlignment("algae_sod2_sp_ctp_leader.hmm").spans(3, 50).between(1, 60)

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

rule_fasta = rule_tangle_curated

is_positive = (
      FieldRegex(r"HMMAlignment.+").all().eq("true")
    & FieldRegex(r"Pfam.matches.+").all().eq("true")
)
