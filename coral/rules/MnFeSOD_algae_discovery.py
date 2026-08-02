from sieve.rules import HMMAlignment, KO, Leader, Pfam, Rules, TFMotifs
from sieve.result_filters import Field, FieldRegex, LeaderCall

rule_fasta = Rules(
    Leader().upstreamOfPfam("PF00081").betweenAA(-70,-15).is_SP(deeploc=True)

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
