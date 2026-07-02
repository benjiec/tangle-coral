from sieve.rules import HMMAlignment, KO, Leader, Pfam, Rules, TFMotifs

rule = Rules(
    Leader.is_mTP()
    & Pfam.matches("PF00081")
    & Pfam.matches("PF02777")
    & KO.matches("K04564")
    & HMMAlignment("k04564.hmm").is_at("H", 61)
    & HMMAlignment("k04564.hmm").is_at("H", 106)
    & HMMAlignment("k04564.hmm").is_at("D", 187)
    & HMMAlignment("k04564.hmm").is_at("H", 191)
    & (
        HMMAlignment("k04564.hmm").is_at("DVWEHAYYLQ", 187)
        | HMMAlignment("k04564.hmm").is_at("DVWEHAYYVQ", 187)
    )
    & TFMotifs.has_within(
        20,
        "GM.5.0.Rel",
        "GM.5.0.bZIP",
        min_score_threshold=8,
    ).in_intron(2)
)
