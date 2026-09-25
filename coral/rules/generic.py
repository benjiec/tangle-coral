from sieve.rules import Rules, Sequence, TFMotifs

nothing = Rules(
    Sequence.length_between(1, 1000000)
)

loci = Rules(
    TFMotifs.has("something")
)
