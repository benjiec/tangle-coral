from sieve.rules import Rules, Sequence

nothing = Rules(
    Sequence.length_between(1, 1000000)
)
