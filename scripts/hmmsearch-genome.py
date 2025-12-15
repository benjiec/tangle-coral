import argparse
from defaults import DefaultPath
from needle.match import read_fasta_as_dict
from needle.detect import hmm_search_genome_sequence, Results

ap = argparse.ArgumentParser()
ap.add_argument("hmm_file")
ap.add_argument("genome_accession")
ap.add_argument("output_file")
ap.add_argument("--target-accession", type=str, default=None)
ap.add_argument("--target-left", type=int, default=None)
ap.add_argument("--target-right", type=int, default=None)
args = ap.parse_args()

win = 50000
win_overlap = 10000
genome_accession = args.genome_accession
fna_file = DefaultPath.ncbi_genome_fna(genome_accession)
genomic_fasta = read_fasta_as_dict(fna_file)

detected = []

for acc, genome_sequence in genomic_fasta.items():
    if args.target_accession and acc != args.target_accession:
        continue

    if args.target_left:
        contig_left_0b = args.target_left - 1
    else:
        contig_left_0b = None
    if args.target_right:
        contig_right_excl_0b = args.target_right
    else:
        contig_right_excl_0b = None

    detected.extend(hmm_search_genome_sequence(
        args.hmm_file, acc, genome_sequence, win, win_overlap,
        contig_left_0b=contig_left_0b,
        contig_right_excl_0b=contig_right_excl_0b
    ))

with open(args.output_file, "w") as f:
    f.write("\t".join(Results.PRODUCER_HEADER)+"\n")
    for d in detected:
        assert len(d) == len(Results.PRODUCER_HEADER)
        f.write("\t".join([str(x) for x in d])+"\n")
