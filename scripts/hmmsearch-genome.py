import argparse
from defaults import DefaultPath
from needle.match import read_fasta_as_dict
from needle.detect import hmm_search_genome, Results

ap = argparse.ArgumentParser()
ap.add_argument("hmm_file")
ap.add_argument("genome_accession")
ap.add_argument("output_file")
ap.add_argument("--target-accession", type=str, default=None)
ap.add_argument("--target-left", type=int, default=None)
ap.add_argument("--target-right", type=int, default=None)
args = ap.parse_args()

genome_accession = args.genome_accession
fna_file = DefaultPath.ncbi_genome_fna(genome_accession)
genomic_fasta = read_fasta_as_dict(fna_file)

hmm_rows = hmm_search_genome(
    args.hmm_file, genome_accession, genomic_fasta,
    target_accession = args.target_accession,
    target_left = args.target_left,
    target_right = args.target_right,
    conditional = True  # setting this argument to true, so we can use c-Evalue, which is more permissive
)

detected = []
for row in hmm_rows:
    out = (
      row["query_accession"],
      row["target_accession"],
      row["dom_evalue"],
      '',
      row["hmm_from"],
      row["hmm_to"],
      row["ali_from"],
      row["ali_to"],
      row["matched_sequence"]
    )
    detected.append(out)

with open(args.output_file, "w") as f:
    f.write("\t".join(Results.PRODUCER_HEADER)+"\n")
    for d in detected:
        assert len(d) == len(Results.PRODUCER_HEADER)
        f.write("\t".join([str(x) for x in d])+"\n")
