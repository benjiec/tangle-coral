#!/usr/bin/env python3
import argparse

from needle.detect import Results
from needle.match import group_matches
from needle.hits import hmm_find_proteins, hmm_clean
from needle.io import export_protein_hits
from needle.hmm import HMMCollection
from defaults import DefaultPath

def main():
    parser = argparse.ArgumentParser(description="Export protein matches from detection results TSV.")
    parser.add_argument("results_tsv", help="Detection results TSV (Results.PRODUCER headers)")
    parser.add_argument("genome_accession", help="Genome accession")
    parser.add_argument("hmm_file", help="HMM file to use for improve search results")
    parser.add_argument("output_dir", help="Path of output directory")
    args = parser.parse_args()

    target_fasta = DefaultPath.ncbi_genome_fna(args.genome_accession)

    res = Results(args.results_tsv, target_fasta_path=target_fasta)
    protein_matches = group_matches(res.matches())
    accession_ids = [pm.query_accession for pm in protein_matches]
    hmm_collection = HMMCollection(args.hmm_file, accession_ids)

    try:

        """
        # If we are using BLAST as the main "detection" method, the following
        # may find more matches, but we are now using hmmsearch on all AA
        # fragments, so the following is likely not adding much

        # use HMM to find more fragments
        protein_matches = hmm_find_proteins(protein_matches, res, hmm_collection)
        """

        pre_filter = len(protein_matches)
        protein_matches = [m for m in protein_matches if m.can_collate()]
        print("filter by collatable", pre_filter, "=>", len(protein_matches))
        cleaned_protein_matches = hmm_clean(protein_matches, hmm_collection)

        export_protein_hits(
            args.genome_accession,
            cleaned_protein_matches,
            args.output_dir+"/proteins.tsv",
            args.output_dir+"/matches.tsv",
            args.output_dir+"/faa"
        )

    finally:
        hmm_collection.clean()


if __name__ == "__main__":
    main()
