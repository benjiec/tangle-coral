import os
import csv
import argparse
from needle.duckdb import load, Clusters
from needle.seq import read_fasta_as_dict

parser = argparse.ArgumentParser()
parser.add_argument("module")
parser.add_argument("ko")
parser.add_argument("output")
args = parser.parse_args()

load(args.module.lower())
clusters = Clusters(args.ko.upper()).clusters()
print(len(clusters), "clusters")

tall_table = []

for cluster in clusters:
    print(cluster.cluster_fn())
    seq_dict = read_fasta_as_dict(f"data/{args.module.lower()}_results/clusters/{cluster.cluster_fn()}.faa")
    proteins = cluster.proteins(seq_dict)
    print(len(proteins), "proteins")

    aln_fn = f"data/{args.module.lower()}_results/alignments/{cluster.cluster_fn()}.faa"
    if not os.path.exists(aln_fn):
        continue
    alignments = read_fasta_as_dict(aln_fn)

    for protein in proteins:
        if protein.protein_accession not in alignments:
            continue
        aligned = alignments[protein.protein_accession].upper()

        protein_sequence = protein.sequence
        aligned_no_gap = aligned.replace("-", "")
        assert aligned_no_gap in protein_sequence
        front_ntrimmed = protein_sequence.index(aligned_no_gap)

        pro_pos = front_ntrimmed
        for aln_pos, aln in enumerate(aligned):
            # don't ignore X here - goal is to report protein positions that were aligned
            if aln != "-":
                pro_pos += 1

                ko_match = protein.ko_match_at(pro_pos)
                pfam_match = protein.pfam_match_at(pro_pos)

                data = dict(
                    cluster_id=cluster.cluster_id,
                    parent_cluster_id=cluster.parent_cluster_id,
                    protein_accession=protein.protein_accession,
                    genome_accession=protein.genome_accession,
                    alignment_pos=aln_pos+1,
                    protein_pos=pro_pos,
                    aa = aln,
                    ko_accession = None if ko_match is None else ko_match['ko_accession'],
                    ko_match = None if ko_match is None else f"{ko_match['ko_hmm_start']}-{ko_match['ko_hmm_end']}",
                    ko_evalue = None if ko_match is None else ko_match['ko_evalue'],
                    ko_score = None if ko_match is None else ko_match['ko_score'],
                    pfam_accession = None if pfam_match is None else pfam_match['pfam_accession'],
                    pfam_match = None if pfam_match is None else f"{pfam_match['pfam_hmm_start']}-{pfam_match['pfam_hmm_end']}",
                    pfam_evalue = None if pfam_match is None else pfam_match['pfam_evalue'],
                )

                tall_table.append(data) 

if len(tall_table):
    header = tall_table[0].keys()
    with open(args.output, "w") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter='\t')
        writer.writeheader()
        for row in tall_table:
            writer.writerow(row)
