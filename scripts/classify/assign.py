import os
import argparse
from pathlib import Path
from scripts.defaults import DefaultPath
from needle.duckdb import load, CandidateClassifiedProteins, write_tsv_from_records

ap = argparse.ArgumentParser()
ap.add_argument("module_id")
ap.add_argument("--additional-genome-accession", type=str, default=None)
args = ap.parse_args()

protein_faa = f"data/{args.module_id}_results/proteins.faa"
output_dir = f"data/{args.module_id}_results/faa"
os.makedirs(output_dir, exist_ok=True)

multiple_protein_faas = []
if args.additional_genome_accession:
    if os.path.exists(args.additional_genome_accession):
        with open(args.additional_genome_accession, "r") as f:
            for acc in f.readlines():
                acc = acc.strip()
                multiple_protein_faas.append(DefaultPath.ncbi_genome_protein_faa(acc))
                print("additional protein sequences from", multiple_protein_faas[-1])
    else:
        multiple_protein_faas.append(DefaultPath.ncbi_genome_protein_faa(args.additional_genome_accession))
        print("additional protein sequences from", multiple_protein_faas[-1])

load(args.module_id)
candidate_proteins = CandidateClassifiedProteins()

# matches
ko_matches = candidate_proteins.ko_matches()
ko_matches = ko_matches.to_dict(orient='records')
pfam_matches = candidate_proteins.pfam_matches()
pfam_matches = pfam_matches.to_dict(orient='records')

match_sorter = lambda d: (d['genome_accession'], d['protein_accession'], d['hmm_db'], d['dom_rank_for_protein'], d['hmm_accession'], d['hmm_start'], d['hmm_end'])
ko_matches = sorted(ko_matches, key=match_sorter)
pfam_matches = sorted(pfam_matches, key=match_sorter)

ko_assignments = candidate_proteins.ko_assignments(0.9, 1.0)
ko_assignments = ko_assignments.to_dict(orient='records')
ko_assignments = {(d['protein_accession'], d['genome_accession']):d['hmm_accession'] for d in ko_assignments}
for ko_match in ko_matches:
    k = (ko_match['protein_accession'], ko_match['genome_accession'])
    if k in ko_assignments and ko_assignments[k] == ko_match['hmm_accession']:
        ko_match['assignment_status'] = 'assigned'
    else:
        ko_match['assignment_status'] = 'putative'

output_ko = f"data/{args.module_id}_results/candidate_ko.tsv"
output_pf = f"data/{args.module_id}_results/candidate_pfam.tsv"
write_tsv_from_records(output_ko, ko_matches)
write_tsv_from_records(output_pf, pfam_matches)

# protein manifest
proteins = candidate_proteins.proteins()
proteins = proteins.to_dict(orient='records')

manifest_sorter = lambda d: (d['genome_accession'], d['major_contig'], d['proteome_type'], d['protein_accession'])
proteins = sorted(proteins, key=manifest_sorter)

output_manifest = f"data/{args.module_id}_results/proteins.tsv"
write_tsv_from_records(output_manifest, proteins)

# fragments
ncbi_fragments = candidate_proteins.ncbi_fragments()
ncbi_fragments = ncbi_fragments.to_dict(orient='records')
detected_fragments = candidate_proteins.detected_fragments()
detected_fragments = detected_fragments.to_dict(orient='records')
fragments = ncbi_fragments+detected_fragments

fragments_sorter = lambda d: (d['genome_accession'], d['target_accession'], d['protein_hit_id'], d['query_accession'], d['target_start'])
fragments = sorted(fragments, key=fragments_sorter)

output_fragments = f"data/{args.module_id}_results/protein_fragments.tsv"
write_tsv_from_records(output_fragments, fragments)



"""
#
# first, non-reference genome, use default scoring threshold
#
assign_ko(classify_rows, ko_hmm_db_name, protein_faa, output_dir, score_to_threshold_ratio=0.9)

#
# then, reference genome, use more stringent scoring threshold
#

# use more stringent scoring threshold
assign_ko(classify_rows, ko_hmm_db_name, multiple_protein_faas, output_dir, score_to_threshold_ratio=1)
"""
