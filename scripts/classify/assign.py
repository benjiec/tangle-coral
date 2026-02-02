import os
import argparse
from pathlib import Path
from scripts.defaults import DefaultPath
from needle.duckdb import load, CandidateClassifiedProteins, write_tsv_from_records, module_ko_ids
from needle.seq import read_fasta_as_dict, write_fasta_from_dict

ap = argparse.ArgumentParser()
ap.add_argument("module_id")
args = ap.parse_args()

load(args.module_id, load_final_assets=False)
module_kos = module_ko_ids(args.module_id)
candidate_proteins = CandidateClassifiedProteins()

#
# first, create manifests and fragments
#

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

#
# reload data, including newly created manifests, before generating assignments, which uses manifests
#

load(args.module_id)

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
faa_ko_assigned = {}  # which KO:proteins to save to FAA files
faa_ko_putative = {}  # which KO:proteins to save to FAA files

for ko_match in ko_matches:
    ko_id = ko_match['hmm_accession']
    k = (ko_match['protein_accession'], ko_match['genome_accession'])
    if k in ko_assignments and ko_assignments[k] == ko_id:
        ko_match['assignment_status'] = 'assigned'
        if ko_id in module_kos:
            faa_ko_assigned.setdefault(ko_id, []).append(k)
    elif k in ko_assignments:
        ko_match['assignment_status'] = 'other'
    else:
        ko_match['assignment_status'] = 'putative'
        if ko_id in module_kos:
            faa_ko_putative.setdefault(ko_id, []).append(k)

output_ko = f"data/{args.module_id}_results/candidate_ko.tsv"
output_pf = f"data/{args.module_id}_results/candidate_pfam.tsv"
write_tsv_from_records(output_ko, ko_matches)
write_tsv_from_records(output_pf, pfam_matches)


#
# write FAA files
#

protein_faa = f"data/{args.module_id}_results/protein_detected.faa"
output_dir = f"data/{args.module_id}_results/faa"
os.makedirs(output_dir, exist_ok=True)

# write two FAAs files for each KO in the module, assigned and putative
genome_with_proteins = {}
protein_sequences = read_fasta_as_dict(protein_faa)

ko_assigned_protein_sequences = {}
for ko_id, items in faa_ko_assigned.items():
    for protein_accession, genome_accession in items:
        if protein_accession in protein_sequences:
            protein_sequence = protein_sequences[protein_accession]
        else:
            if genome_accession not in genome_with_proteins:
                protein_fn = DefaultPath.ncbi_genome_protein_faa(genome_accession)
                print("reading sequences from", protein_fn)
                protein_sequences.update(read_fasta_as_dict(protein_fn))
            protein_sequence = protein_sequences[protein_accession]
        ko_assigned_protein_sequences.setdefault(ko_id, {})[protein_accession] = protein_sequence

ko_putative_protein_sequences = {}
for ko_id, items in faa_ko_putative.items():
    for protein_accession, genome_accession in items:
        if protein_accession in protein_sequences:
            protein_sequence = protein_sequences[protein_accession]
        else:
            if genome_accession not in genome_with_proteins:
                protein_fn = DefaultPath.ncbi_genome_protein_faa(genome_accession)
                print("reading sequences from", protein_fn)
                protein_sequences.update(read_fasta_as_dict(protein_fn))
            protein_sequence = protein_sequences[protein_accession]
        ko_putative_protein_sequences.setdefault(ko_id, {})[protein_accession] = protein_sequence

for ko_id, sequence_dict in ko_assigned_protein_sequences.items():
    faa_fn = f"data/{args.module_id}_results/faa/{ko_id}.faa"
    write_fasta_from_dict(sequence_dict, faa_fn)
    print(faa_fn)

for ko_id, sequence_dict in ko_putative_protein_sequences.items():
    faa_fn = f"data/{args.module_id}_results/faa/{ko_id}-putative.faa"
    write_fasta_from_dict(sequence_dict, faa_fn)
    print(faa_fn)
