import time
import argparse
from Bio.Seq import Seq
from defaults import DefaultPath
from needle.match import read_fasta_as_dict, extract_subsequence
from needle.hits import compute_three_frame_translations
from needle.hmm import hmmsearch_sequence_dict

ap = argparse.ArgumentParser()
ap.add_argument("hmm_file")
ap.add_argument("genome_accession")
ap.add_argument("--target-accession", type=str, default=None)
ap.add_argument("--target-left", type=int, default=None)
ap.add_argument("--target-right", type=int, default=None)
args = ap.parse_args()

win = 50000
win_overlap = 10000
genome_accession = args.genome_accession
fna_file = DefaultPath.ncbi_genome_fna(genome_accession)
genomic_fasta = read_fasta_as_dict(fna_file)


def to_dna_coordinate(frame_dna_start, frame_dna_end, ali_from, ali_to):
    if frame_dna_end > frame_dna_start: # fwd strand
        dna_ali_from = frame_dna_start+(ali_from-1)*3
        dna_ali_to = frame_dna_start+ali_to*3-1
    else: # rev strand 
        dna_ali_from = frame_dna_start-(ali_from-1)*3
        dna_ali_to = frame_dna_start-ali_to*3+1
    return dna_ali_from, dna_ali_to

total_time_translation = 0
total_time_hmmsearch = 0

for acc, genome_sequence in genomic_fasta.items():
    if args.target_accession and acc != args.target_accession:
        continue

    contig_start = 0
    contig_end = len(genome_sequence)
    if args.target_left:
        contig_start = args.target_left - 1
    if args.target_right:
        contig_end = args.target_right

    for win_i in range(contig_start, contig_end, win):
        translated_fasta = {}
        subs = genome_sequence[win_i:win_i+win+win_overlap]
        t0 = time.time()
        translations_fwd = compute_three_frame_translations(subs, 1, len(subs))
        translations_rev = compute_three_frame_translations(subs, len(subs), 1)
        total_time_translation += time.time()-t0
        translated_fasta[f"{acc}_fwd_0"] = translations_fwd[0][2]
        translated_fasta[f"{acc}_fwd_1"] = translations_fwd[1][2]
        translated_fasta[f"{acc}_fwd_2"] = translations_fwd[2][2]
        translated_fasta[f"{acc}_rev_0"] = translations_rev[0][2]
        translated_fasta[f"{acc}_rev_1"] = translations_rev[1][2]
        translated_fasta[f"{acc}_rev_2"] = translations_rev[2][2]
        for k,v in translated_fasta.items():
            # print(f">{k}\n{v}")
            pass

        t0 = time.time()
        hmm_rows = hmmsearch_sequence_dict(args.hmm_file, translated_fasta)
        total_time_hmmsearch += time.time()-t0

        for row in hmm_rows:
            frame = int(row["target_name"][-1])
            if "_fwd_" in row["target_name"]:
                frame_dna_start = translations_fwd[frame][0]
                frame_dna_end = translations_fwd[frame][1]
                aa_seq = extract_subsequence(translations_fwd[frame][2], row["ali_from"], row["ali_to"])
            else:
                frame_dna_start = translations_rev[frame][0]
                frame_dna_end = translations_rev[frame][1]
                aa_seq = extract_subsequence(translations_rev[frame][2], row["ali_from"], row["ali_to"])

            assert (abs(frame_dna_end-frame_dna_start)+1)%3 == 0

            dna_ali_from, dna_ali_to = to_dna_coordinate(frame_dna_start+win_i, frame_dna_end+win_i, row["ali_from"], row["ali_to"])
            contig_accession = row["target_name"][:-6]

            out = [
              row["query_accession"],
              contig_accession,
              str(row["evalue"]),
              "",
              str(row["hmm_from"]),
              str(row["hmm_to"]),
              str(dna_ali_from),
              str(dna_ali_to),
              aa_seq
            ]
            print("\t".join(out))

print(total_time_translation)
print(total_time_hmmsearch)
