import needle.duckdb
from scripts.defaults import DefaultPath
from needle.seq import read_fasta_as_dict
from needle.detect import hmm_search_genome
from needle.hmm import HMMCollection

        
def search_gaps(hmm_file, hmm_name, genome_accession, target_accession, strand, gap_coords):

    hmm_c = HMMCollection(hmm_file, [hmm_name])
    hmm_file = hmm_c.get(hmm_name)

    fna_file = DefaultPath.ncbi_genome_fna(genome_accession)
    genomic_fasta = read_fasta_as_dict(fna_file)

    for query_start, query_end, target_start, target_end in gap_coords:
        if target_start is None and target_end is None:
            continue
        if target_start is None:
            if strand > 0:
                target_start = max(target_end - 20000, 1)
            else:
                target_start = target_end + 20000
        if target_end is None:
            if strand > 0:
                target_end = target_start + 20000
            else:
                target_end = max(target_start - 20000, 1)

        print("searching", target_start, target_end)
        hmm_rows = hmm_search_genome(
            hmm_file, genome_accession, genomic_fasta,
            target_accession = target_accession,
            target_left = min(target_start, target_end),
            target_right = max(target_start, target_end),
            strand = strand,
            conditional = True
        )
        for row in hmm_rows:
            print(row)

    hmm_c.clean()


def find_gaps(min_query_start, max_query_end, strand, existing_query_target_coords, query_gap_threshold=8):

    existing_query_target_coords = sorted(existing_query_target_coords)
    bookended_coords = [(None, min_query_start-1, None, None)] + existing_query_target_coords + [(max_query_end+1, None, None, None)]
    gap_coords = []

    for i in range(len(bookended_coords)-1):
        this = bookended_coords[i]
        next = bookended_coords[i+1]

        print("    ", this)
        if next[0]-this[1] > query_gap_threshold:  # there is a gap
            gap_coord = (this[1]+1, next[0]-1, this[3]+strand if this[3] is not None else None, next[2]-strand if next[2] is not None else None)
            gap_coords.append(gap_coord)
            print("        ", gap_coord)

    print("    ", next)

    return gap_coords


def process_candidates(df, hmm_file, hmm_name):

    df['query_match_len'] = df['query_end'] - df['query_start'] + 1
    min_query_start = df['query_start'].min()
    max_query_end = df['query_end'].max()

    for (protein_accession, genome_accession), group_df in df.groupby(['protein_accession', 'genome_accession']):
        if protein_accession not in ("K00164_CAXAMN010000001.1_9623db44", "K00164_CAXAMN010021374.1_607cb229", "K00164_CAXAMN010016668.1_1e5d987e"):
            continue
        print(protein_accession)

        target_accession = group_df['target_accession'].unique()
        assert len(target_accession) == 1
        target_accession = target_accession[0]

        target_start_first = group_df.iloc[0]['target_start']
        target_end_first = group_df.iloc[0]['target_end']
        if target_start_first <= target_end_first:
            strand = 1
        else:
            strand = -1

        query_target_coords = list(zip(group_df['query_start'], group_df['query_end'], group_df['target_start'], group_df['target_end']))

        print(protein_accession, genome_accession, target_accession, target_start_first, target_end_first, strand)
        gap_coords = find_gaps(min_query_start, max_query_end, strand, query_target_coords)
        search_gaps(hmm_file, hmm_name, genome_accession, target_accession, strand, gap_coords)


needle.duckdb.load("m00009")
needle.duckdb.assignments()
df = needle.duckdb.candidates("K00164")
process_candidates(df, "data/m00009_ko.hmm", "K00164")
