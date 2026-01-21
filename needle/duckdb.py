import duckdb
import pandas as pd


def load(module_id):
    duckdb.execute("CREATE SCHEMA needle")

    duckdb.execute("CREATE TABLE needle.pfams AS SELECT * FROM 'data/Pfam-A.clans.tsv'")
    duckdb.execute("CREATE TABLE needle.module_steps AS SELECT * FROM 'data/module_defs.csv'")
    duckdb.execute("CREATE TABLE needle.modules AS SELECT * FROM 'data/modules.tsv'")
    duckdb.execute("CREATE TABLE needle.ko AS SELECT * FROM read_csv_auto('data/ko.tsv', normalize_names=TRUE)")
    duckdb.execute("CREATE TABLE needle.genomes AS SELECT * FROM read_csv_auto('data/genomes.tsv', normalize_names=TRUE)")

    duckdb.execute(f"CREATE TABLE needle.ko_match AS SELECT * FROM 'data/{module_id}_results/classify.tsv' WHERE hmm_db = '{module_id}_ko'")
    duckdb.execute(f"CREATE TABLE needle.pfam_match AS SELECT * FROM 'data/{module_id}_results/classify.tsv' WHERE hmm_db = 'Pfam-A'")
    duckdb.execute(f"CREATE TABLE needle.protein_names AS SELECT * FROM 'data/{module_id}_results/protein_names.tsv'")
    duckdb.execute(f"CREATE TABLE needle.detected AS SELECT * FROM 'data/{module_id}_results/proteins.tsv'")
    duckdb.execute(f"CREATE TABLE needle.ncbi_exons AS SELECT * FROM 'data/{module_id}_results/protein_ncbi.tsv'")


def assignments():
    con = duckdb.connect(":default:")
    df = con.sql("""
        SELECT ko_match.protein_accession, ko_match.genome_accession, ko_match.hmm_accession, ko_match.dom_score, ko_match.score_threshold
          FROM needle.ko_match
         WHERE ko_match.dom_rank_for_protein == 1
           AND ((ko_match.protein_accession LIKE 'K%' AND ko_match.dom_score / ko_match.score_threshold >= 0.9) OR
                (ko_match.dom_score / ko_match.score_threshold >= 1.0))
  """).df()

    res = df.to_dict(orient='records')
    print(len(res))
    print(df['protein_accession'].nunique())

    df_all_duplicates = df[df.duplicated(subset=['protein_accession'], keep=False)]
    print(df_all_duplicates)


def candidates(ko_number, evalue_threshold=1e-20):

    con = duckdb.connect(":default:")
    df = con.sql("""
        SELECT ko_match.protein_accession, ko_match.genome_accession, ko_match.hmm_accession,
               detected.target_accession, detected.target_start, detected.target_end,
               detected.query_accession, detected.query_start, detected.query_end
          FROM needle.ko_match
          JOIN needle.detected ON ko_match.protein_accession = detected.protein_hit_id AND ko_match.hmm_accession = detected.query_accession
         WHERE ko_match.dom_rank_for_protein == 1
           AND ko_match.hmm_accession = '%s'
           AND ko_match.dom_evalue < %s
  """ % (ko_number, evalue_threshold)).df()

    print(df)
    return df


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


def group_candidates_by_protein(df):

    df['query_match_len'] = df['query_end'] - df['query_start'] + 1
    min_query_start = df['query_start'].min()
    max_query_end = df['query_end'].max()

    for (protein_accession, genome_accession), group_df in df.groupby(['protein_accession', 'genome_accession']):

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
