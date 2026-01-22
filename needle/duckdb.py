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
