import duckdb
import pandas as pd


def load(module_id):

    duckdb.execute("CREATE SCHEMA needle")

    duckdb.execute("CREATE TABLE needle.pfam_names AS SELECT * FROM 'data/Pfam-A.clans.tsv'")
    duckdb.execute("CREATE TABLE needle.module_defs AS SELECT * FROM 'data/module_defs.csv'")
    duckdb.execute("CREATE TABLE needle.module_names AS SELECT * FROM 'data/modules.tsv'")
    duckdb.execute("CREATE TABLE needle.ko AS SELECT * FROM read_csv_auto('data/ko.tsv', normalize_names=TRUE)")
    duckdb.execute("CREATE TABLE needle.genomes AS SELECT * FROM read_csv_auto('data/genomes.tsv', normalize_names=TRUE)")

    duckdb.execute(f"CREATE TABLE needle.classified AS SELECT * FROM 'data/{module_id}_results/classify.tsv' WHERE hmm_db = '{module_id}_ko'")
    duckdb.execute(f"CREATE TABLE needle.pfam_match AS SELECT * FROM 'data/{module_id}_results/classify.tsv' WHERE hmm_db = 'Pfam-A'")
    duckdb.execute(f"CREATE TABLE needle.protein_names AS SELECT * FROM 'data/{module_id}_results/protein_names.tsv'")
    duckdb.execute(f"CREATE TABLE needle.detected AS SELECT * FROM 'data/{module_id}_results/proteins.tsv'")
    duckdb.execute(f"CREATE TABLE needle.ncbi_exons AS SELECT * FROM 'data/{module_id}_results/protein_ncbi.tsv'")


def example():
    con = duckdb.connect(":default:")
    res_df = con.sql("""
        SELECT *
        FROM needle.classified
        JOIN needle.ko ON needle.classified.hmm_accession = needle.ko.ortholog_id
        JOIN needle.pfam_match ON needle.classified.protein_accession = needle.pfam_match.protein_accession
        LIMIT 10
  """).df()

    res = res_df.to_dict(orient='records')
    print(res)
