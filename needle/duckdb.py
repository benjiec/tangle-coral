import csv
import duckdb
import pandas as pd


def load(module_id):
    duckdb.execute("CREATE SCHEMA needle")

    # used in Tableau and scripts

    # generic
    duckdb.execute("CREATE TABLE needle.genomes AS SELECT * FROM read_csv_auto('data/genomes.tsv', normalize_names=TRUE)")
    duckdb.execute("CREATE TABLE needle.ko AS SELECT * FROM read_csv_auto('data/ko.tsv', normalize_names=TRUE)")
    duckdb.execute("CREATE TABLE needle.module_steps AS SELECT * FROM 'data/module_defs.csv'")
    duckdb.execute("CREATE TABLE needle.modules AS SELECT * FROM 'data/modules.tsv'")
    duckdb.execute("CREATE TABLE needle.pfams AS SELECT * FROM 'data/Pfam-A.clans.tsv'")

    # module specific
    duckdb.execute(f"CREATE TABLE needle.proteins AS SELECT * FROM 'data/{module_id}_results/proteins.tsv'")
    duckdb.execute(f"CREATE TABLE needle.ko_match AS SELECT * FROM 'data/{module_id}_results/candidate_ko.tsv'")
    duckdb.execute(f"CREATE TABLE needle.pfam_match AS SELECT * FROM 'data/{module_id}_results/candidate_pfam.tsv'")
    duckdb.execute(f"CREATE TABLE needle.clusters AS SELECT * FROM read_csv_auto('data/{module_id}_results/clusters.tsv', normalize_names=TRUE)")
    duckdb.execute(f"CREATE TABLE needle.protein_fragments AS SELECT * FROM 'data/{module_id}_results/protein_fragments.tsv'")

    # used to generate the above tables

    duckdb.execute(f"CREATE TABLE needle.classify AS SELECT * FROM 'data/{module_id}_results/classify.tsv'")
    duckdb.execute(f"CREATE TABLE needle.protein_names AS SELECT * FROM 'data/{module_id}_results/protein_names.tsv'")
    duckdb.execute(f"CREATE TABLE needle.detected AS SELECT * FROM 'data/{module_id}_results/protein_detected.tsv'")
    duckdb.execute(f"CREATE TABLE needle.ncbi_exons AS SELECT * FROM 'data/{module_id}_results/protein_ncbi.tsv'")


def write_tsv_from_records(fn, records):
    fieldnames = list(records[0].keys())
    with open(fn, "w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)


class CandidateClassifiedProteins(object):

    def _selection_sql(self, evalue_threshold):
        return """SELECT DISTINCT protein_accession, genome_accession, hmm_accession
                    FROM needle.classify
                   WHERE hmm_db = 'ko'
                   GROUP BY protein_accession, genome_accession, hmm_accession
                  HAVING MIN(dom_evalue) < %s""" % evalue_threshold

    def __init__(self, evalue_threshold=1E-80):
        self.con = duckdb.connect(":default:")
        self.selection_sql = self._selection_sql(evalue_threshold)

    def protein_genome_accessions(self):
        df = self.con.sql(self.selection_sql).df()
        unique_pairs = df[['protein_accession','genome_accession']].drop_duplicates()
        return list(unique_pairs.itertuples(index=False, name=None))

    def pfam_matches(self):
        sql = """
            SELECT DISTINCT classify.protein_accession,
                   classify.genome_accession,
                   hmm_db, classify.hmm_accession, hmm_start, hmm_end, protein_start, protein_end,
                   dom_evalue_cond, dom_evalue, dom_score, score_threshold, dom_rank_for_protein
              FROM needle.classify
              JOIN (%s) as filtered ON classify.protein_accession = filtered.protein_accession AND classify.genome_accession = filtered.genome_accession
             WHERE hmm_db = 'Pfam-A'""" % self.selection_sql

        df = self.con.sql(sql).df()
        return df

    def ko_matches(self, incl_evalue_threshold=1E-10):
        sql = """
            SELECT DISTINCT classify.protein_accession,
                   classify.genome_accession,
                   hmm_db, classify.hmm_accession, hmm_start, hmm_end, protein_start, protein_end,
                   dom_evalue_cond, dom_evalue, dom_score, score_threshold, dom_rank_for_protein
              FROM needle.classify
              JOIN (%s) as filtered ON classify.protein_accession = filtered.protein_accession AND classify.genome_accession = filtered.genome_accession
             WHERE hmm_db = 'ko'
               AND (classify.hmm_accession = filtered.hmm_accession OR classify.dom_evalue < %s)
          """ % (self.selection_sql, incl_evalue_threshold)

        df = self.con.sql(sql).df()
        return df

    def ko_assignments(self, score_to_threshold_ratio_detected, score_to_threshold_ratio_reference):
        sql = """
            SELECT DISTINCT classify.protein_accession, classify.genome_accession, classify.hmm_accession
              FROM needle.classify
              JOIN (%s) as filtered ON classify.protein_accession = filtered.protein_accession AND classify.genome_accession = filtered.genome_accession
             WHERE hmm_db = 'ko'
               AND dom_rank_for_protein = 1
               AND (score_threshold = '-' OR
                    (LEFT(classify.protein_accession, 6) == classify.hmm_accession AND CAST(dom_score AS FLOAT) / CAST(score_threshold AS FLOAT) >= %s) OR
                    (LEFT(classify.protein_accession, 6) != classify.hmm_accession AND CAST(dom_score AS FLOAT) / CAST(score_threshold AS FLOAT) >= %s))
          """ % (self.selection_sql, score_to_threshold_ratio_detected, score_to_threshold_ratio_reference)

        df = self.con.sql(sql).df()
        return df

    def proteins(self):
        sql = """
            SELECT DISTINCT filtered.protein_accession,
                   filtered.genome_accession,
                   MAX(IFNULL(ncbi_exons.target_accession,detected.target_accession)) as 'major_contig',
                   MAX(CASE WHEN ncbi_exons.target_accession IS NULL THEN 'hmm-detected' ELSE 'ncbi-reference' END) as 'proteome_type',
                   MAX(protein_names.name) as 'protein_name'
              FROM (%s) as filtered
         LEFT JOIN needle.ncbi_exons ON filtered.protein_accession = ncbi_exons.protein_hit_id AND filtered.genome_accession = ncbi_exons.genome_accession
         LEFT JOIN needle.detected ON filtered.protein_accession = detected.protein_hit_id AND filtered.genome_accession = detected.genome_accession
         LEFT JOIN needle.protein_names ON filtered.protein_accession = protein_names.protein_accession
             GROUP BY filtered.protein_accession, filtered.genome_accession
          """ % self.selection_sql

        df = self.con.sql(sql).df()
        return df

    def ncbi_fragments(self):
        sql = """
            SELECT DISTINCT ncbi_exons.protein_hit_id,
                   ncbi_exons.genome_accession,
                   target_accession,
                   target_start,
                   target_end,
                   query_accession,
                   query_start,
                   query_end  
              FROM (%s) as filtered
              JOIN needle.ncbi_exons ON filtered.protein_accession = ncbi_exons.protein_hit_id AND filtered.genome_accession = ncbi_exons.genome_accession
          """ % self.selection_sql

        df = self.con.sql(sql).df()
        return df

    def detected_fragments(self):
        sql = """
            SELECT DISTINCT detected.protein_hit_id,
                   detected.genome_accession,
                   target_accession,
                   target_start,
                   target_end,
                   query_accession,
                   query_start,
                   query_end  
              FROM (%s) as filtered
              JOIN needle.detected ON filtered.protein_accession = detected.protein_hit_id AND filtered.genome_accession = detected.genome_accession
          """ % self.selection_sql

        df = self.con.sql(sql).df()
        return df


class ProteinMatches(object):

    def __init__(self, ko_id, protein_accession, genome_accession, ko_match_df, pfam_match_df, sequence):
        self.ko_id = ko_id
        self.protein_accession = protein_accession
        self.genome_accession = genome_accession
        self.sequence = sequence
        self.ko_match_df = ko_match_df
        self.pfam_match_df = pfam_match_df
    
    def ko_match_at(self, protein_pos):
        df = self.ko_match_df[(self.ko_match_df['ko_protein_start'] <= protein_pos) & \
                              (self.ko_match_df['ko_protein_end'] >= protein_pos)]

        if df.size == 0:
            return None
        res = df.sort_values(by='ko_evalue', ascending=True).iloc[0]
        return res

    def pfam_match_at(self, protein_pos):
        df = self.pfam_match_df[(self.pfam_match_df['pfam_protein_start'] <= protein_pos) & \
                                (self.pfam_match_df['pfam_protein_end'] >= protein_pos)]

        if df.size == 0:
            return None
        res = df.sort_values(by='pfam_evalue', ascending=True).iloc[0]
        return res


class ClusterProteins(object):

    def __init__(self, ko_id, cluster_id, parent_cluster_id, group_df):
        self.ko_id = ko_id
        self.cluster_id = cluster_id
        self.parent_cluster_id = parent_cluster_id
        self.df = group_df

    def cluster_fn(self):
        return f"{self.parent_cluster_id}_{self.cluster_id}"

    def proteins(self, protein_sequence_dict):

        con = duckdb.connect(":default:")
        pfam_df = con.sql("""
            SELECT pfam_match.protein_accession,
                   pfam_match.hmm_accession as 'pfam_accession',
                   pfam_match.hmm_start as 'pfam_hmm_start',
                   pfam_match.hmm_end as 'pfam_hmm_end',
                   pfam_match.protein_start as 'pfam_protein_start',
                   pfam_match.protein_end as 'pfam_protein_end',
                   pfam_match.dom_evalue as 'pfam_evalue'
              FROM needle.clusters
              JOIN needle.pfam_match ON pfam_match.protein_accession = clusters.member_accession
             WHERE clusters.cluster_id = '%s'
        """ % self.cluster_id).df()

        proteins = []
        for (protein_accession, genome_accession), protein_ko_match_df in self.df.groupby(['protein_accession', 'genome_accession']):
            if protein_sequence_dict is None or protein_accession not in protein_sequence_dict:
                seq = None
            else:
                seq = protein_sequence_dict[protein_accession]

            protein_pfam_match_df = pfam_df[pfam_df['protein_accession'] == protein_accession]
            proteins.append(ProteinMatches(self.ko_id, protein_accession, genome_accession, protein_ko_match_df, protein_pfam_match_df, seq))
        return proteins


class AssignedClusters(object):

    def __init__(self, ko_id):
        self.ko_id = ko_id

        con = duckdb.connect(":default:")

	# do the join with ko_match, so we can filter to clusters for a
	# specific KO ID. and might as well get the KO matches now.

        self.df = con.sql("""
            SELECT clusters.cluster_id,
                   clusters.parent_cluster_id,
                   ko_match.protein_accession,
                   ko_match.genome_accession,
                   ko_match.hmm_accession as 'ko_accession',
                   ko_match.hmm_start as 'ko_hmm_start',
                   ko_match.hmm_end as 'ko_hmm_end',
                   ko_match.protein_start as 'ko_protein_start',
                   ko_match.protein_end as 'ko_protein_end',
                   ko_match.dom_evalue as 'ko_evalue',
                   ko_match.dom_score as 'ko_score'
              FROM needle.clusters
              JOIN needle.ko_match ON ko_match.protein_accession = clusters.member_accession AND ko_match.hmm_accession = '%s'
        """ % ko_id).df()

    def clusters(self):

        clusters = []
        for (cluster_id, parent_cluster_id), group_df in self.df.groupby(['cluster_id', 'parent_cluster_id']):
            clusters.append(ClusterProteins(self.ko_id, cluster_id, parent_cluster_id, group_df))
        return clusters


class ToUpdate_PutativeProteins(object):

    def __init__(self, ko_id, evalue_threshold=1E-100):
        self.ko_id = ko_id

        con = duckdb.connect(":default:")
        self.df = con.sql("""
            SELECT ko_match.protein_accession,
                   ko_match.genome_accession
              FROM needle.ko_match
             WHERE ko_match.hmm_accession = '%s'
             GROUP BY ko_match.protein_accession, ko_match.genome_accession
             HAVING MIN(ko_match.dom_evalue) < %s
        """ % (ko_id, evalue_threshold)).df()

    def protein_genome_accessions(self):
        unique_pairs = self.df[['protein_accession','genome_accession']].drop_duplicates()
        return list(unique_pairs.itertuples(index=False, name=None))
