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
    duckdb.execute(f"CREATE TABLE needle.clusters AS SELECT * FROM read_csv_auto('data/{module_id}_results/clusters.tsv', normalize_names=TRUE)")


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
                   ko_match.dom_evalue as 'ko_evalue'
              FROM needle.clusters
              JOIN needle.ko_match ON ko_match.protein_accession = clusters.member_accession AND ko_match.hmm_accession = '%s'
        """ % ko_id).df()

    def clusters(self):

        clusters = []
        for (cluster_id, parent_cluster_id), group_df in self.df.groupby(['cluster_id', 'parent_cluster_id']):
            clusters.append(ClusterProteins(self.ko_id, cluster_id, parent_cluster_id, group_df))
        return clusters


class PutativeProteins(object):

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
