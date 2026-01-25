import os
from pathlib import Path


class DefaultPath(object):

    @staticmethod
    def ncbi_download_dir():
        return "ncbi-downloads"

    @staticmethod
    def ncbi_genome_dir(genome_accession):
        return os.path.join(DefaultPath.ncbi_download_dir(), "ncbi_dataset/data", genome_accession)

    @staticmethod
    def ncbi_genome_protein_faa(genome_accession):
        return os.path.join(DefaultPath.ncbi_genome_dir(genome_accession), "protein.faa")

    @staticmethod
    def ncbi_genome_gff(genome_accession):
        return os.path.join(DefaultPath.ncbi_genome_dir(genome_accession), "genomic.gff")

    @staticmethod
    def ncbi_genome_fna(genome_accession):
        ncbi_download_dir = os.path.join(DefaultPath.ncbi_genome_dir(genome_accession))
        fna_file = None
        for pattern in ["*.fna", "*.fasta"]:
            matches = list(Path(ncbi_download_dir).glob(pattern))
            if matches:
                fna_file = str(matches[0])
        return fna_file

    @staticmethod
    def pfam_hmm():
        return "pfam-downloads/Pfam-A.hmm"

    @staticmethod
    def ko_hmm():
        return "kegg-downloads/ko.hmm"

    @staticmethod
    def module_detected_proteins(module):
        return f"data/{module.lower()}_results/proteins.faa"

    @staticmethod
    def module_ko_assigned_proteins(module, ko):
        return f"data/{module.lower()}_results/faa/{ko.upper()}.faa"

    @staticmethod
    def module_ko_putative_proteins(module, ko):
        return f"data/{module.lower()}_results/faa/{ko.upper()}-putative.faa"

    @staticmethod
    def module_cluster_dir(module):
        return f"data/{module.lower()}_results/clusters"

    @staticmethod
    def module_cluster_tsv(module):
        return f"data/{module.lower()}_results/clusters.tsv"
