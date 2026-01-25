import os
import uuid
import argparse
import subprocess
from needle.cluster import append_to_cluster_tsv, create_cluster_tsv, create_cluster_faa_files


def run_mmseqs_cluster(faa_file, output_result_prefix):
    cmd = ["scripts/cluster/mmseqs-cluster", faa_file, output_result_prefix]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result


def parse_mmseqs_cluster_results(mmseqs_cluster_tsv):
    cluster_reps = {}
    if not os.path.exists(mmseqs_cluster_tsv):
        return cluster_reps
    with open(mmseqs_cluster_tsv, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if line:
                rep, mem = line.split("\t")
                cluster_reps.setdefault(rep, []).append(mem)
    return cluster_reps


def cleanup_mmseqs_outputs(output_result_prefix):
    suffixes = ["_cluster.tsv", "_rep_seq.fasta", "_all_seqs.fasta"]
    for sf in suffixes:
        fn = output_result_prefix+sf
        if os.path.exists(fn):
            os.remove(fn)


def cluster_faa(faa_file, cluster_dir, cluster_tsv_fn):
    """
    Given faa_file, e.g. data/m00009_results/faa/K00024.faa, cluster sequences
    in the FAA file using mmseqs, then create/append to TSV of clustering
    metadata, and create new faa files with cluster ID suffix in cluster_dir.
    """

    create_cluster_tsv(cluster_tsv_fn)
    output_result_prefix = faa_file+"_"+uuid.uuid4().hex[:16] 
    run_mmseqs_cluster(faa_file, output_result_prefix)
    cluster_reps = parse_mmseqs_cluster_results(output_result_prefix+"_cluster.tsv")
    cluster_ids = append_to_cluster_tsv(faa_file, cluster_reps, cluster_tsv_fn)
    create_cluster_faa_files(faa_file, cluster_dir, cluster_ids)
    cleanup_mmseqs_outputs(output_result_prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("module")
    parser.add_argument("query_fasta")
    args = parser.parse_args()

    cluster_dir = DefaultPath.module_cluster_dir(args.module)
    cluster_tsv_fn = DefaultPath.module_cluster_tsv(args.module)
    cluster_faa(args.query_fasta, cluster_dir, cluster_tsv_fn)
