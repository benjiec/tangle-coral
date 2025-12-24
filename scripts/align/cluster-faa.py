import os
import csv
import uuid
import hashlib
import argparse
import subprocess
from pathlib import Path
from needle.match import read_fasta_as_dict


COL_PARENT = "Parent Cluster ID"
COL_CLUSTER_ID = "Cluster ID"
COL_REP_ACC = "Representative Accession"
COL_MEM_ACC = "Member Accession"
TSV_FIELDS = [COL_PARENT, COL_CLUSTER_ID, COL_REP_ACC, COL_MEM_ACC]


def create_cluster_tsv(tsv_fn):

    if os.path.exists(tsv_fn):
        return

    with open(tsv_fn, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=TSV_FIELDS, delimiter='\t')
        writer.writeheader()


def format_cluster_tsv_row(faa_file, rep_acc, mem_acc):

    faa_name = Path(faa_file).stem
    cluster_id = str(hashlib.sha1((faa_name+"_"+rep_acc).encode()).hexdigest())[:8]

    return {
        COL_PARENT: faa_name,
        COL_CLUSTER_ID: cluster_id,
        COL_REP_ACC: rep_acc,
        COL_MEM_ACC: mem_acc
    }


def run_mmseqs_cluster(faa_file, output_result_prefix):

    cmd = ["scripts/mmseqs-cluster", faa_file, output_result_prefix]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result


def parse_mmseqs_cluster_results(mmseqs_cluster_tsv):

    cluster_reps = {}
    with open(mmseqs_cluster_tsv, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if line:
                rep, mem = line.split("\t")
                cluster_reps.setdefault(rep, []).append(mem)

    return cluster_reps


def append_to_cluster_tsv(faa_file, output_result_prefix, cluster_tsv_fn):

    cluster_reps = parse_mmseqs_cluster_results(output_result_prefix+"_cluster.tsv")

    data = []
    cluster_ids = {}

    for rep, members in cluster_reps.items():
        for member in members:
            d = format_cluster_tsv_row(faa_file, rep, member)
            data.append(d)
            cluster_ids.setdefault(d[COL_CLUSTER_ID], []).append(d)

    with open(cluster_tsv_fn, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=TSV_FIELDS, delimiter='\t')
        for d in data:
            writer.writerow(d)

    return cluster_ids


def create_cluster_faa_files(faa_file, cluster_ids):
    sequences = read_fasta_as_dict(faa_file)

    for cluster_id, members in cluster_ids.items():
        p = Path(faa_file)
        fn = str(p.parent)+'/'+p.stem+"."+cluster_id+".cluster.faa"
        with open(fn, "w") as f:
            for member in members:
                member_seq = sequences[member[COL_MEM_ACC]]
                f.write(">"+member[COL_MEM_ACC]+"\n"+member_seq+"\n")


def cleanup_mmseqs_outputs(output_result_prefix):
    suffixes = ["_cluster.tsv", "_rep_seq.fasta", "_all_seqs.fasta"]
    for sf in suffixes:
        fn = output_result_prefix+sf
        if os.path.exists(fn):
            os.remove(fn)


def cluster_faa(faa_file, cluster_tsv_fn):
    """
    Given faa_file, e.g. data/m00009_results/faa/K00024.faa, cluster sequences
    in the FAA file using mmseqs, then create/append to TSV of clustering
    metadata, and create new faa files with cluster ID suffix.
    """

    create_cluster_tsv(cluster_tsv_fn)
    output_result_prefix = faa_file+"_"+uuid.uuid4().hex[:16] 
    run_mmseqs_cluster(faa_file, output_result_prefix)
    cluster_ids = append_to_cluster_tsv(faa_file, output_result_prefix, cluster_tsv_fn)
    create_cluster_faa_files(faa_file, cluster_ids)
    cleanup_mmseqs_outputs(output_result_prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("query_fasta")
    parser.add_argument("cluster_tsv_fn")
    args = parser.parse_args()

    cluster_faa(args.query_fasta, args.cluster_tsv_fn)
