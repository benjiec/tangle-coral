import os
import csv
import hashlib
from pathlib import Path
from needle.seq import read_fasta_as_dict


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


def append_to_cluster_tsv(faa_file, cluster_reps, cluster_tsv_fn, parent_cluster_id=None):
    data = []
    cluster_ids = {}

    for rep, members in cluster_reps.items():
        for member in members:
            d = format_cluster_tsv_row(faa_file, rep, member)
            if parent_cluster_id is not None:
                d[COL_PARENT] = parent_cluster_id
            data.append(d)
            cluster_ids.setdefault(d[COL_CLUSTER_ID], []).append(d)

    with open(cluster_tsv_fn, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=TSV_FIELDS, delimiter='\t')
        for d in data:
            writer.writerow(d)

    return cluster_ids


def create_cluster_faa_files(faa_file, cluster_dir, cluster_ids, parent_cluster_id=None):
    sequences = read_fasta_as_dict(faa_file)

    for cluster_id, members in cluster_ids.items():
        p = Path(faa_file)
        if parent_cluster_id is None:
            fn = cluster_dir+'/'+p.stem+"_"+cluster_id+".faa"
        else:
            fn = cluster_dir+'/'+parent_cluster_id+"_"+cluster_id+".faa"
        with open(fn, "w") as f:
            for member in members:
                member_seq = sequences[member[COL_MEM_ACC]]
                f.write(">"+member[COL_MEM_ACC]+"\n"+member_seq+"\n")
