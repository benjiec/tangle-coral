import os
import tempfile
import argparse
from scripts.defaults import DefaultPath
from needle.seq import read_fasta_as_dict, write_fasta_from_dict
from needle.cluster import append_to_cluster_tsv, create_cluster_faa_files


def cluster_fn(module, ko, cluster_id):
    fn = f"data/{module}_results/clusters/{ko}_{cluster_id}.faa"
    if os.path.exists(fn):
        return fn
    # print(f"{fn} does not exist")
    fn = f"data/{module}_results/clusters/{ko}-putative_{cluster_id}.faa"
    if os.path.exists(fn):
        return fn
    # print(f"{fn} does not exist")
    raise Exception("Cannot find cluster FAA file")


def combine(module, ko, cluster_id_1, cluster_id_2, cluster_dir, cluster_tsv_fn):
    assert os.path.exists(cluster_tsv_fn)

    combined = read_fasta_as_dict(cluster_fn(module, ko, cluster_id_1))
    cluster2 = read_fasta_as_dict(cluster_fn(module, ko, cluster_id_2))
    combined.update(cluster2)

    tmpf = tempfile.NamedTemporaryFile(delete=False, suffix=".faa", mode="w")
    tmpf.close()
    write_fasta_from_dict(combined, tmpf.name)

    combined_id = f"{cluster_id_1}_{cluster_id_2}"
    cluster_reps = {
        combined_id: list(combined.keys())
    }

    cluster_ids = append_to_cluster_tsv(tmpf.name, cluster_reps, cluster_tsv_fn, parent_cluster_id=f"{ko}-manual")
    create_cluster_faa_files(tmpf.name, cluster_dir, cluster_ids, parent_cluster_id=f"{ko}-manual")
    os.remove(tmpf.name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("module")
    parser.add_argument("ko")
    parser.add_argument("cluster_id_1")
    parser.add_argument("cluster_id_2")
    args = parser.parse_args()

    combine(args.module, args.ko, args.cluster_id_1, args.cluster_id_2,
            DefaultPath.module_cluster_dir(args.module),
            DefaultPath.module_cluster_tsv(args.module))
