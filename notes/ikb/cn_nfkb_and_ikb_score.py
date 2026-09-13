#!/usr/bin/env python3

# Requirements:
# Usage: python3 play/cn_nfkb_and_ikb_score.py DIRECTORY
# Read *.zip directly inside DIRECTORY (not recursively). Other paths are fixed
# relative to this script: cn_nfkb_and_ikb_pairs.tsv, res-nfkb-ikb/input.faa,
# and output/cache cn_nfkb_and_ikb_scores.tsv.
# Match fold_ACCESSION[_ACCESSION].zip against original TSV accessions, lowercased
# with periods replaced by underscores. Reject unmatched/ambiguous names and
# multiple archives for one TSV row. Genome identity comes from that unique row.
# Paired folds must have NFkB copies on A/B and IkB on C. Select the inclusive
# TSV contact window on A/B; score all three full and regional chain pairs.
# Standalone folds must have only chain A. Compare its TSV window to the inclusive
# interval [TSV end + 20, original FASTA length]. Reject empty/absent intervals.
# Reuse sieve ZIP parsing, residue/PAE mapping, contact detection and pDockQ2.
# Same-chain scoring applies the existing formula to disjoint residue selections;
# this does not establish scientific calibration for intramolecular interactions.
# Require five models. Each output score is the maximum across all five models
# of the existing maximum of the two directional pDockQ2 scores (8-Angstrom cutoff).
# Preserve TSV rows/order and original columns, appending six paired score columns,
# A_region_A_Cterm_pDockQ2_max, and zipfile/zipfile_sha256/scoring_signature.
# Leave inapplicable scores and rows without a current ZIP blank (zero means a
# computed zero). Rebuild from current inputs, replacing output atomically only
# after success. Cache hits require ZIP content hash and scoring signature to
# match, plus complete finite applicable scores. Signature includes intervals,
# chain layout, model count, cutoff and explicit scoring version. Bump VERSION
# when scoring semantics change, including changes in the shared calculation.
# Cache hashing and atomic TSV I/O are shared with alphafold-pdockq2-directory.py.

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import zipfile

from Bio import SeqIO

# Use the adjacent source checkout, including its shared scoring extensions.
BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parent / "sieve"))
from sieve.alphafold_cache import read_cache, sha256_file, write_cache
from sieve.alphafold_pdockq2 import score_zip


PAIRS = BASE / "cn_nfkb_and_ikb_pairs.tsv"
FASTA = BASE / "res-nfkb-ikb/input.faa"
OUTPUT = BASE / "cn_nfkb_and_ikb_scores.tsv"
VERSION = 1
INPUT_COLUMNS = ["Genome accession", "NFkB candidate", "IkB candidate",
                 "NFkB chain contact coordinate"]
PAIRED_COLUMNS = ["A_B_pDockQ2_max", "A_C_pDockQ2_max", "B_C_pDockQ2_max",
                  "A_region_B_region_pDockQ2_max", "A_region_C_pDockQ2_max",
                  "B_region_C_pDockQ2_max"]
SINGLE_COLUMN = "A_region_A_Cterm_pDockQ2_max"
SCORE_COLUMNS = PAIRED_COLUMNS + [SINGLE_COLUMN]
META_COLUMNS = ["zipfile", "zipfile_sha256", "scoring_signature"]


def filename(row):
    accessions = [row["NFkB candidate"], row["IkB candidate"]]
    return "fold_" + "_".join(a.lower().replace(".", "_") for a in accessions if a) + ".zip"


def configuration(row, lengths):
    try:
        start, end = map(int, row[INPUT_COLUMNS[3]].split(":"))
    except (AttributeError, ValueError):
        raise ValueError(f"invalid contact coordinate: {row[INPUT_COLUMNS[3]]!r}") from None
    if not 1 <= start <= end:
        raise ValueError(f"invalid contact window: {start}:{end}")
    nfkb = row["NFkB candidate"]
    if nfkb not in lengths or end > lengths[nfkb]:
        raise ValueError(f"{nfkb}: missing sequence or contact exceeds original length")
    config = {"version": VERSION, "cutoff": 8.0, "expected_models": 5}
    if row["IkB candidate"]:
        config.update(regions={"A": [(start, end)], "B": [(start, end)]},
                      expected_chains=["A", "B", "C"])
    else:
        cterm = end + 20
        if cterm > lengths[nfkb]:
            raise ValueError(f"{nfkb}: C-terminal interval {cterm}:{lengths[nfkb]} is empty")
        config.update(comparisons=[("A", [(start, end)], "A", [(cterm, lengths[nfkb])])],
                      expected_chains=["A"])
    return config


def aggregate(rows, paired):
    scores = {}
    for row in rows:
        if not paired:
            column = SINGLE_COLUMN
        else:
            first, second = sorted((row["chain 1"], row["chain 2"]))
            if row["scope"] == "regional":
                first += "_region" if first in {"A", "B"} else ""
                second += "_region" if second in {"A", "B"} else ""
            column = f"{first}_{second}_pDockQ2_max"
        value = row["pDockQ2 max"]
        if not math.isfinite(value):
            raise ValueError("non-finite pDockQ2 score")
        scores[column] = max(scores.get(column, value), value)
    expected = PAIRED_COLUMNS if paired else [SINGLE_COLUMN]
    if set(scores) != set(expected):
        raise ValueError(f"unexpected score columns: {sorted(scores)}")
    return {column: f"{value:.6f}" for column, value in scores.items()}


def run(directory):
    columns, pairs = read_cache(PAIRS, INPUT_COLUMNS)
    if columns is None:
        raise ValueError(f"missing or empty pairs TSV: {PAIRS}")
    if set(columns) & set(SCORE_COLUMNS + META_COLUMNS):
        raise ValueError("input pairs TSV already contains score/cache columns")
    lengths = {}
    for record in SeqIO.parse(FASTA, "fasta"):
        if record.id in lengths or not record.seq:
            raise ValueError(f"duplicate or empty FASTA entry: {record.id}")
        lengths[record.id] = len(record.seq)
    names = {}
    identities = set()
    for index, row in enumerate(pairs):
        if any(row.get(c) is None for c in columns) or any(not row[c] for c in INPUT_COLUMNS if c != "IkB candidate"):
            raise ValueError(f"incomplete pairs TSV row {index + 2}")
        identity = tuple(row[c] for c in INPUT_COLUMNS[:3])
        if identity in identities:
            raise ValueError(f"duplicate pairs TSV row: {identity}")
        identities.add(identity)
        names.setdefault(filename(row), []).append(index)
    paths = sorted(directory.glob("*.zip"))
    if not paths:
        raise ValueError(f"directory contains no ZIP files: {directory}")
    jobs = []
    matched = set()
    for path in paths:
        indices = names.get(path.name.lower(), [])
        if len(indices) != 1:
            raise ValueError(f"{path.name}: expected one matching TSV row, found {len(indices)}")
        index = indices[0]
        if index in matched:
            raise ValueError(f"multiple ZIPs match TSV row {index + 2}")
        matched.add(index)
        jobs.append((path, index, configuration(pairs[index], lengths)))
    output_columns = columns + SCORE_COLUMNS + META_COLUMNS
    cached_columns, cached_rows = read_cache(OUTPUT, output_columns)
    cache = {}
    for row in cached_rows:
        key = (row["zipfile_sha256"], row["scoring_signature"])
        if all(key):
            cache[key] = row
    output_rows = [dict(row) for row in pairs]
    for path, index, config in jobs:
        digest = sha256_file(path)
        signature = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
        cached = cache.get((digest, signature), {})
        paired = bool(pairs[index]["IkB candidate"])
        expected = PAIRED_COLUMNS if paired else [SINGLE_COLUMN]
        try:
            valid = all(math.isfinite(float(cached.get(column, ""))) for column in expected)
        except (ValueError, TypeError):
            valid = False
        if valid:
            print(f"Skipping {path.name}: cached", file=sys.stderr)
            scores = {column: cached[column] for column in expected}
        else:
            print(f"Processing {path.name}", file=sys.stderr)
            options = {key: value for key, value in config.items() if key != "version"}
            scores = aggregate(score_zip(path, **options), paired)
        output_rows[index].update(scores, zipfile=path.name, zipfile_sha256=digest,
                                  scoring_signature=signature)
    write_cache(OUTPUT, output_columns, output_rows)
    print(f"Wrote {OUTPUT}: {len(jobs)} folds, {len(pairs) - len(jobs)} rows without ZIPs")


def main():
    parser = argparse.ArgumentParser(description="Score NFkB/IkB AlphaFold ZIPs using TSV contact regions")
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    try:
        run(args.directory)
    except (OSError, ValueError, zipfile.BadZipFile) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
