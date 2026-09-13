#!/usr/bin/env python3

# Generate genome-specific NFkB/IkB candidates for AlphaFold testing.
# All paths are fixed relative to this script, independent of the working directory.
# Inputs:
#   cn_nfkb_and_ikb_loci.tsv: protein accession, genome accession, contig accession,
#       locus start 1b, locus end 1b; repeated feature rows describe the same locus.
#   res-nfkb-ikb/input.faa: original protein sequences, keyed by FASTA accession.
#   cn_truncations.faa: optional replacements named ORIGINAL_ACCESSION_trunc.
#   protein_pfam.tsv: query_accession/query_database identify protein/genome;
#       target_accession identifies Pfam; query_start/query_end are inclusive AA
#       coordinates. Ignore Pfam version suffixes and stream this large table.
# Deduplication:
#   Collapse repeated locus rows and normalize reverse-strand coordinates.
#   Group by genome and contig, sort by locus start, and merge inclusive overlaps
#   transitively in one pass. Strand does not affect overlap. Retain the longest
#   ORIGINAL protein sequence per merged locus; break ties by accession.
# Classification (after locus deduplication):
#   p105-like: PF00554 and PF12796; p50-like: PF00554 without PF12796;
#   IkB: PF12796 without PF00554. Candidates with neither are not emitted.
# Pairing:
#   Emit each p105 alone, every p105/IkB pair, and every p50/IkB pair within a
#   genome. Do not emit standalone p50 or IkB candidates or cross-genome pairs.
# Sequences:
#   Standalone p105 uses its original sequence/accession.
#   Paired p105 ends immediately AFTER its first GGG (retain GGG), with accession
#   ORIGINAL_Cterm_GGG_trunc. Raise an error if that paired p105 lacks GGG.
#   p50 uses its original sequence/accession. IkB uses ORIGINAL_trunc from the
#   truncation FASTA when present, otherwise its original sequence/accession.
# Contacts:
#   NFkB must have exactly one distinct PF16179 interval. Its contact window is
#   one inclusive, 1-based AA interval: domain start - 15 through domain end + 25.
#   Do not clip the window. Error if it extends outside the selected sequence,
#   including a paired p105 truncated at GGG. Duplicate identical hits are allowed.
# Outputs:
#   cn_nfkb_and_ikb_pairs.tsv: Genome accession, NFkB candidate, IkB candidate,
#       NFkB chain contact coordinate (formatted start:end). TSV accessions are
#       always ORIGINAL; the IkB field is blank for standalone p105.
#   cn_nfkb_and_ikb_pairs.faa: a sequence-free >pair description containing genome
#       and selected accession(s), followed by one (standalone) or two separate
#       FASTA records. Never concatenate chains. Truncation suffixes appear here.
# Validation:
#   Reject missing columns/values, locus error rows, conflicting locus metadata,
#   invalid coordinates, missing original sequences, and empty/duplicate FASTA
#   records. Validate all generated groups before opening either output file.
#   Output ordering is deterministic (genome, NFkB accession, IkB accession),
#   with each standalone p105 preceding its pairs. Existing outputs are replaced.
# No automated tests are included; the user will test this script.

import csv
from pathlib import Path


BASE = Path(__file__).resolve().parent
LOCI = BASE / "cn_nfkb_and_ikb_loci.tsv"
PROTEINS = BASE / "res-nfkb-ikb/input.faa"
TRUNCATIONS = BASE / "cn_truncations.faa"
PFAM = BASE / "protein_pfam.tsv"
OUTPUT_TSV = BASE / "cn_nfkb_and_ikb_pairs.tsv"
OUTPUT_FASTA = BASE / "cn_nfkb_and_ikb_pairs.faa"


def read_fasta(path):
    sequences = {}
    accession = None
    parts = []

    def store():
        if accession is not None:
            if accession in sequences or not parts:
                raise ValueError(f"{path}: duplicate or empty FASTA record {accession}")
            sequences[accession] = "".join(parts)

    with path.open() as stream:
        for raw in stream:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                if not line[1:].strip():
                    raise ValueError(f"{path}: empty FASTA header")
                accession = line[1:].split()[0]
                parts = []
            else:
                if accession is None:
                    raise ValueError(f"{path}: sequence before FASTA header")
                parts.append("".join(line.split()).upper())
    store()
    return sequences


def rows(path, required):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for row in reader:
            if any(not row.get(field) for field in required):
                raise ValueError(f"{path}:{reader.line_num}: missing required value")
            yield row


def interval(start, end):
    start, end = sorted((int(start), int(end)))
    if start < 1:
        raise ValueError(f"Invalid 1-based coordinates: {start}, {end}")
    return start, end


def retained_candidates(sequences):
    loci = {}
    required = ["protein accession", "genome accession", "contig accession",
                "locus start 1b", "locus end 1b"]
    for row in rows(LOCI, required):
        accession, genome = row["protein accession"], row["genome accession"]
        if row.get("error"):
            raise ValueError(f"{accession}: locus error: {row['error']}")
        if accession not in sequences:
            raise ValueError(f"{genome}/{accession}: missing original sequence")
        start, end = interval(row["locus start 1b"], row["locus end 1b"])
        key = (genome, accession)
        locus = (row["contig accession"], start, end)
        if key in loci and loci[key] != locus:
            raise ValueError(f"{key}: conflicting locus metadata")
        loci[key] = locus

    contigs = {}
    for (genome, accession), (contig, start, end) in loci.items():
        contigs.setdefault((genome, contig), []).append((start, end, accession))
    retained = set()
    for (genome, contig), candidates in sorted(contigs.items()):
        cluster_end = 0
        best = None
        for start, end, accession in sorted(candidates):
            if best is None or start > cluster_end:
                if best is not None:
                    retained.add((genome, best))
                best, cluster_end = accession, end
            else:
                cluster_end = max(cluster_end, end)
                if (-len(sequences[accession]), accession) < (-len(sequences[best]), best):
                    best = accession
        if best is not None:
            retained.add((genome, best))
    return retained


def read_domains(candidates):
    domains = {key: {} for key in candidates}
    required = ["query_accession", "query_database", "target_accession",
                "query_start", "query_end"]
    for row in rows(PFAM, required):
        key = (row["query_database"], row["query_accession"])
        pfam = row["target_accession"].split(".")[0]
        if key in domains and pfam in {"PF00554", "PF12796", "PF16179"}:
            region = interval(row["query_start"], row["query_end"])
            domains[key].setdefault(pfam, set()).add(region)
    return domains


def generate_groups(sequences, truncations, domains):
    genomes = {}
    for (genome, accession), hits in sorted(domains.items()):
        classes = genomes.setdefault(genome, {"nfkb": [], "ikb": []})
        if "PF00554" in hits:
            classes["nfkb"].append(accession)
        elif "PF12796" in hits:
            classes["ikb"].append(accession)

    groups = []
    for genome, classes in sorted(genomes.items()):
        for nfkb in classes["nfkb"]:
            hits = domains[(genome, nfkb)]
            contacts = hits.get("PF16179", set())
            if len(contacts) != 1:
                raise ValueError(f"{genome}/{nfkb}: expected one PF16179 interval, got {sorted(contacts)}")
            start, end = next(iter(contacts))
            start, end = start - 15, end + 25
            p105 = "PF12796" in hits
            partners = ([None] if p105 else []) + classes["ikb"]
            for ikb in partners:
                nfkb_id, nfkb_seq = nfkb, sequences[nfkb]
                if p105 and ikb is not None:
                    cut = nfkb_seq.find("GGG")
                    if cut == -1:
                        raise ValueError(f"{genome}/{nfkb}: paired p105 lacks GGG")
                    nfkb_id += "_Cterm_GGG_trunc"
                    nfkb_seq = nfkb_seq[:cut + 3]
                if start < 1 or end > len(nfkb_seq):
                    raise ValueError(f"{genome}/{nfkb_id}: contact {start}:{end} outside sequence length {len(nfkb_seq)}")
                records = [(nfkb_id, nfkb_seq)]
                if ikb is not None:
                    truncated_id = ikb + "_trunc"
                    records.append((truncated_id, truncations[truncated_id])
                                   if truncated_id in truncations else (ikb, sequences[ikb]))
                groups.append(((genome, nfkb, ikb or "", f"{start}:{end}"), records))
    return groups


def main():
    sequences = read_fasta(PROTEINS)
    truncations = read_fasta(TRUNCATIONS)
    candidates = retained_candidates(sequences)
    domains = read_domains(candidates)
    groups = generate_groups(sequences, truncations, domains)
    with OUTPUT_TSV.open("w", newline="") as tsv, OUTPUT_FASTA.open("w") as fasta:
        writer = csv.writer(tsv, delimiter="\t", lineterminator="\n")
        writer.writerow(["Genome accession", "NFkB candidate", "IkB candidate",
                         "NFkB chain contact coordinate"])
        for row, records in groups:
            writer.writerow(row)
            fasta.write(f">pair {row[0]} {'-'.join(accession.replace(".","_") for accession, _ in records)}\n")
            for accession, sequence in records:
                fasta.write(f">{accession}\n{sequence}\n")
    print(f"Wrote {len(groups)} groups from {len(candidates)} retained candidates")
    print(OUTPUT_TSV)
    print(OUTPUT_FASTA)


if __name__ == "__main__":
    main()
