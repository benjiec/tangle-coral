import csv
from typing import Dict, List, Optional
from Bio.Seq import Seq

from .match import Match
from .seq import extract_subsequence, extract_subsequence_strand_sensitive, read_fasta_as_dict
from .seq import to_dna_coordinate, compute_three_frame_translations
from .hmm import hmmsearch_sequence_dict


class Results:
    # NCBI-style canonical headers (preferred)
    H_QSEQID = "qseqid"
    H_SSEQID = "sseqid"
    H_EVALUE = "evalue"
    H_PIDENT = "pident"
    H_QSTART = "qstart"
    H_QEND = "qend"
    H_SSTART = "sstart"
    H_SEND = "send"
    H_SSEQ = "sseq"

    # Producer header order (used by blast-genome.py). Keep minimal set we rely on.
    PRODUCER_HEADER = [
        H_QSEQID,
        H_SSEQID,
        H_EVALUE,
        H_PIDENT,
        H_QSTART,
        H_QEND,
        H_SSTART,
        H_SEND,
        H_SSEQ,
    ]

    # Raw BLAST outfmt fields (as configured in run_blast_search)
    RAW_OUTFMT_FIELDS = [
        "qseqid",
        "sseqid",
        "evalue",
        "pident",
        "qstart",
        "qend",
        "stitle",  # not included in output header, but present in raw
        "sseq",
        "sstart",
        "send",
    ]

    # Mapping from raw BLAST field names to our header names
    RAW_TO_HEADER = {
        "qseqid": H_QSEQID,
        "sseqid": H_SSEQID,
        "evalue": H_EVALUE,
        "pident": H_PIDENT,
        "qstart": H_QSTART,
        "qend": H_QEND,
        "sseq": H_SSEQ,
        "sstart": H_SSTART,
        "send": H_SEND,
    }
    # Inverse mapping: header -> raw key
    HEADER_TO_RAW = {v: k for k, v in RAW_TO_HEADER.items()}

    def __init__(
        self,
        results_tsv_path: str,
        target_fasta_path: Optional[str] = None
    ) -> None:
        self._results_tsv_path = results_tsv_path
        self._target_fasta_path = target_fasta_path
        self._matches: Optional[List[Match]] = None
        self._target_sequences_by_accession: Optional[Dict[str, str]] = None

    def matches(self) -> List[Match]:
        if self._matches is None:
            self._parse_once()
        # mypy: _matches is now set
        return self._matches or []

    # ----- Internal helpers -----

    def _parse_once(self) -> None:
        if self._target_fasta_path:
            self._target_sequences_by_accession = read_fasta_as_dict(self._target_fasta_path)

        with open(self._results_tsv_path, "r") as tsv_file:
            reader = csv.reader(tsv_file, delimiter="\t")
            header_row = next(reader, None)
            if header_row is None:
                self._matches = []
                return
            header_index = self._build_header_index(header_row)

            matches: List[Match] = []
            for row in reader:
                if not row or all(not cell for cell in row):
                    continue
                match = self._row_to_match(row, header_index)
                # Attach sequences if requested and available
                if match is not None:
                    if self._target_sequences_by_accession is not None:
                        seq = extract_subsequence_strand_sensitive(
                            self._target_sequences_by_accession.get(match.target_accession, None),
                            match.target_start,
                            match.target_end,
                        )
                        match.target_sequence = seq
                    matches.append(match)

            self._matches = matches

    def _build_header_index(self, header_row: List[str]) -> Dict[str, int]:
        # Normalize header names by stripping whitespace
        normalized = [h.strip() for h in header_row]
        index: Dict[str, int] = {name: i for i, name in enumerate(normalized)}

        # Build a mapping for canonical names only; sseq is optional
        mapping: Dict[str, Optional[int]] = {
            self.H_QSEQID: index.get(self.H_QSEQID),
            self.H_SSEQID: index.get(self.H_SSEQID),
            self.H_EVALUE: index.get(self.H_EVALUE),
            self.H_PIDENT: index.get(self.H_PIDENT),
            self.H_QSTART: index.get(self.H_QSTART),
            self.H_QEND: index.get(self.H_QEND),
            self.H_SSTART: index.get(self.H_SSTART),
            self.H_SEND: index.get(self.H_SEND),
            self.H_SSEQ: index.get(self.H_SSEQ),
        }

        # Ensure required fields exist
        required = [
            self.H_QSEQID,
            self.H_SSEQID,
            self.H_EVALUE,
            self.H_PIDENT,
            self.H_QSTART,
            self.H_QEND,
            self.H_SSTART,
            self.H_SEND,
        ]
        missing = [key for key in required if mapping.get(key) is None]
        if missing:
            raise ValueError(f"Missing required columns in results TSV: {', '.join(missing)}")

        # Convert Optional[int] to int where present
        return {k: v for k, v in mapping.items() if v is not None}

    def _row_to_match(self, row: List[str], header_index: Dict[str, int]) -> Optional[Match]:
        qacc = row[header_index[self.H_QSEQID]].strip()
        sacc = row[header_index[self.H_SSEQID]].strip()
        evalue_str = row[header_index[self.H_EVALUE]].strip()
        pident_str = row[header_index[self.H_PIDENT]].strip()
        qstart_str = row[header_index[self.H_QSTART]].strip()
        qend_str = row[header_index[self.H_QEND]].strip()
        sstart_str = row[header_index[self.H_SSTART]].strip()
        send_str = row[header_index[self.H_SEND]].strip()

        matched_seq = None
        if self.H_SSEQ in header_index and header_index[self.H_SSEQ] < len(row):
            cell = row[header_index[self.H_SSEQ]].strip()
            matched_seq = cell if cell != "" else None

        qstart = int(qstart_str)
        qend = int(qend_str)
        if qstart > qend:
            raise ValueError(f"qstart ({qstart}) must be <= qend ({qend}) in results TSV row: {row}")

        sstart = int(sstart_str)
        send = int(send_str)

        match = Match(
            query_accession=qacc,
            target_accession=sacc,
            query_start=qstart,
            query_end=qend,
            target_start=sstart,
            target_end=send,
            e_value=float(evalue_str.replace(",", "")) if evalue_str else None,
            identity=float(pident_str.replace(",", "")) if pident_str else None,
            matched_sequence=matched_seq,
        )
        return match


def hmm_search_genome_sequence(
    hmm_file, target_accession, target_sequence,
    win, win_overlap,
    contig_left_0b = None, contig_right_excl_0b = None
):
    """
    Returns array of lists, each list has values in Results.PRODUCER_HEADER order
    """

    if contig_left_0b is None:
        contig_left_0b = 0
    if contig_right_excl_0b is None:
        contig_right_excl_0b = len(target_sequence)

    detected = []

    for win_i in range(contig_left_0b, contig_right_excl_0b, win):
        translated_fasta = {}

        subs = target_sequence[win_i:min(win_i+win+win_overlap, contig_right_excl_0b)]
        translations_fwd = compute_three_frame_translations(subs, 1, len(subs))
        translations_rev = compute_three_frame_translations(subs, len(subs), 1)

        if translations_fwd[0][2]:
            translated_fasta[f"{target_accession}_fwd_0"] = translations_fwd[0][2]
        if translations_fwd[1][2]:
            translated_fasta[f"{target_accession}_fwd_1"] = translations_fwd[1][2]
        if translations_fwd[2][2]:
            translated_fasta[f"{target_accession}_fwd_2"] = translations_fwd[2][2]
        if translations_rev[0][2]:
            translated_fasta[f"{target_accession}_rev_0"] = translations_rev[0][2]
        if translations_rev[1][2]:
            translated_fasta[f"{target_accession}_rev_1"] = translations_rev[1][2]
        if translations_rev[2][2]:
            translated_fasta[f"{target_accession}_rev_2"] = translations_rev[2][2]
        
        for k,v in translated_fasta.items():
            # print(f">{k}\n{v}")
            pass
        if len(translated_fasta.keys()) == 0:
            continue

        hmm_rows = hmmsearch_sequence_dict(hmm_file, translated_fasta)
        hmm_rows = [row for row in hmm_rows if row["evalue"] <= 0.001]

        for row in hmm_rows:
            query_accession = row["query_accession"]
            if not query_accession or query_accession.strip() == "-":
                query_accession = row["query_name"]

            frame = int(row["target_name"][-1])
            if "_fwd_" in row["target_name"]:
                frame_dna_start = translations_fwd[frame][0]
                frame_dna_end = translations_fwd[frame][1]
                aa_seq = extract_subsequence(translations_fwd[frame][2], row["ali_from"], row["ali_to"])
            else:
                frame_dna_start = translations_rev[frame][0]
                frame_dna_end = translations_rev[frame][1]
                aa_seq = extract_subsequence(translations_rev[frame][2], row["ali_from"], row["ali_to"])

            assert (abs(frame_dna_end-frame_dna_start)+1)%3 == 0

            dna_ali_from, dna_ali_to = to_dna_coordinate(frame_dna_start+win_i, frame_dna_end+win_i, row["ali_from"], row["ali_to"])

            out = [
              query_accession,
              target_accession,
              row["evalue"],
              '',
              row["hmm_from"],
              row["hmm_to"],
              dna_ali_from,
              dna_ali_to,
              aa_seq
            ]

            detected.append(out)

    return detected
