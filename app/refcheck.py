"""Reference FASTA parsing and VCF-vs-reference verification.

Acceptance rules for a reference upload (size limiting happens at the API
layer):

* exactly one FASTA record (a second ``>`` header is rejected),
* the sequence contains the upper-case letters A/C/G/T only,
* no blank/whitespace-only sequence lines,
* the sequence is non-empty.

Verification walks every record of the built BGZF VCF (never a partial scan)
and checks, for each variant:

1. the VCF chromosome name equals the FASTA record name;
2. the REF span [POS, POS+len(REF)-1] lies inside the reference sequence;
3. the reference substring at that span equals the VCF REF allele.

Every failure becomes a :class:`Problem` carrying the variant position and a
machine-readable reason, so the caller can report *where* and *why*.
"""
from __future__ import annotations

import gzip
import hashlib
from dataclasses import dataclass, field

FASTA_NAME = "reference.fasta"

# Bounded error payloads: a 10 MiB VCF can hold far more bad records than any
# client wants in one response.
MAX_PROBLEMS = 100


class FASTAError(ValueError):
    """The uploaded reference FASTA violates the acceptance rules."""


@dataclass(frozen=True)
class FastaInfo:
    name: str
    sequence: str
    sha256: str
    length: int = field(init=False)

    def __post_init__(self):
        object.__setattr__(self, "length", len(self.sequence))


@dataclass(frozen=True)
class Problem:
    record: int          # 1-based index of the variant record in the VCF
    chrom: str
    position: int        # VCF POS (1-based)
    reason: str          # chromosome_mismatch | position_out_of_range | ref_mismatch
    detail: str
    ref: str | None = None
    expected: str | None = None

    def as_dict(self) -> dict:
        return {
            "record": self.record,
            "chrom": self.chrom,
            "position": self.position,
            "reason": self.reason,
            "detail": self.detail,
            "ref": self.ref,
            "expected": self.expected,
        }


@dataclass(frozen=True)
class VerifyResult:
    problems: list
    records_checked: int
    truncated: bool


def parse_fasta(path: str) -> FastaInfo:
    """Strictly parse a single-record, ACGT-only FASTA file."""
    try:
        with open(path, "rb") as fh:
            raw = fh.read()
    except OSError as exc:  # pragma: no cover - defensive
        raise FASTAError(f"cannot read upload: {exc}")
    digest = hashlib.sha256(raw).hexdigest()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise FASTAError(f"upload is not valid UTF-8 text: {exc}")

    name = None
    parts: list[str] = []
    for lineno, line in enumerate(text.splitlines(), 1):
        if not line:
            raise FASTAError(f"line {lineno}: blank lines are not allowed")
        if line.startswith(">"):
            if name is not None:
                raise FASTAError(
                    f"line {lineno}: only a single chromosome/sequence is allowed"
                )
            header = line[1:].strip()
            if not header:
                raise FASTAError(f"line {lineno}: missing sequence name after '>'")
            name = header.split()[0]
            continue
        if name is None:
            raise FASTAError(f"line {lineno}: sequence data before any '>' header")
        stripped = line.strip()
        if not stripped or any(ch not in "ACGT" for ch in stripped):
            raise FASTAError(
                f"line {lineno}: sequence lines may contain A/C/G/T only "
                "(upper case, no gaps or N)"
            )
        parts.append(stripped)

    if name is None:
        raise FASTAError("missing '>' header line")
    sequence = "".join(parts)
    if not sequence:
        raise FASTAError("FASTA contains no sequence")
    return FastaInfo(name=name, sequence=sequence, sha256=digest)


def iter_vcf_records(gz_path: str):
    """Yield (record_index, chrom, pos, ref_allele) for every VCF data line.

    Verification is a sequential whole-file scan, so this uses the standard
    library gzip reader rather than htslib: BGZF is a valid multi-member gzip
    stream, and the stdlib implementation is safe to use concurrently from the
    two app instances' worker threads on different files.
    """
    index = 0
    with gzip.open(gz_path, "rt", encoding="utf-8", newline="") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            index += 1
            yield index, fields[0], int(fields[1]), fields[3]


def verify(fasta_name: str, sequence: str, gz_path: str, vcf_chrom: str) -> VerifyResult:
    """Verify a built VCF against a reference sequence.

    The whole file is walked; the result is an exhaustive (but capped at
    MAX_PROBLEMS) list of concrete problems rather than a single boolean.
    """
    problems: list[Problem] = []
    truncated = False
    checked = 0
    ref_len = len(sequence)

    def add(problem: Problem) -> bool:
        nonlocal truncated
        if len(problems) < MAX_PROBLEMS:
            problems.append(problem)
            return True
        truncated = True
        return False

    for index, chrom, pos, allele in iter_vcf_records(gz_path):
        checked += 1
        if chrom != fasta_name or chrom != vcf_chrom:
            add(
                Problem(
                    record=index,
                    chrom=chrom,
                    position=pos,
                    reason="chromosome_mismatch",
                    detail=(
                        f"VCF chromosome '{chrom}' does not match reference "
                        f"sequence '{fasta_name}'"
                    ),
                )
            )
            # no positional checks are meaningful against another chromosome
            continue
        end = pos + len(allele) - 1
        if pos < 1 or pos > ref_len or end > ref_len:
            add(
                Problem(
                    record=index,
                    chrom=chrom,
                    position=pos,
                    reason="position_out_of_range",
                    detail=(
                        f"REF span [{pos}, {end}] is outside the reference "
                        f"sequence '{fasta_name}' of length {ref_len}"
                    ),
                    ref=allele,
                )
            )
            continue
        expected = sequence[pos - 1:end]
        if expected != allele:
            add(
                Problem(
                    record=index,
                    chrom=chrom,
                    position=pos,
                    reason="ref_mismatch",
                    detail=(
                        f"REF allele at position {pos} does not match the "
                        f"reference fragment"
                    ),
                    ref=allele,
                    expected=expected,
                )
            )

    return VerifyResult(
        problems=[p.as_dict() for p in problems],
        records_checked=checked,
        truncated=truncated,
    )
