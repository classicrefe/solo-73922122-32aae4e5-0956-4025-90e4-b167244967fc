"""Shared fixtures/helpers for the test suite."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

from app.main import Config, create_app

META = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
CHROM_LINE = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
HEADER = META + CHROM_LINE

# Reference sequence used by the helpers: chr1, 200 kb of A/C/G/T, fixed seed.
REF_LENGTH = 200_000


def vcf(records, header=HEADER):
    """Build a VCF body from (pos, ref, alt) tuples or raw strings."""
    lines = []
    for r in records:
        if isinstance(r, str):
            lines.append(r)
        else:
            pos, ref, alt = r
            lines.append(f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.")
    return header + "\n".join(lines) + "\n"


# chr1:100 deletion spanning 100-103, SNP at 105, insertion at 110, del 200-204
SAMPLE = vcf([
    (100, "ACGT", "A"),
    (105, "G", "T"),
    (110, "C", "CT"),
    (200, "TAAAA", "T"),
])


def base_sequence(length=REF_LENGTH, seed=12345):
    """Deterministic ACGT background sequence."""
    bases = "ACGT"
    out = []
    state = seed
    for _ in range(length):
        # xorshift32: deterministic, dependency-free pseudo-random bases.
        state ^= (state << 13) & 0xFFFFFFFF
        state ^= state >> 17
        state ^= (state << 5) & 0xFFFFFFFF
        state &= 0xFFFFFFFF
        out.append(bases[state % 4])
    return "".join(out)


def sequence_matching(body: str, length=REF_LENGTH):
    """Build a chr1 sequence whose fragments match every parseable VCF REF."""
    seq = list(base_sequence(length))
    for line in body.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        try:
            pos = int(fields[1])
        except ValueError:
            continue
        ref = fields[3]
        if pos < 1 or pos + len(ref) - 1 > length:
            continue
        if any(ch not in "ACGT" for ch in ref):
            continue
        for i, ch in enumerate(ref):
            seq[pos - 1 + i] = ch
    return "".join(seq)


def fasta(sequence, name="chr1"):
    """Wrap a sequence into a single-record FASTA body (80-col lines)."""
    body = f">{name}\n"
    for i in range(0, len(sequence), 80):
        body += sequence[i:i + 80] + "\n"
    return body


def make_config(tmp_path, **kw):
    kw.setdefault("instance_id", "test-inst")
    kw.setdefault("reaper_interval", 0)    # deterministic: tests call reap_once()
    kw.setdefault("heartbeat_interval", 0)
    return Config(data_dir=str(tmp_path / "data"), **kw)


@pytest.fixture()
def app(tmp_path):
    return create_app(make_config(tmp_path))


@pytest.fixture()
def client(app):
    from fastapi.testclient import TestClient

    with TestClient(app) as c:
        yield c


# data_dir -> VCF bodies seen by the helpers; keeps reference swaps consistent
# when two test instances share one data directory.
_SEEN_BODIES: dict[str, list[str]] = {}


def upload(client, body, base_version):
    # The service gates VCF publishing on an enabled reference and verifies
    # every replacement VCF in full. The helper therefore keeps a reference
    # whose fragments match the union of every VCF uploaded against this data
    # directory, swapping it in before each upload so the tests can focus on
    # their original validation/concurrency/query behaviour.
    text = body.decode("utf-8", "ignore") if isinstance(body, bytes) else body
    data_dir = client.app.state.cfg.data_dir
    seen = _SEEN_BODIES.setdefault(data_dir, [])
    union = sequence_matching("\n".join(seen + [text]))
    if client.get("/reference").json().get("active") is False:
        r = upload_reference(client, fasta(union))
    else:
        r = upload_reference(client, fasta(union))
    assert r.status_code == 201, r.text
    seen.append(text)
    return client.post(
        "/vcf",
        params={"base_version": base_version},
        files={"file": ("variants.vcf", body if isinstance(body, bytes) else body.encode())},
    )


def upload_reference(client, body, filename="reference.fasta"):
    return client.post(
        "/reference",
        files={"file": (filename, body if isinstance(body, bytes) else body.encode())},
    )


def enable_reference(client, sequence=None, name="chr1"):
    """Enable a reference (matching SAMPLE by default) and return its ref_id."""
    r = upload_reference(client, fasta(sequence or sequence_matching(SAMPLE), name))
    assert r.status_code == 201, r.text
    return r.json()["ref_id"]
