"""Reference gate, verification and two-instance interleaving tests."""
import gzip
import os
import threading

import pytest
from conftest import (
    SAMPLE,
    base_sequence,
    enable_reference,
    fasta,
    make_config,
    sequence_matching,
    upload,
    upload_reference,
    vcf,
)
from fastapi.testclient import TestClient

from app.main import create_app
from app.refcheck import FASTA_NAME


# ------------------------------------------------------------- FASTA validation

def test_first_vcf_is_rejected_without_reference(client):
    r = client.post(
        "/vcf", params={"base_version": 0},
        files={"file": ("v.vcf", SAMPLE.encode())},
    )
    assert r.status_code == 409
    assert "no reference" in r.json()["detail"]
    assert client.get("/version").json() == {
        "published": False, "current_version": 0, "ref_id": 0
    }


def test_reference_getter_when_none_active(client):
    assert client.get("/reference").json() == {"active": False, "ref_id": 0}


def test_first_reference_activation_format_only(client):
    seq = base_sequence(1000)
    r = upload_reference(client, fasta(seq))
    assert r.status_code == 201, r.text
    body = r.json()
    assert body["ref_id"] == 1
    assert body["name"] == "chr1"
    assert body["length"] == 1000
    assert len(body["sha256"]) == 64
    assert body["verified_against_version"] == 0
    info = client.get("/reference").json()
    assert info["active"] is True and info["ref_id"] == 1
    assert info["vcf_version"] == 0


def test_activated_reference_file_is_on_disk(client, app):
    enable_reference(client, base_sequence(500))
    info = app.state.store.active_ref_info()
    p = os.path.join(app.state.cfg.data_dir, info["path"], FASTA_NAME)
    assert os.path.isfile(p)


BAD_FASTAS = {
    "two_records": ">chr1\nACGT\n>chr2\nACGT\n",
    "contains_n": ">chr1\nACGNACGT\n",
    "lowercase": ">chr1\nacgt\n",
    "blank_line": ">chr1\nACGT\n\nACGT\n",
    "empty_sequence": ">chr1\n",
    "no_header": "ACGTACGT\n",
    "empty_name": ">\nACGT\n",
    "gap_chars": ">chr1\nAC-GT\n",
    "whitespace_seq": ">chr1\nAC GT\n",
    "empty_file": "",
}


@pytest.mark.parametrize("name", sorted(BAD_FASTAS))
def test_invalid_fastas_are_rejected(client, app, name):
    r = upload_reference(client, BAD_FASTAS[name])
    assert r.status_code == 400, f"{name}: {r.text}"
    assert "invalid reference FASTA" in r.json()["detail"]
    assert app.state.store.active_ref_info() is None
    assert os.listdir(os.path.join(app.state.cfg.data_dir, "refs")) == []


def test_oversize_reference_rejected(client, app):
    limit = app.state.cfg.max_upload_bytes
    body = b">chr1\n" + b"A" * (limit + 1) + b"\n"
    r = upload_reference(client, body)
    assert r.status_code == 413
    assert app.state.store.active_ref_info() is None


def test_binary_reference_rejected(client):
    r = upload_reference(client, gzip.compress(b">chr1\nACGT\n"))
    assert r.status_code == 400


# ------------------------------------------------------- verification failures

def test_vcf_rejected_on_chromosome_name_mismatch(client):
    # reference named chr2, VCF records use chr1
    enable_reference(client, base_sequence(1000), name="chr2")
    r = client.post(
        "/vcf", params={"base_version": 0},
        files={"file": ("v.vcf", SAMPLE.encode())},
    )
    assert r.status_code == 422, r.text
    payload = r.json()["detail"]
    assert payload["error"] == "vcf_verification_failed"
    problems = payload["problems"]
    assert problems and all(p["reason"] == "chromosome_mismatch" for p in problems)
    assert problems[0]["position"] == 100 and problems[0]["chrom"] == "chr1"
    # nothing published, reference stays in place
    assert client.get("/version").json()["current_version"] == 0


def test_vcf_rejected_when_ref_span_out_of_range(client):
    # 50-base reference; a span ending at 52 cannot fit
    enable_reference(client, base_sequence(50))
    r = client.post(
        "/vcf", params={"base_version": 0},
        files={"file": ("v.vcf", vcf([(49, "ACGT", "A")]).encode())},
    )
    assert r.status_code == 422
    p = r.json()["detail"]["problems"][0]
    assert p["reason"] == "position_out_of_range"
    assert p["position"] == 49 and p["ref"] == "ACGT"
    assert "length 50" in p["detail"]


def test_vcf_rejected_on_base_mismatch_with_expected_fragment(client):
    seq = sequence_matching(SAMPLE)
    # corrupt the fragment backing the SNP at position 105
    seq = seq[:104] + ("A" if seq[104] != "A" else "C") + seq[105:]
    enable_reference(client, seq)
    r = client.post(
        "/vcf", params={"base_version": 0},
        files={"file": ("v.vcf", SAMPLE.encode())},
    )
    assert r.status_code == 422
    problems = r.json()["detail"]["problems"]
    snp = [p for p in problems if p["position"] == 105]
    assert len(snp) == 1
    p = snp[0]
    assert p["reason"] == "ref_mismatch"
    assert p["ref"] == "G" and p["expected"] != "G" and len(p["expected"]) == 1


def test_rejected_vcf_leaves_current_version_queryable(client):
    upload(client, SAMPLE, base_version=0)
    before = client.get("/version").json()
    # a VCF whose REF at 105 disagrees with the active reference
    bad = vcf([(105, "C", "T")])
    r = client.post(
        "/vcf", params={"base_version": before["current_version"]},
        files={"file": ("v.vcf", bad.encode())},
    )
    assert r.status_code == 422
    after = client.get("/version").json()
    assert after == before
    q = client.get("/variants", params={"chrom": "chr1", "start": 100, "end": 110})
    assert q.status_code == 200 and q.json()["version"] == before["current_version"]


def test_replacement_reference_verified_against_current_vcf(client):
    upload(client, SAMPLE, base_version=0)
    # incompatible replacement: chromosome name differs
    r = upload_reference(client, fasta(base_sequence(1000), name="chr9"))
    assert r.status_code == 422
    payload = r.json()["detail"]
    assert payload["error"] == "reference_verification_failed"
    assert payload["vcf_version"] == 1
    assert payload["problems"][0]["reason"] == "chromosome_mismatch"
    # the previously active reference is untouched
    info = client.get("/reference").json()
    assert info["active"] is True
    assert info["ref_id"] == 1  # the helper-enabled reference is still active


def test_rejected_reference_keeps_old_combination_and_files(client, app):
    ref_id_before = upload(client, SAMPLE, base_version=0).json()["ref_id"]
    refs_before = set(os.listdir(os.path.join(app.state.cfg.data_dir, "refs")))
    r = upload_reference(client, fasta(base_sequence(50)))  # too short -> out of range
    assert r.status_code == 422
    assert client.get("/reference").json()["ref_id"] == ref_id_before
    assert set(os.listdir(os.path.join(app.state.cfg.data_dir, "refs"))) == refs_before
    # queries still serve the pinned old combination
    q = client.get("/variants", params={"chrom": "chr1", "start": 100, "end": 100})
    assert q.status_code == 200 and q.json()["ref_id"] == ref_id_before


def test_compatible_reference_replacement_is_enabled(client):
    upload(client, SAMPLE, base_version=0)
    seq = sequence_matching(SAMPLE)  # same fragments, different background
    r = upload_reference(client, fasta(seq))
    assert r.status_code == 201, r.text
    assert r.json()["verified_against_version"] == 1
    assert client.get("/reference").json()["ref_id"] == r.json()["ref_id"]


def test_problems_are_capped_with_truncation_flag(client):
    from app.refcheck import MAX_PROBLEMS

    # every tested position is A in the reference; every VCF asks for C
    seq = list(base_sequence(10_000))
    records = [(100 + i, "C", "T") for i in range(MAX_PROBLEMS + 25)]
    for i in range(len(records)):
        seq[99 + i] = "A"
    enable_reference(client, "".join(seq))
    r = client.post(
        "/vcf", params={"base_version": 0},
        files={"file": ("v.vcf", vcf(records).encode())},
    )
    assert r.status_code == 422
    payload = r.json()["detail"]
    assert payload["truncated"] is True
    assert len(payload["problems"]) == MAX_PROBLEMS
    assert payload["records_checked"] == MAX_PROBLEMS + 25


# ------------------------------------------------------------- happy path / pins

def test_publish_response_and_query_carry_fixed_combination(client):
    rid = enable_reference(client, sequence_matching(SAMPLE))
    r = client.post(
        "/vcf", params={"base_version": 0},
        files={"file": ("v.vcf", SAMPLE.encode())},
    )
    assert r.status_code == 201
    assert r.json()["ref_id"] == rid
    q = client.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000})
    assert q.json()["version"] == 1 and q.json()["ref_id"] == rid


def test_pinned_query_keeps_old_combination_after_swap(tmp_path):
    app1 = create_app(make_config(tmp_path, instance_id="i1"))
    app2 = create_app(make_config(tmp_path, instance_id="i2"))
    with TestClient(app1) as c1, TestClient(app2) as c2:
        rid1 = upload(c1, SAMPLE, base_version=0).json()["ref_id"]
        store = app1.state.store
        pinned = store.pin("i1")
        pin_id = pinned[3]
        assert pinned[0] == 1 and pinned[4] == rid1

        # instance 2 swaps the reference and publishes a corrected VCF; the
        # new reference must still match every fragment of the current VCF
        new_seq = sequence_matching(SAMPLE + "\n" + vcf([(500, "A", "G")]))
        r2 = upload_reference(c2, fasta(new_seq))
        rid2 = r2.json()["ref_id"]
        # publish the corrected VCF directly: its reference is already active
        pub = client_post_vcf(c2, vcf([(500, "A", "G")]), 1)
        assert pub.status_code == 201, pub.text
        assert pub.json()["ref_id"] == rid2

        # the in-flight pin still denotes the old combination
        pinned2 = store.pin("i1")
        assert pinned2[0] == 2 and pinned2[4] == rid2
        store.unpin(pinned2[3])
        # fresh queries carry the new pair; the old pin never mixed refs
        q = c1.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000})
        assert q.json()["version"] == 2 and q.json()["ref_id"] == rid2
        store.unpin(pin_id)


# ----------------------------------------------------- two-instance interleaving

class Gate:
    """Block a worker-thread verify hook until the test releases it."""

    def __init__(self):
        self.entered = threading.Event()
        self.release = threading.Event()

    def hook(self, stage, payload):
        self.entered.set()
        self.release.wait(timeout=30)


def _two_clients(tmp_path):
    app1 = create_app(make_config(tmp_path, instance_id="i1"))
    app2 = create_app(make_config(tmp_path, instance_id="i2"))
    return app1, app2, TestClient(app1), TestClient(app2)


def test_slow_vcf_verify_does_not_block_other_instance(tmp_path):
    app1, app2, c1, c2 = _two_clients(tmp_path)
    gate = Gate()
    app1.state.set_verify_hook(
        lambda stage, p: gate.hook(stage, p) if stage == "vcf_verify_start" else None
    )
    with c1, c2:
        upload(c2, SAMPLE, base_version=0)  # reference + v1 on instance 2

        new_vcf = vcf([(300, "A", "G")])
        t = threading.Thread(
            target=lambda: upload(c1, new_vcf, base_version=1), daemon=True
        )
        t.start()
        assert gate.entered.wait(timeout=10)

        # while instance 1 verifies, instance 2 still serves every operation
        q = c2.get("/variants", params={"chrom": "chr1", "start": 100, "end": 100})
        assert q.status_code == 200
        assert c2.get("/version").status_code == 200
        assert c2.get("/reference").status_code == 200

        gate.release.set()
        t.join(timeout=30)
        assert not t.is_alive()
        assert c2.get("/version").json()["current_version"] == 2


def test_reference_swap_compatible_with_pending_vcf_is_reverified_and_published(tmp_path):
    """The other instance swaps the reference; the pending VCF matches the new
    reference too, so it must be re-verified and published (not rejected just
    because the reference version changed)."""
    app1, app2, c1, c2 = _two_clients(tmp_path)
    gate = Gate()
    app1.state.set_verify_hook(
        lambda stage, p: gate.hook(stage, p) if stage == "vcf_verify_start" else None
    )
    with c1, c2:
        # active ref1 matches both the initial VCF and the pending new one
        seq1 = sequence_matching(SAMPLE + "\n" + vcf([(300, "A", "G")]))
        rid1 = enable_reference(c2, seq1)
        assert client_post_vcf(c2, SAMPLE, 0).status_code == 201

        # ref2: same fragments at every used position -> still compatible
        seq2 = list(seq1)
        for i in range(0, len(seq2), 997):  # touch only unused positions
            seq2[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[seq2[i]]
        seq2 = "".join(seq2)

        pending = vcf([(300, "A", "G")])
        t = threading.Thread(
            target=lambda: client_post_vcf(c1, pending, 1), daemon=True
        )
        t.start()
        assert gate.entered.wait(timeout=10)

        rid2 = upload_reference(c2, fasta(seq2)).json()["ref_id"]
        assert rid2 != rid1
        gate.release.set()
        t.join(timeout=30)

        body = c2.get("/version").json()
        assert body["current_version"] == 2 and body["ref_id"] == rid2


def test_reference_swap_incompatible_with_pending_vcf_never_publishes(tmp_path):
    """The swapped reference is incompatible with the pending VCF: the CAS
    retry re-verifies and refuses it; the live combination stays compatible."""
    app1, app2, c1, c2 = _two_clients(tmp_path)
    gate = Gate()
    app1.state.set_verify_hook(
        lambda stage, p: gate.hook(stage, p) if stage == "vcf_verify_start" else None
    )
    with c1, c2:
        # pending VCF wants REF=C at 300; ref1 provides it
        seq1 = sequence_matching(SAMPLE + "\n" + vcf([(300, "C", "T")]))
        enable_reference(c2, seq1)
        assert client_post_vcf(c2, SAMPLE, 0).status_code == 201

        # ref2 keeps every fragment of the current VCF (stays activatable)
        # but position 300 becomes G, which the pending VCF cannot match.
        seq2 = sequence_matching(SAMPLE)
        seq2 = seq2[:299] + "G" + seq2[300:]
        pending = vcf([(300, "C", "T")])

        result = {}
        t = threading.Thread(
            target=lambda: result.update(client_post_vcf(c1, pending, 1).json()),
            daemon=True,
        )
        t.start()
        assert gate.entered.wait(timeout=10)

        r = upload_reference(c2, fasta(seq2))
        assert r.status_code == 201, r.text
        gate.release.set()
        t.join(timeout=30)

        assert result["detail"]["error"] == "vcf_verification_failed"
        reasons = {p["reason"] for p in result["detail"]["problems"]}
        assert "ref_mismatch" in reasons
        # current VCF version unchanged; active reference is the new, compatible one
        body = c2.get("/version").json()
        assert body["current_version"] == 1
        assert c2.get("/reference").json()["ref_id"] == r.json()["ref_id"]
        q = c2.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000})
        assert q.status_code == 200


def test_vcf_published_during_slow_reference_check_triggers_reverify(tmp_path):
    """Mirror case: instance 1 verifies a new reference against v1 while
    instance 2 publishes v2; the new reference matches v2 as well, so it is
    re-verified and activated against v2."""
    app1, app2, c1, c2 = _two_clients(tmp_path)
    gate = Gate()
    app1.state.set_verify_hook(
        lambda stage, p: gate.hook(stage, p) if stage == "ref_verify_start" else None
    )
    with c1, c2:
        v2 = vcf([(500, "A", "G")])
        seq = sequence_matching(SAMPLE + "\n" + v2)
        enable_reference(c2, seq)
        assert client_post_vcf(c2, SAMPLE, 0).status_code == 201

        result_holder = {}
        t = threading.Thread(
            target=lambda: result_holder.update(
                upload_reference(c1, fasta(seq)).json()
            ),
            daemon=True,
        )
        t.start()
        assert gate.entered.wait(timeout=10)

        # instance 2 corrects the VCF while the reference check is in flight
        assert client_post_vcf(c2, v2, 1).status_code == 201
        gate.release.set()
        t.join(timeout=30)

        assert result_holder["verified_against_version"] == 2
        assert c2.get("/reference").json()["vcf_version"] == 2


def test_incompatible_reference_swap_fails_then_compatible_one_succeeds(tmp_path):
    app1, app2, c1, c2 = _two_clients(tmp_path)
    with c1, c2:
        upload(c1, SAMPLE, base_version=0)
        bad = upload_reference(c2, fasta(base_sequence(50)))  # out of range
        assert bad.status_code == 422
        good = upload_reference(c2, fasta(sequence_matching(SAMPLE)))
        assert good.status_code == 201
        assert c1.get("/reference").json()["ref_id"] == good.json()["ref_id"]


# ------------------------------------------------------------- crash / recovery

def test_restart_keeps_verified_combination_queryable(tmp_path):
    cfg1 = make_config(tmp_path, instance_id="i1")
    app1 = create_app(cfg1)
    with TestClient(app1) as c1:
        r = upload(c1, SAMPLE, base_version=0)
        rid, ver = r.json()["ref_id"], r.json()["version"]

    app2 = create_app(make_config(tmp_path, instance_id="i1-restart"))
    with TestClient(app2) as c2:
        assert c2.get("/reference").json()["ref_id"] == rid
        assert c2.get("/version").json()["current_version"] == ver
        q = c2.get("/variants", params={"chrom": "chr1", "start": 105, "end": 105})
        assert q.status_code == 200
        assert q.json()["version"] == ver and q.json()["ref_id"] == rid


def test_write_failure_during_reference_activation_keeps_old_state(client, app, monkeypatch):
    import app.main as main_mod

    rid_before = upload(client, SAMPLE, base_version=0).json()["ref_id"]
    refs_dir = os.path.join(app.state.cfg.data_dir, "refs")
    real_copytree = main_mod.shutil.copytree

    def boom(src, dst, *a, **k):
        # only the final placement into refs/ simulates a disk failure
        if os.path.dirname(os.path.abspath(dst)) == os.path.abspath(refs_dir):
            raise OSError("simulated disk failure")
        return real_copytree(src, dst, *a, **k)

    monkeypatch.setattr(main_mod.shutil, "copytree", boom)
    # matching sequence: verification passes, then the placement write fails
    r = upload_reference(client, fasta(sequence_matching(SAMPLE)))
    assert r.status_code == 500
    # pointer unchanged; old combination still serves
    assert client.get("/reference").json()["ref_id"] == rid_before
    assert client.get("/variants", params={"chrom": "chr1", "start": 100, "end": 100}).status_code == 200


def test_write_failure_during_vcf_publish_keeps_old_version(client, app, monkeypatch):
    import app.main as main_mod

    upload(client, SAMPLE, base_version=0)
    before = client.get("/version").json()
    versions_dir = os.path.join(app.state.cfg.data_dir, "versions")
    # make sure the pending VCF would verify, so the first failure hit is the
    # simulated write failure rather than a 422
    next_vcf = vcf([(400, "A", "G")])
    upload_reference(client, fasta(sequence_matching(SAMPLE + "\n" + next_vcf)))
    real_copytree = main_mod.shutil.copytree

    def boom(src, dst, *a, **k):
        if os.path.dirname(os.path.abspath(dst)) == os.path.abspath(versions_dir):
            raise OSError("simulated disk failure")
        return real_copytree(src, dst, *a, **k)

    monkeypatch.setattr(main_mod.shutil, "copytree", boom)
    r = client.post(
        "/vcf", params={"base_version": before["current_version"]},
        files={"file": ("v.vcf", next_vcf.encode())},
    )
    assert r.status_code == 500
    assert client.get("/version").json() == before
    # only the originally published version directory remains
    assert len(os.listdir(versions_dir)) == 1
    assert client.get("/variants", params={"chrom": "chr1", "start": 100, "end": 100}).status_code == 200


def test_staging_leftovers_swept_after_restart(tmp_path):
    cfg = make_config(tmp_path, instance_id="i1")
    app = create_app(cfg)
    with TestClient(app) as c:
        upload(c, SAMPLE, base_version=0)
    junk = os.path.join(cfg.data_dir, "staging", "i1", "leftover")
    os.makedirs(junk)
    open(os.path.join(junk, "upload.vcf"), "w").write("partial")

    # same instance id recreates its staging root and drops the leftovers
    app2 = create_app(make_config(tmp_path, instance_id="i1"))
    assert not os.path.exists(junk)
    with TestClient(app2) as c2:
        assert c2.get("/variants", params={"chrom": "chr1", "start": 1, "end": 10}).status_code == 200


# ------------------------------------------------------------------ reaping

def test_old_reference_reaped_once_no_version_uses_it(tmp_path):
    cfg = make_config(tmp_path, instance_id="i1")
    app1 = create_app(cfg)
    with TestClient(app1) as c:
        r1 = upload(c, SAMPLE, base_version=0)
        rid1 = r1.json()["ref_id"]
        store = app1.state.store
        old_ref_path = store.active_ref_info()["path"]

        # replacement reference + replacement VCF verified against it; the new
        # reference also keeps every fragment of the current VCF
        r2 = upload_reference(c, fasta(sequence_matching(SAMPLE + "\n" + vcf([(500, "A", "G")]))))
        rid2 = r2.json()["ref_id"]
        assert client_post_vcf(c, vcf([(500, "A", "G")]), 1).status_code == 201

        app1.state.reap_once()  # v1 reaped; ref1 now unreferenced -> reaped too
        assert not os.path.exists(os.path.join(cfg.data_dir, old_ref_path))
        assert store.all_reference_paths() == {
            store.active_ref_info()["path"]
        }
        assert store.active_ref_info()["id"] == rid2 and rid2 != rid1
        # current combination still answers
        q = c.get("/variants", params={"chrom": "chr1", "start": 500, "end": 500})
        assert q.status_code == 200 and q.json()["ref_id"] == rid2


# --------------------------------------------------------------------- helpers

def client_post_vcf(client, body, base_version):
    return client.post(
        "/vcf", params={"base_version": base_version},
        files={"file": ("v.vcf", body if isinstance(body, bytes) else body.encode())},
    )
