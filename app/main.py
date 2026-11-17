"""FastAPI application: reference-gated, versioned VCF upload + index queries.

Layout of the shared data directory:

    <data_dir>/meta.sqlite3              metadata, pointers, pins
    <data_dir>/versions/<uuid>/          published VCF sets (variants.vcf.gz + .tbi)
    <data_dir>/refs/<uuid>/              activated references (reference.fasta)
    <data_dir>/staging/<instance>/...    in-progress builds (both kinds)

Gate rules:

* no VCF may be published until a reference has been enabled (``POST
  /reference``); the first reference is validated for format only because no
  VCF exists yet;
* every later reference upload is verified against the *current* VCF —
  chromosome name and every REF fragment — and is refused (the active
  reference stays untouched) on any mismatch or out-of-range span;
* every later VCF upload is verified in full against the *active* reference
  before publishing;
* both pointer flips are compare-and-swap on the (reference id, VCF version)
  pair. If the other instance changes one input while a slow verification is
  running, the combination is re-verified instead of blindly rejected or
  published; a genuinely incompatible combination never goes live.

Verification runs in worker threads, so a slow check never blocks the other
instance's uploads, reference swaps or queries, and no global serialization is
held during it. Publishing itself stays atomic: completed file sets are
renamed into place and only then is the pointer flipped inside one SQLite
transaction.
"""
from __future__ import annotations

import asyncio
import os
import shutil
import sqlite3
import time
import uuid
from contextlib import asynccontextmanager
from dataclasses import dataclass

import pysam
from fastapi import FastAPI, File, HTTPException, Query, UploadFile
from fastapi.responses import RedirectResponse

from .db import StateChanged, Store
from .refcheck import FASTA_NAME, FASTAError, parse_fasta, verify
from .vcfbuild import GZ_NAME, VCFError, validate_and_build

DEFAULT_MAX_UPLOAD = 10 * 1024 * 1024  # 10 MiB
# Bounded CAS retries: each retry means the other instance changed the state
# underneath an in-flight verification, which requires a full re-check.
MAX_STATE_RETRIES = 10


@dataclass
class Config:
    data_dir: str
    instance_id: str
    max_upload_bytes: int = DEFAULT_MAX_UPLOAD
    reaper_interval: float = 2.0      # seconds between reaper runs (0 disables)
    heartbeat_interval: float = 2.0   # seconds between heartbeats (0 disables)
    dead_after: float = 15.0          # instance considered dead after this silence
    orphan_grace: float = 60.0        # min age before orphan dirs are swept


def _config_from_env() -> Config:
    e = os.environ
    return Config(
        data_dir=e.get("VCF_DATA_DIR", "/data"),
        instance_id=e.get("VCF_INSTANCE_ID") or uuid.uuid4().hex,
        max_upload_bytes=int(e.get("VCF_MAX_UPLOAD_BYTES", DEFAULT_MAX_UPLOAD)),
        reaper_interval=float(e.get("VCF_REAPER_INTERVAL", "2")),
        heartbeat_interval=float(e.get("VCF_HEARTBEAT_INTERVAL", "2")),
        dead_after=float(e.get("VCF_DEAD_AFTER", "15")),
        orphan_grace=float(e.get("VCF_ORPHAN_GRACE", "60")),
    )


def _fsync_dir(path: str) -> None:
    try:
        fd = os.open(path, os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
    except OSError:
        pass


def create_app(config: Config | None = None) -> FastAPI:
    cfg = config or _config_from_env()
    store = Store(os.path.join(cfg.data_dir, "meta.sqlite3"))
    versions_dir = os.path.join(cfg.data_dir, "versions")
    refs_dir = os.path.join(cfg.data_dir, "refs")
    staging_base = os.path.join(cfg.data_dir, "staging")
    staging_root = os.path.join(staging_base, cfg.instance_id)
    for d in (versions_dir, refs_dir):
        os.makedirs(d, exist_ok=True)
    # Leftover builds of a previous incarnation of this instance can never be
    # published any more; drop them before serving.
    shutil.rmtree(staging_root, ignore_errors=True)
    os.makedirs(staging_root, exist_ok=True)

    # Test/observation hooks: invoked from worker threads around slow steps.
    # Production leaves them at None.
    app_state = {"verify_hook": None}

    def _hook(stage: str, payload: dict) -> None:
        hook = app_state["verify_hook"]
        if hook is not None:
            hook(stage, payload)

    def reap_once() -> None:
        now = time.time()
        # 1) DB: drop dead instances' pins, delete reclaimable rows.
        reaped = store.reap_metadata(now, cfg.dead_after)
        for _version, rel in reaped["versions"]:
            shutil.rmtree(os.path.join(cfg.data_dir, rel), ignore_errors=True)
        for _ref_id, rel in reaped["references"]:
            shutil.rmtree(os.path.join(cfg.data_dir, rel), ignore_errors=True)
        # 2) FS: sweep dirs with no DB row (crashed publish/reaper).
        cutoff = now - cfg.orphan_grace
        known_v = store.all_version_paths()
        known_r = store.all_reference_paths()
        for base, known in ((versions_dir, known_v), (refs_dir, known_r)):
            for name in os.listdir(base):
                path = os.path.join(base, name)
                rel = os.path.relpath(path, cfg.data_dir)
                try:
                    old = os.path.getmtime(path) < cutoff
                except OSError:
                    continue
                if rel not in known and old:
                    shutil.rmtree(path, ignore_errors=True)
        # 3) FS: sweep staging dirs of dead instances.
        live = store.live_instances(now, cfg.dead_after)
        for inst in os.listdir(staging_base):
            path = os.path.join(staging_base, inst)
            try:
                old = os.path.getmtime(path) < cutoff
            except OSError:
                continue
            if inst not in live and old:
                shutil.rmtree(path, ignore_errors=True)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        store.heartbeat(cfg.instance_id, time.time())
        stopped = asyncio.Event()

        async def heartbeat_loop():
            while not stopped.is_set():
                await asyncio.sleep(cfg.heartbeat_interval)
                try:
                    store.heartbeat(cfg.instance_id, time.time())
                except Exception:
                    pass

        async def reaper_loop():
            while not stopped.is_set():
                await asyncio.sleep(cfg.reaper_interval)
                try:
                    await asyncio.to_thread(reap_once)
                except Exception:
                    pass

        tasks = []
        if cfg.heartbeat_interval > 0:
            tasks.append(asyncio.create_task(heartbeat_loop()))
        if cfg.reaper_interval > 0:
            tasks.append(asyncio.create_task(reaper_loop()))
        try:
            yield
        finally:
            stopped.set()
            for t in tasks:
                t.cancel()
            await asyncio.gather(*tasks, return_exceptions=True)
            try:
                store.unregister(cfg.instance_id)
            except Exception:
                pass

    app = FastAPI(
        title="Reference-gated versioned VCF variant store",
        version="2.0.0",
        lifespan=lifespan,
        description=(
            "Upload single-chromosome biallelic VCFs (max 10 MiB, no sample "
            "columns, REF/ALT in ACGT); records are sorted, BGZF-compressed and "
            "Tabix-indexed. A reference FASTA (single sequence, ACGT only, max "
            "10 MiB) must be enabled first: every VCF is verified in full "
            "against the active reference (chromosome name and each REF "
            "fragment) before it can be published, and a replacement reference "
            "is verified against the current VCF before it can be enabled. "
            "Query 1-based closed intervals via the Tabix index; every response "
            "carries the fixed (reference id, VCF version) combination it read."
        ),
    )
    app.state.cfg = cfg
    app.state.store = store
    app.state.reap_once = reap_once
    app.state.verify_hook = None

    def _set_hook(value):  # pragma: no cover - trivial setter used by tests
        app_state["verify_hook"] = value

    app.state.set_verify_hook = _set_hook

    async def _spool_upload(file: UploadFile, staging: str, filename: str) -> str:
        """Stream an upload to staging, enforcing the 10 MiB cap. Returns path."""
        path = os.path.join(staging, filename)
        size = 0
        with open(path, "wb") as out:
            while True:
                chunk = await file.read(1 << 20)
                if not chunk:
                    break
                size += len(chunk)
                if size > cfg.max_upload_bytes:
                    raise HTTPException(
                        413, f"upload exceeds the {cfg.max_upload_bytes}-byte limit"
                    )
                out.write(chunk)
        return path

    @app.get("/", include_in_schema=False)
    def root():
        return RedirectResponse("/docs")

    # ------------------------------------------------------------- references

    @app.post("/reference", status_code=201)
    async def upload_reference(file: UploadFile = File(..., description="Single-record ACGT-only FASTA, max 10 MiB")):
        """Validate and enable a reference FASTA.

        With a VCF already published, the FASTA is verified against the whole
        current VCF; any chromosome mismatch, out-of-range span or REF
        mismatch rejects the upload with the offending positions, leaving the
        active reference and every VCF version untouched. With no VCF yet,
        format validation is enough (first-time activation).
        """
        staging = os.path.join(staging_root, "ref-" + uuid.uuid4().hex)
        os.makedirs(staging)
        try:
            upload_path = await _spool_upload(file, staging, "upload.fasta")
            try:
                info = await asyncio.to_thread(parse_fasta, upload_path)
            except FASTAError as exc:
                raise HTTPException(400, f"invalid reference FASTA: {exc}")
            # Give the immutable file set its canonical name before moving it.
            os.rename(upload_path, os.path.join(staging, FASTA_NAME))

            # CAS loop: re-verify whenever the other instance changes the VCF
            # underneath this verification.
            attempts = 0
            while True:
                attempts += 1
                snap = await asyncio.to_thread(store.state_snapshot)
                vcf = snap["vcf"]
                if vcf is not None:
                    gz_path = os.path.join(cfg.data_dir, vcf["path"], GZ_NAME)
                    await asyncio.to_thread(
                        _hook, "ref_verify_start",
                        {"instance": cfg.instance_id, "vcf_version": vcf["version"]},
                    )
                    result = await asyncio.to_thread(
                        verify, info.name, info.sequence, gz_path, vcf["chrom"]
                    )
                    await asyncio.to_thread(
                        _hook, "ref_verify_done",
                        {"instance": cfg.instance_id, "problems": len(result.problems)},
                    )
                    if result.problems:
                        raise HTTPException(
                            422,
                            {
                                "error": "reference_verification_failed",
                                "message": (
                                    "reference does not match the current VCF; "
                                    "the active reference was left unchanged"
                                ),
                                "vcf_version": vcf["version"],
                                "records_checked": result.records_checked,
                                "truncated": result.truncated,
                                "problems": result.problems,
                            },
                        )
                # Copy the validated file set into place; staging is kept so a
                # lost compare-and-swap can retry. The pointer flip is the last
                # step; a lost flip only leaves an orphan for the reaper.
                final_rel = os.path.join("refs", uuid.uuid4().hex)
                final_path = os.path.join(cfg.data_dir, final_rel)
                try:
                    shutil.copytree(staging, final_path)
                    _fsync_dir(refs_dir)
                except OSError as exc:
                    shutil.rmtree(final_path, ignore_errors=True)
                    raise HTTPException(500, f"failed to persist reference: {exc}")
                try:
                    ref_id = await asyncio.to_thread(
                        store.activate_reference,
                        info.name, info.sha256, info.length, final_rel,
                        snap["current_version"], cfg.instance_id,
                        snap["current_version"], time.time(),
                    )
                except sqlite3.Error as exc:
                    shutil.rmtree(final_path, ignore_errors=True)
                    raise HTTPException(500, f"failed to activate reference: {exc}")
                except StateChanged:
                    shutil.rmtree(final_path, ignore_errors=True)
                    if attempts >= MAX_STATE_RETRIES:
                        raise HTTPException(
                            503,
                            "state kept changing during verification; retry the upload",
                        )
                    continue
                return {
                    "ref_id": ref_id,
                    "name": info.name,
                    "length": info.length,
                    "sha256": info.sha256,
                    "verified_against_version": snap["current_version"],
                }
        finally:
            if os.path.isdir(staging):
                shutil.rmtree(staging, ignore_errors=True)

    @app.get("/reference")
    def get_reference():
        info = store.active_ref_info()
        if info is None:
            return {"active": False, "ref_id": 0}
        return {
            "active": True,
            "ref_id": info["id"],
            "name": info["name"],
            "length": info["length"],
            "sha256": info["sha256"],
            "vcf_version": info["vcf_version"],
            "activated_by": info["activated_by"],
            "created_at": info["created_at"],
        }

    # ------------------------------------------------------------------ VCFs

    @app.post("/vcf", status_code=201)
    async def upload_vcf(
        file: UploadFile = File(...),
        base_version: int = Query(..., ge=0, description="Version the client last saw (0 when empty)"),
    ):
        # Gate: a reference must be enabled before the first VCF is accepted.
        snap0 = await asyncio.to_thread(store.state_snapshot)
        if snap0["active_ref_id"] == 0:
            raise HTTPException(
                409,
                "no reference has been enabled yet; upload and verify a "
                "reference FASTA via POST /reference before publishing a VCF",
            )

        staging = os.path.join(staging_root, uuid.uuid4().hex)
        os.makedirs(staging)
        try:
            upload_path = await _spool_upload(file, staging, "upload.vcf")
            try:
                chrom, n_records = await asyncio.to_thread(
                    validate_and_build, upload_path, staging
                )
            except VCFError as exc:
                raise HTTPException(400, f"invalid VCF: {exc}")
            os.remove(upload_path)
            gz_staging = os.path.join(staging, GZ_NAME)

            # Whole-file verification against the active reference, with a CAS
            # loop so a swapped reference triggers a re-check of the new
            # combination rather than a blind accept/reject.
            attempts = 0
            loaded_ref_id = -1
            sequence = ""
            while True:
                attempts += 1
                snap = await asyncio.to_thread(store.state_snapshot)
                cur = snap["current_version"]
                active_ref_id = snap["active_ref_id"]
                if cur != base_version:
                    # Preserve the original optimistic-concurrency behaviour.
                    raise HTTPException(
                        409,
                        f"base_version {base_version} does not match current version "
                        f"{cur}; re-read and retry",
                    )
                ref_row = snap["ref"]
                if ref_row is None:  # pragma: no cover - guarded above
                    raise HTTPException(409, "no active reference")
                if active_ref_id != loaded_ref_id:
                    sequence = await asyncio.to_thread(
                        _load_sequence,
                        os.path.join(cfg.data_dir, ref_row["path"], FASTA_NAME),
                    )
                    loaded_ref_id = active_ref_id
                await asyncio.to_thread(
                    _hook, "vcf_verify_start",
                    {"instance": cfg.instance_id, "ref_id": active_ref_id},
                )
                result = await asyncio.to_thread(
                    verify, ref_row["name"], sequence, gz_staging, chrom,
                )
                await asyncio.to_thread(
                    _hook, "vcf_verify_done",
                    {"instance": cfg.instance_id, "problems": len(result.problems)},
                )
                if result.problems:
                    raise HTTPException(
                        422,
                        {
                            "error": "vcf_verification_failed",
                            "message": (
                                "VCF does not match the active reference; it "
                                "was not published and the current version is unchanged"
                            ),
                            "ref_id": active_ref_id,
                            "records_checked": result.records_checked,
                            "truncated": result.truncated,
                            "problems": result.problems,
                        },
                    )

                # Copy the verified staging set into place; staging itself is
                # never moved so a lost compare-and-swap can re-verify against
                # the intact files and retry. The DB flip is the last step; a
                # lost flip only leaves an orphaned directory for the reaper.
                final_rel = os.path.join("versions", uuid.uuid4().hex)
                final_path = os.path.join(cfg.data_dir, final_rel)
                try:
                    shutil.copytree(staging, final_path)
                    _fsync_dir(versions_dir)
                except OSError as exc:
                    shutil.rmtree(final_path, ignore_errors=True)
                    raise HTTPException(500, f"failed to persist VCF set: {exc}")
                try:
                    new_version = await asyncio.to_thread(
                        store.publish_vcf,
                        base_version, active_ref_id, chrom, final_rel,
                        n_records, time.time(),
                    )
                except sqlite3.Error as exc:
                    shutil.rmtree(final_path, ignore_errors=True)
                    raise HTTPException(500, f"failed to publish VCF: {exc}")
                except StateChanged as changed:
                    shutil.rmtree(final_path, ignore_errors=True)
                    if attempts >= MAX_STATE_RETRIES:
                        raise HTTPException(
                            503,
                            "state kept changing during verification; retry the upload",
                        )
                    # Either the reference or the VCF pointer changed. Re-check
                    # the whole latest combination; the VCF's optimistic
                    # version guard still applies on the next iteration.
                    base_version = changed.current_version
                    continue
                return {
                    "version": new_version,
                    "ref_id": active_ref_id,
                    "chrom": chrom,
                    "n_records": n_records,
                }
        finally:
            if os.path.isdir(staging):
                shutil.rmtree(staging, ignore_errors=True)

    @app.get("/variants")
    def get_variants(
        chrom: str = Query(..., description="Chromosome name as in the VCF"),
        start: int = Query(..., ge=1, description="1-based interval start (inclusive)"),
        end: int = Query(..., ge=1, description="1-based interval end (inclusive)"),
    ):
        if end < start:
            raise HTTPException(400, "end must be >= start")
        pinned = store.pin(cfg.instance_id, time.time())
        if pinned is None:
            raise HTTPException(404, "no VCF has been published yet")
        # The (version, ref_id) pair is fixed for this query inside the pin
        # transaction; a later swap or publish cannot mix into the response.
        version, vchrom, rel, pin_id, ref_id = pinned
        try:
            if chrom != vchrom:
                raise HTTPException(
                    400, f"the published dataset contains chromosome '{vchrom}' only"
                )
            gz_path = os.path.join(cfg.data_dir, rel, GZ_NAME)
            try:
                tb = pysam.TabixFile(gz_path)
                try:
                    # 1-based closed [start, end] -> 0-based half-open fetch.
                    lines = tb.fetch(vchrom, start - 1, end)
                    records = []
                    for ln in lines:
                        f = ln.split("\t", 4)
                        pos = int(f[1])
                        # Exact span check on top of the index iterator.
                        if pos <= end and pos + len(f[3]) - 1 >= start:
                            records.append(ln)
                finally:
                    tb.close()
            except HTTPException:
                raise
            except Exception as exc:
                raise HTTPException(500, f"indexed fetch failed: {exc}")
            return {
                "version": version,
                "ref_id": ref_id,
                "chrom": vchrom,
                "start": start,
                "end": end,
                "count": len(records),
                "records": records,
            }
        finally:
            store.unpin(pin_id)

    @app.get("/version")
    def get_version():
        info = store.version_info()
        if info is None:
            return {"published": False, "current_version": 0, "ref_id": 0}
        return {
            "published": True,
            "current_version": info["version"],
            "ref_id": info["ref_id"],
            "chrom": info["chrom"],
            "n_records": info["n_records"],
            "created_at": info["created_at"],
        }

    return app


def _load_sequence(fasta_path: str) -> str:
    """Read the sequence of an activated, already-validated reference."""
    return parse_fasta(fasta_path).sequence


# Run with: uvicorn --factory app.main:create_app
