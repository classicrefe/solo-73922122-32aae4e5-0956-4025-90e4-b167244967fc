"""SQLite metadata layer for the versioned VCF + reference store.

All multi-step state changes happen inside BEGIN IMMEDIATE transactions so
that two app instances sharing the same database file are fully serialized.
The "current combination" is a pair:

    (active reference id, current VCF version)

Both reference activation and VCF publishing are compare-and-swap operations
on this pair: the expensive file validation runs without any lock, and only
the final pointer flip takes the write lock, re-checking that the inputs the
validation was based on are still current. A change of either input therefore
forces a re-validation of the latest combination instead of being blindly
rejected or, worse, published as a mismatched pair.
"""
from __future__ import annotations

import os
import sqlite3
import time
from contextlib import contextmanager

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    k TEXT PRIMARY KEY,
    v TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS versions (
    version   INTEGER PRIMARY KEY,
    ref_id    INTEGER NOT NULL,
    chrom     TEXT NOT NULL,
    path      TEXT NOT NULL,      -- directory with variants.vcf.gz + .tbi, relative to data dir
    n_records INTEGER NOT NULL,
    created_at REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS reference_sets (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    name       TEXT NOT NULL,     -- FASTA sequence/chromosome name
    sha256     TEXT NOT NULL,
    length     INTEGER NOT NULL,
    path       TEXT NOT NULL,     -- directory with reference.fasta, relative to data dir
    vcf_version INTEGER NOT NULL, -- VCF version this reference was verified against (0 = none)
    activated_by TEXT NOT NULL,
    created_at REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS refs (
    id        INTEGER PRIMARY KEY AUTOINCREMENT,
    version   INTEGER NOT NULL,
    instance  TEXT NOT NULL,
    created_at REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS instances (
    id        TEXT PRIMARY KEY,
    last_seen REAL NOT NULL
);
"""


class Conflict(Exception):
    """Raised when base_version does not match the current version."""

    def __init__(self, current: int):
        super().__init__(f"current version is {current}")
        self.current = current


class StateChanged(Exception):
    """The (reference, VCF) combination changed while validation was running.

    The caller must re-validate against the state carried on this exception
    and retry the pointer flip.
    """

    def __init__(self, current_version: int, active_ref_id: int):
        super().__init__(
            f"state changed: current VCF version {current_version}, "
            f"active reference {active_ref_id}"
        )
        self.current_version = current_version
        self.active_ref_id = active_ref_id


class Store:
    def __init__(self, db_path: str):
        self.db_path = db_path
        parent = os.path.dirname(db_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        con = sqlite3.connect(db_path, timeout=30.0)
        try:
            # WAL: readers never block the single writer and vice versa.
            con.execute("PRAGMA journal_mode=WAL")
            con.execute("PRAGMA synchronous=NORMAL")
            con.executescript(SCHEMA)
            self._migrate(con)
            con.execute("INSERT OR IGNORE INTO meta(k, v) VALUES ('current_version', '0')")
            con.execute("INSERT OR IGNORE INTO meta(k, v) VALUES ('active_ref_id', '0')")
            con.commit()
        finally:
            con.close()

    @staticmethod
    def _migrate(con: sqlite3.Connection) -> None:
        """Upgrade databases created before the reference feature."""
        cols = {r[1] for r in con.execute("PRAGMA table_info(versions)")}
        if "ref_id" not in cols:
            con.execute("ALTER TABLE versions ADD COLUMN ref_id INTEGER NOT NULL DEFAULT 0")

    def _connect(self) -> sqlite3.Connection:
        con = sqlite3.connect(self.db_path, timeout=30.0)
        con.execute("PRAGMA busy_timeout=30000")
        con.row_factory = sqlite3.Row
        return con

    @contextmanager
    def _txn(self):
        """BEGIN IMMEDIATE transaction: serialized against all other writers."""
        con = self._connect()
        try:
            con.execute("BEGIN IMMEDIATE")
            yield con
            con.commit()
        except BaseException:
            con.rollback()
            raise
        finally:
            con.close()

    # ------------------------------------------------------------------ reads

    def _state(self, con: sqlite3.Connection) -> tuple[int, int]:
        cur = int(
            con.execute("SELECT v FROM meta WHERE k = 'current_version'").fetchone()["v"]
        )
        ref = int(
            con.execute("SELECT v FROM meta WHERE k = 'active_ref_id'").fetchone()["v"]
        )
        return cur, ref

    def state_snapshot(self) -> dict:
        """Current (reference, VCF) combination plus its description rows."""
        con = self._connect()
        try:
            cur, ref = self._state(con)
            vrow = con.execute(
                "SELECT version, ref_id, chrom, path, n_records, created_at "
                "FROM versions WHERE version = ?",
                (cur,),
            ).fetchone()
            rrow = con.execute(
                "SELECT id, name, sha256, length, path, vcf_version, activated_by, created_at "
                "FROM reference_sets WHERE id = ?",
                (ref,),
            ).fetchone()
            return {
                "current_version": cur,
                "active_ref_id": ref,
                "vcf": dict(vrow) if vrow is not None else None,
                "ref": dict(rrow) if rrow is not None else None,
            }
        finally:
            con.close()

    def current_version(self) -> int:
        con = self._connect()
        try:
            return int(
                con.execute("SELECT v FROM meta WHERE k = 'current_version'").fetchone()["v"]
            )
        finally:
            con.close()

    def version_info(self):
        """Row of the currently published version, or None if nothing published."""
        con = self._connect()
        try:
            return con.execute(
                "SELECT v.version, v.ref_id, v.chrom, v.path, v.n_records, v.created_at "
                "FROM meta m JOIN versions v ON v.version = CAST(m.v AS INTEGER) "
                "WHERE m.k = 'current_version'"
            ).fetchone()
        finally:
            con.close()

    def active_ref_info(self):
        """Row of the active reference, or None when no reference is enabled."""
        con = self._connect()
        try:
            return con.execute(
                "SELECT r.id, r.name, r.sha256, r.length, r.path, r.vcf_version, "
                "r.activated_by, r.created_at FROM meta m JOIN reference_sets r "
                "ON r.id = CAST(m.v AS INTEGER) WHERE m.k = 'active_ref_id'"
            ).fetchone()
        finally:
            con.close()

    def all_version_paths(self) -> set:
        con = self._connect()
        try:
            return {r["path"] for r in con.execute("SELECT path FROM versions")}
        finally:
            con.close()

    def all_reference_paths(self) -> set:
        con = self._connect()
        try:
            return {r["path"] for r in con.execute("SELECT path FROM reference_sets")}
        finally:
            con.close()

    def live_instances(self, now: float, dead_after: float) -> set:
        con = self._connect()
        try:
            cutoff = now - dead_after
            return {
                r["id"]
                for r in con.execute("SELECT id FROM instances WHERE last_seen >= ?", (cutoff,))
            }
        finally:
            con.close()

    # ----------------------------------------------------------- state changes

    def activate_reference(self, name: str, sha256: str, length: int, path: str,
                           verified_against_version: int, instance_id: str,
                           expected_version: int, now: float | None = None) -> int:
        """Enable a verified reference iff the VCF version is still as expected.

        Inserts the reference row and flips ``active_ref_id`` in one
        transaction. Raises :class:`StateChanged` when another instance
        published a VCF in the meantime (the caller re-verifies and retries).
        """
        now = time.time() if now is None else now
        with self._txn() as con:
            cur, _ref = self._state(con)
            if cur != expected_version:
                raise StateChanged(cur, _ref)
            cur_row = con.execute(
                "INSERT INTO reference_sets(name, sha256, length, path, vcf_version, "
                "activated_by, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)",
                (name, sha256, length, path, verified_against_version,
                 instance_id, now),
            )
            ref_id = cur_row.lastrowid
            con.execute(
                "UPDATE meta SET v = ? WHERE k = 'active_ref_id'", (str(ref_id),)
            )
            return ref_id

    def publish_vcf(self, base_version: int, base_ref_id: int, chrom: str, path: str,
                    n_records: int, now: float | None = None) -> int:
        """Publish a fully verified VCF as the new current version.

        The whole-file reference check ran against ``(base_version,
        base_ref_id)``; this transaction re-checks that exact combination is
        still current, so a reference swapped in by the other instance can
        never be paired with an unverified VCF.
        """
        now = time.time() if now is None else now
        with self._txn() as con:
            cur, active_ref = self._state(con)
            if cur != base_version or active_ref != base_ref_id:
                raise StateChanged(cur, active_ref)
            new_version = cur + 1
            con.execute(
                "INSERT INTO versions(version, ref_id, chrom, path, n_records, created_at) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (new_version, base_ref_id, chrom, path, n_records, now),
            )
            con.execute(
                "UPDATE meta SET v = ? WHERE k = 'current_version'", (str(new_version),)
            )
            return new_version

    # ------------------------------------------------------------------ refs

    def pin(self, instance_id: str, now: float | None = None):
        """Pin the current (reference, VCF) combination for an in-flight query.

        Returns (version, chrom, path, pin_row_id, ref_id) or None when no VCF
        has been published. The pointer read and pin insert are one
        transaction, so the reaper can never delete the version mid-query.
        """
        now = time.time() if now is None else now
        with self._txn() as con:
            cur, ref = self._state(con)
            if cur == 0:
                return None
            row = con.execute(
                "SELECT chrom, path FROM versions WHERE version = ?", (cur,)
            ).fetchone()
            if row is None:  # pragma: no cover - defensive
                return None
            c = con.execute(
                "INSERT INTO refs(version, instance, created_at) VALUES (?, ?, ?)",
                (cur, instance_id, now),
            )
            return cur, row["chrom"], row["path"], c.lastrowid, ref

    def unpin(self, pin_id: int) -> None:
        with self._txn() as con:
            con.execute("DELETE FROM refs WHERE id = ?", (pin_id,))

    # ------------------------------------------------------------------ heartbeat

    def heartbeat(self, instance_id: str, now: float | None = None) -> None:
        now = time.time() if now is None else now
        with self._txn() as con:
            con.execute(
                "INSERT INTO instances(id, last_seen) VALUES (?, ?) "
                "ON CONFLICT(id) DO UPDATE SET last_seen = excluded.last_seen",
                (instance_id, now),
            )

    def unregister(self, instance_id: str) -> None:
        with self._txn() as con:
            con.execute("DELETE FROM refs WHERE instance = ?", (instance_id,))
            con.execute("DELETE FROM instances WHERE id = ?", (instance_id,))

    # ------------------------------------------------------------------ reaper

    def reap_metadata(self, now: float, dead_after: float):
        """Drop refs of dead instances; delete reclaimable versions and refs.

        A VCF version is reclaimable when it is not the current one and no
        live pin points at it (the reference it was verified against stays in
        the DB as a historical record and does not keep VCF files alive).
        A reference is reclaimable when it is neither active nor referenced
        by a surviving VCF version. Returns
        ``{"versions": [(version, path)], "references": [(id, path)]}``; the
        caller deletes the files after commit, so a crash only leaves
        sweepable orphan directories.
        """
        with self._txn() as con:
            cutoff = now - dead_after
            con.execute(
                "DELETE FROM refs WHERE instance NOT IN "
                "(SELECT id FROM instances WHERE last_seen >= ?)",
                (cutoff,),
            )
            cur, active_ref = self._state(con)
            vrows = con.execute(
                "SELECT version, path FROM versions "
                "WHERE version != ? "
                "AND NOT EXISTS (SELECT 1 FROM refs WHERE refs.version = versions.version)",
                (cur,),
            ).fetchall()
            con.executemany(
                "DELETE FROM versions WHERE version = ?", [(r["version"],) for r in vrows]
            )
            rrows = con.execute(
                "SELECT id, path FROM reference_sets "
                "WHERE id != ? "
                "AND NOT EXISTS (SELECT 1 FROM versions WHERE versions.ref_id = reference_sets.id)",
                (active_ref,),
            ).fetchall()
            con.executemany(
                "DELETE FROM reference_sets WHERE id = ?", [(r["id"],) for r in rrows]
            )
            return {
                "versions": [(r["version"], r["path"]) for r in vrows],
                "references": [(r["id"], r["path"]) for r in rrows],
            }
