# Reference-Gated Versioned VCF Variant Store

FastAPI + pysam + SQLite service for uploading single-chromosome VCF files and
querying them by 1-based closed chromosome intervals, with atomic versioned
replacement, reference-sequence verification, and safe concurrent access from
two app instances sharing one data volume.

## Run

```bash
docker compose up --build
```

Two instances come up sharing one named data volume:

- instance 1: http://localhost:8001/docs
- instance 2: http://localhost:8002/docs

Both serve the same published (reference, VCF) combination; uploads and
reference swaps may go to either instance. The interactive Swagger UI at
`/docs` on both ports can change the reference, upload a VCF and run queries.

## Reference gate (supply-time interception)

A reference FASTA must be enabled before any VCF can be published, and every
new input is verified in full against the input it will be paired with:

- **First reference** (`POST /reference`, store empty): format validation
  only. A VCF cannot be published until this succeeds — the first VCF upload is
  rejected with `409` while no reference is active.
- **Replacement reference**: verified against the *current* VCF — chromosome
  name and every `REF` fragment. Any chromosome mismatch, out-of-range span or
  base mismatch rejects the upload (`422`) with the offending record positions
  and reasons; the active reference and every VCF version stay untouched.
- **Replacement VCF**: verified in full against the *active* reference before
  publishing. A mismatch/out-of-range result rejects it (`422`); the current
  version is unchanged.

Reference FASTA acceptance rules:

- at most 10 MiB,
- exactly one FASTA record (single chromosome),
- sequence characters from `A/C/G/T` only (upper case; no `N`, gaps,
  whitespace or blank lines),
- non-empty sequence with a named header.

Verification problems are reported as a list (capped at 100, with a
`truncated` flag) carrying the 1-based VCF record index, chromosome, `POS`, a
machine-readable `reason` (`chromosome_mismatch`, `position_out_of_range`,
`ref_mismatch`), a human `detail`, and — for base mismatches — the expected
reference fragment.

## API

### `POST /reference` (multipart field `file`)

Validate and enable a reference FASTA. Returns `201` with
`{"ref_id", "name", "length", "sha256", "verified_against_version"}`.
`400` for malformed FASTA, `413` over the size limit, `422` when verification
against the current VCF fails (active reference left unchanged).

### `GET /reference`

The active reference (`ref_id`, name, length, sha256, the VCF version it was
verified against, who activated it), or `{"active": false, "ref_id": 0}`.

### `POST /vcf?base_version=<n>` (multipart field `file`)

Accepts a VCF that is:

- at most 10 MiB,
- single-chromosome (all records share one CHROM),
- exactly 8 columns (`#CHROM POS ID REF ALT QUAL FILTER INFO`, no FORMAT/sample
  columns),
- biallelic with `REF`/`ALT` containing only `A/C/G/T`.

Records are sorted by position, BGZF-compressed and Tabix-indexed, then the
whole file is verified against the active reference. The upload carries the
version the client last saw (`base_version`, `0` when the store is empty).
Returns `201` with `{"version", "ref_id", "chrom", "n_records"}`. `409` when
the reference gate is not satisfied or `base_version` is stale, `422` when the
VCF fails reference verification, `400` for malformed VCF, `413` over the limit.

### `GET /variants?chrom=<c>&start=<s>&end=<e>`

1-based closed interval `[start, end]`. Returns every raw VCF record whose
reference span `[POS, POS+len(REF)-1]` intersects the interval — including
deletions that start before the interval. Fetching goes through the Tabix
index (`TabixFile.fetch`), never a full scan. The response carries the fixed
pair it actually read:

```json
{"version": 3, "ref_id": 6, "chrom": "chr1", "start": 1, "end": 1000,
 "count": 11, "records": ["chr1\t1000\t..."]}
```

The `(ref_id, version)` combination is pinned for the whole query, so a
reference swap or VCF replacement that lands during the query can never mix a
new reference into the older response, and vice versa.

### `GET /version`

Current published VCF version, its `ref_id`, chromosome and record count.

## Consistency design

- **Atomic publish of the pair**: the current state is the pair
  `(active_ref_id, current_version)`. Files are fully built/verified in a
  per-instance staging directory, copied into `versions/` or `refs/`, and only
  then is the relevant pointer flipped inside a single `BEGIN IMMEDIATE`
  SQLite transaction. Readers never see a half-built or unverified file set.
- **Compare-and-swap, re-verify instead of blind reject**: the expensive
  verification runs without holding any lock. The final flip re-checks the
  exact `(reference, VCF)` combination the check was based on. If the other
  instance changed one input while a slow check was running, the latest
  combination is re-verified and retried: a still-compatible combination is
  published, an incompatible one is refused, and an unmatched pair can never
  go live. A VCF whose optimistic `base_version` moved still gets `409`.
- **No global serialization / non-blocking checks**: verification runs in
  worker threads (`asyncio.to_thread`) and uses the standard-library gzip
  reader rather than htslib for the sequential scan, so a slow check on one
  instance never blocks the other instance's corrections, reference swaps or
  queries. SQLite's single-writer lock is held only for the short pointer
  flip; WAL keeps readers and the writer independent.
- **Old queries finish on their own pair**: a query pins its
  `(version, ref_id)` before opening files. A background reaper deletes a
  superseded VCF's data/index (and an unreferenced reference) only when no pin
  references it. Instances heartbeat; pins of crashed instances expire.
- **Crash / write-failure recovery**: pointer flips are the last step, so a
  crash or a failed disk write leaves either the complete old combination or a
  fully verified new one. Orphaned `versions/`/`refs/` directories and
  per-instance staging leftovers are swept (with a grace period); the
  previously published combination keeps serving across restarts.

## Data directory layout

```
<data_dir>/meta.sqlite3              metadata, pointers, pins (WAL)
<data_dir>/versions/<uuid>/          published VCF sets (variants.vcf.gz + .tbi)
<data_dir>/refs/<uuid>/              activated references (reference.fasta)
<data_dir>/staging/<instance>/...    in-progress builds
```

## Local development

```bash
pip install -r requirements.txt
python -m pytest tests/          # 73 tests
uvicorn --factory app.main:create_app   # dev instance (uses /data or $VCF_DATA_DIR)
```

Configuration via environment: `VCF_DATA_DIR`, `VCF_INSTANCE_ID`,
`VCF_MAX_UPLOAD_BYTES`, `VCF_REAPER_INTERVAL`, `VCF_HEARTBEAT_INTERVAL`,
`VCF_DEAD_AFTER`, `VCF_ORPHAN_GRACE`.
