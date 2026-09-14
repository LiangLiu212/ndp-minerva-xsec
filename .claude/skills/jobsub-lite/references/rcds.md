# RCDS — shipping code payloads via CVMFS

RCDS (Rapid Code Distribution Service) is how FIFE distributes job payloads.
When you submit with `--tar_file_name dropbox:///abs/path.tar`, jobsub_lite
uploads the tarball to RCDS, which **unpacks it** into one of four CVMFS repos:

    /cvmfs/fifeuser{1,2,3,4}.opensciencegrid.org/sw/<group>/<sha256-hash>/

- The repo number (1–4) is randomly assigned per upload; the hash is printed in
  the jobsub_submit output as `Publishing hash <group>/<hex>`.
- The content is **unpacked** — there is no `.tar` on CVMFS. Workers read files
  (and run binaries) directly out of that directory; CVMFS caches blocks on
  each worker node, so 1000 jobs share one distribution instead of 1000
  transfers.
- Identical tarball content re-uses the same hash dir (content-addressed) —
  re-submitting with an unchanged tarball is cheap.
- On the worker, `$INPUT_TAR_FILE` points into the unpacked dir (dirname = the
  payload root; the skill's worker skeleton resolves this).

## Publish-once / reuse-many (what `jobsub.py publish` automates)

For campaign use you want the CVMFS path *itself*, so many later submissions
can reference the same published payload without re-uploading:

1. Submit **one sentinel job** (`templates/publish_only.sh`) with
   `--tar_file_name dropbox://<tar>` — this triggers the RCDS upload.
2. Parse `Publishing hash <group>/<hex>` from the submit output.
3. When the sentinel drains, `jobsub_fetchlog` it and grep
   `PUBLISH_SENTINEL_CVMFS_DIR=/cvmfs/...` from its stdout (fallback: probe
   `/cvmfs/fifeuserN.../sw/<group>/<hash>` for N=1..4 on the login node).
4. Store `label → cvmfs_dir` in `.jobsub/catalog.json`.

Later submissions pass the dir straight to the worker
(`submit --tar-label L -- -R @TAR_DIR@ ...`) — no `--tar_file_name`, no upload,
instant starts.

## Garbage collection — payloads expire

RCDS prunes unused uploads after roughly **30 days**. Consequences:

- `jobsub.py verify --label L` before reusing an old label: `ok` (<21 d),
  `warn` (21–28 d), `republish` (>28 d or dir already gone).
- `submit --tar-label` refuses when the dir is missing and warns when stale.
- Republishing the identical tarball restores the same hash; the catalog entry
  is refreshed with `publish --overwrite`.
- A previously-published-but-uncataloged tarball can be adopted with
  `jobsub.py adopt --label L --jobid <old jobid>` (parses the hash from the
  recorded submit log).

## Choosing a payload route

| Situation | Route |
|---|---|
| A few config/input files, possibly per-job | `-f /pnfs/...` or `-f file:///abs` (lands in `$CONDOR_DIR_INPUT`) |
| Software install / env reused by many jobs or campaigns | tarball → `publish` → `--tar-label` (RCDS + catalog) |
| One-off submission, no reuse planned | `submit --tar-file X.tar` (dropbox upload tied to that cluster) |
| Very large inputs (events, flux files ≫ GB, per-job unique) | stage to `/pnfs`, `ifdh cp` inside the worker |

Size guidance: RCDS handles multi-GB tarballs, but publication time grows with
size and jobsub may reject ≳10 GB. Trim aggressively (drop `.git`, test
outputs, anything not needed at runtime); the `tarball build` excludes exist
for exactly this. Multi-GB builds: use `--background`.

## Practical notes

- The tarball is unpacked flat: pack with the *tree structure the worker
  expects* (`tar` arcnames become paths under the hash dir). `tarball build
  --include a b c` keeps each entry as a top-level dir inside the payload.
- CVMFS is read-only on workers — anything the job writes goes to
  `$_CONDOR_SCRATCH_DIR` (its CWD), then `ifdh cp` out.
- CVMFS propagation to a login node can lag a few minutes behind the publish;
  `verify` may briefly report `missing` on a *fresh* publish — the workers
  themselves will see it.
- Login-node visibility of `fifeuserN` repos is needed for the hash-probe
  fallback; if none are mounted, `verify` reports `unknown` rather than lying.
