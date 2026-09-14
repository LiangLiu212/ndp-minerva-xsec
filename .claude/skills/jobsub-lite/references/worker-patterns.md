# Worker script patterns

Start every worker from `templates/worker_skeleton.sh`. Its skeleton encodes
conventions the monitoring machinery relies on — keep them.

## Anatomy (and why each part exists)

1. `set -e` + `date; hostname; uname -r` — fail fast, and stamp enough context
   to debug "works on one node, dies on another" (CPU/kernel differences).
2. **getopts, not positional args** — jobsub appends worker args after the
   `file://` executable; flags keep them order-independent. Reserve `-R`
   (payload dir) and `-O` (PNFS output base); add project flags to the string.
3. **Seed = `CLUSTER * 100000 + PROCESS`** — unique across processes *and*
   across submissions, reproducible, no coordination needed. Plain
   `CLUSTER + PROCESS` is unique only within one cluster: cluster ids are
   schedd-global sequential ints, so two submissions overlap whenever
   Δcluster < N — combining their outputs then double-counts RNG streams.
   Caveat: cluster ids run ~1e7–1e8, so these seeds reach ~1e12–1e13, past
   int32. Bash arithmetic is 64-bit, but if your generator truncates seeds
   to 32 bits, reduce deliberately (and record how) instead of overflowing
   silently. Print the seed (`echo SEED=...`) so a fetched log fully
   determines the job.
4. **ifdh setup** — `source /cvmfs/larsoft.opensciencegrid.org/setup-env.sh &&
   spack load --first ifdhc` is the standard FIFE pattern and safe on DUNE/
   uboone workers. Load it *before* the payload env so payload libs can't
   shadow ifdh's Python.
5. **Payload resolution** — `-R <cvmfs dir>` override first (published
   payload), else derive from `$INPUT_TAR_FILE` (dropbox payload). Never
   hard-code a fifeuser path in the worker; it changes per publish.
6. **Run in `$PWD`** — that's `$_CONDOR_SCRATCH_DIR`, private and wiped.
   Programs that scatter files into CWD are fine here; never write to the
   payload (CVMFS is read-only) or to `/pnfs` directly.
7. **Outputs** — per-process subdir `<pnfs-out>/<%04d PROCESS>/` via
   `ifdh mkdir_p` + `ifdh cp -D`. The %04d prefix is what `pull` and the
   PNFS-side completion count walk. Copy outputs *before* any optional/fragile
   trailing steps.
8. **`DONE` sentinel** — a bare `DONE` line as the last statement. Completion
   counting prefers grepping fetched logs for it; only jobs that got all the
   way through — including copying ≥1 output file (zero glob matches exits
   4) — print it (`set -e` guarantees that).

## Log hygiene

Grid logs are fetched per process — a chatty program at `-N 100` is gigabytes
of text. Silence verbose generators/frameworks (rate limits, laconic modes,
`>/dev/null` for known-noisy steps) but keep the skeleton's own `echo` markers:
they are the grep targets for debugging and completion.

## Relocatable payloads (conda/pixi/anything with RPATHs)

Installs built under an absolute prefix carry that prefix in RPATHs and
scripts. On a worker the prefix doesn't exist, which is *survivable*:

- The dynamic loader tries RPATH first, fails (path absent), then falls back
  to `LD_LIBRARY_PATH` — so export `LD_LIBRARY_PATH` with the payload's lib
  dirs and compiled binaries run fine.
- Re-create the install's activation by hand in the env section, with the
  payload dir substituted for the original root (`PATH`, `LD_LIBRARY_PATH`,
  and the package-specific vars the software needs).
- Python entry-point *scripts* with baked shebangs (`#!/orig/prefix/bin/python`)
  do NOT survive — invoke them as `"$PAYLOAD/env/bin/python3" script` instead.
- Sanity-check on the login node: `ldd` the main binaries against a scrubbed
  env (`env -i LD_LIBRARY_PATH=... ldd ...`) before burning grid time.

## CPU-architecture portability

Binaries compiled with `-march=native` (or on a newer microarch than the
worker) die with SIGILL. Two fixes, use either or both:

- Rebuild with a floor: `-march=x86-64-v2` (or `-v3`) instead of `native`.
- Constrain workers to match the build host:
  `--append_condor_requirements '(TARGET.Microarch>="x86_64-v3")'`
  (set it once in `.jobsub/config.json` as `append_condor_requirements`).

conda-forge/system packages are generic; the risk is exclusively locally
compiled code.

## Resource requests

- `--memory`: measure locally (`/usr/bin/time -v`, MaxRSS) + ~30% headroom.
  Underasking = held jobs; gross overasking = fewer slots.
- `--expected-lifetime`: per-*process* walltime, not the campaign. Shorter
  classes schedule faster. Prefer more processes × shorter jobs, but keep each
  process ≳15–30 min so container startup and payload caching amortize.
- `--disk`: outputs + scratch of one process.

## Debug loop

1. `jobsub.py status <jobid>` — held? check `jobsub_q --long` `HoldReason`.
2. `jobsub.py fetchlog <jobid>` — read the `.out`/`.err` of the failing
   process; the skeleton's echo markers bracket each phase.
3. Rerun a single process's command **locally** with the same seed/flags in a
   scratch dir (payload from CVMFS if published) — most failures reproduce.
4. Fix, `-N 1` smoke, then resubmit the campaign.
