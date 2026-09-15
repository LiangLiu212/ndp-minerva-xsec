# grid/ — streamed processing of the MINERvA Open Data AnaTuples on the FNAL grid

One job process streams a few AnaTuples over XRootD (`root://fndcadoor.fnal.gov:1095//pnfs/...`),
runs `ndp.grid.process_file` per file (= `ndp data cache` + the channel's truth signal cutflow +
reco selection cutflow + derived-column skims + a `manifest_<tag>.json` sidecar) in the payload's
slim pixi environment, and `ifdh cp`s the outputs to PNFS. The login node plans, tracks, harvests
and merges (`python -m ndp grid ...`, `python -m ndp data merge ...`); consumers then read playlist
products through `ndp.products` (channel `data.playlists`).

Measured on the EAF node (2026-09-14, run 110040 / 10066 streamed with the grid env): MC file
70 s wall, 1.26 GB peak RSS, 332 MB requested of 21.5 GB; data file 4 s, 3.7 MB. Outputs per MC
file: truth 197 MB + reco 40 MB + reco-truth 62 MB (archive) + skims 14 + 10 MB + sidecar; every
array byte-identical to the local `ndp data cache` build.

## Files
- `worker.sh` — the HTCondor worker (jobsub-lite skeleton conventions; flags `-R -O -W -K -C -n -t`).
- `ndp_grid.yaml` — site config for workers (`NDP_CONFIG`).
- `stage_payload.sh` — assembles the RCDS payload (`.pixi/envs/grid` + ndp/ channels/ measurements/
  resources/ grid/ + `ndp_version.json`) in `/exp/dune/data/users/$USER/ndp-grid-payload` and checks it
  imports with a scrubbed environment.
- `campaigns/<name>/` — `campaign.json`, `worklists/*.txt`, `resubmit_*.txt`, `sidecars/` (gitignored).

## Campaign recipe (all commands from the platform root; the skill CLI with system python3)
```bash
kinit                                             # kerberos; the bearer token is fetched by jobsub_lite at submit
pixi install -e grid && grid/stage_payload.sh     # slim env + code -> payload dir
SKILL=.claude/skills/jobsub-lite/scripts/jobsub.py
python3 $SKILL init --group dune                  # once: .jobsub/config.json
python3 $SKILL tarball build --build-dir /exp/dune/data/users/$USER/ndp-grid-payload \
    --include env ndp channels measurements resources grid worklists ndp_version.json --exclude-component __pycache__ --name-prefix ndp-stream
python3 $SKILL publish --tarball .jobsub/tarballs/ndp-stream_<hash>.tar --label ndp-stream-v1   # one sentinel job, ~5-20 min
python3 $SKILL verify --label ndp-stream-v1

python -m ndp grid harvest-pot                    # per-file POT table (resources/minerva/opendata_pot_per_file.tsv)
python -m ndp grid plan fhc_2026-09 --channel minerva_me_ccqelike_1mu1p --beams FHC   # worklists per (beam, kind, playlist)
python -m ndp grid stage-worklists fhc_2026-09          # upload worklists to PNFS scratch (-f /pnfs/... job inputs; -f file:// fails off-site)
python -m ndp grid submit-cmd fhc_2026-09 FHC_mc_1A --max-processes 20     # prints the jobsub-lite submit line; run it
python3 $SKILL status <jobid>; python3 $SKILL fetchlog <jobid>
python -m ndp grid status fhc_2026-09             # sidecars on PNFS -> per-file done/failed/incomplete (POT + entry checks)
python -m ndp grid resubmit fhc_2026-09           # resubmit_N_<wl>.txt of what is not done; stage-worklists --files <it>, then submit-cmd --file <it>
python -m ndp grid harvest fhc_2026-09 --playlists 1A [--no-archive]   # PNFS -> <data_dir>/products/FHC/1A/files/<tag>/ (no-archive: skims+reco only, ~65 MB per MC file)
python -m ndp data merge --beam FHC --playlist 1A --kind data,mc   # playlist products + pot_1A_mc.json / pot_1A_data.json
```
Then point the channel at the products (`data: {beam: FHC, playlists: {mc: [1A], data: [1A]}}`) and run
`ndp signal` / `ndp selection` / `ndp surrogate build` / `ndp run` as usual.

Concurrency: stay at or below 20 processes for the first playlist, then at most 50, never more than
100 simultaneous streams from the OpenData dCache pools (their capacity is unknown; read the
per-job MB/s from the sidecars' `timings_s` / `bytes_requested` before scaling up).

Harvest and status need PNFS access from the login node: `ifdh` (CVMFS larsoft spack) or
`xrdfs`/`xrdcp` against `root://fndca1.fnal.gov:1094` with a DUNE bearer token
(`htgettoken -a htvaultprod.fnal.gov -i dune`).

## GiBUU generation campaigns (`ndp gibuu ...`, worker `gibuu_worker.sh`)

A generator sample for a model spec of `kind: gibuu` (e.g. `models/gibuu_2025_me_fhc_c12.yaml`) is
N independent GiBUU jobs, one seed each, merged by `ndp.theory.gibuu.merge_jobs` into
`runs/_generator_cache/gibuu_<fingerprint>/truth.npz` (weights prescaled by 1/N so they still sum to
the mean flux-averaged cross section per nucleon). The payload is the in-repo `GiBUU.x`, the four
shared libraries it needs from the pixi env, the `buuinput2025` tables (65 MB) and the card template
whose only run-time placeholders are the input path, the flux file and the seed.

```bash
python -m ndp gibuu smoke models/gibuu_2025_me_fhc_c12.yaml --channel minerva_me_ccqelike_1mu1p --ensembles 100   # local job, ~30 s
python -m ndp gibuu plan gibuu_me_c12_2026-09 models/gibuu_2025_me_fhc_c12.yaml --channel minerva_me_ccqelike_1mu1p [--n-jobs N --ensembles E]
grid/stage_gibuu_payload.sh runs/_generator_cache/gibuu_<fingerprint>          # -> /exp/dune/data/users/$USER/ndp-gibuu-payload (+ scrubbed ldd check)
python3 $SKILL tarball build --build-dir /exp/dune/data/users/$USER/ndp-gibuu-payload --include GiBUU.x lib buuinput cards gibuu_payload.json --name-prefix ndp-gibuu
python3 $SKILL publish --tarball .jobsub/tarballs/ndp-gibuu_<hash>.tar --label ndp-gibuu-v1 && python3 $SKILL verify --label ndp-gibuu-v1
python -m ndp gibuu submit-cmd gibuu_me_c12_2026-09 --tar-label ndp-gibuu-v1 [--n 1]   # prints the submit line; run it (-N 1 first)
python -m ndp gibuu record gibuu_me_c12_2026-09 --jobid <jobid> --cluster <cluster> --tar-label ndp-gibuu-v1
python -m ndp gibuu status gibuu_me_c12_2026-09                                  # sidecars on PNFS -> done / failed per process
python -m ndp gibuu submit-cmd gibuu_me_c12_2026-09 --tar-label ndp-gibuu-v1 --processes 0007,0042   # rerun failures with fresh seeds
python -m ndp gibuu harvest gibuu_me_c12_2026-09                                 # PNFS -> runs/_generator_cache/gibuu_<fp>/jobs/<%04d>/
python -m ndp gibuu merge gibuu_me_c12_2026-09 --channel minerva_me_ccqelike_1mu1p   # -> truth.npz (per-job sigma check vs the absorption file)
python -m ndp run models/gibuu_2025_me_fhc_c12.yaml --channel minerva_me_ccqelike_1mu1p --measurement all --modes folded --efficiency-run runs/<eff run>
```
Job outputs per process: `FinalEvents.dat.gz`, `neutrino_absorption_cross_section_ALL.dat` (the
merge refuses a job whose summed weights disagree with it), `neutrino_initialized_energyFlux.dat`,
`job.card`, `gibuu.log.gz`, `manifest_<%04d>.json`. Sizing and the weighted-vs-equal-weights choice:
`.claude/skills/gibuu/references/running.md` (measured on the EAF node, 2026-09-15).
