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
