#!/bin/bash
# stage_payload.sh — assemble the RCDS payload directory for the grid campaign:
#   env/            the slim pixi "grid" environment (pixi install -e grid), copied (relocatable: $ORIGIN rpaths)
#   ndp/ channels/ measurements/ resources/ grid/   the platform code + manifests
#   ndp_version.json  git state of the checkout the payload was built from
# Usage: grid/stage_payload.sh [staging_dir]   (default: /exp/dune/data/users/$USER/ndp-grid-payload)
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
STAGE=${1:-/exp/dune/data/users/${USER}/ndp-grid-payload}
ENV_DIR="${HERE}/.pixi/envs/grid"
if [ ! -x "${ENV_DIR}/bin/python3" ]; then
  echo "grid environment missing: run 'pixi install -e grid' in ${HERE} first" >&2; exit 1
fi
mkdir -p "${STAGE}"
echo "staging payload in ${STAGE}"
rsync -a --delete --exclude '__pycache__' --exclude '*.pyc' "${ENV_DIR}/" "${STAGE}/env/"
for d in ndp channels measurements resources grid; do
  rsync -a --delete --exclude '__pycache__' --exclude '*.pyc' --exclude 'campaigns' "${HERE}/${d}/" "${STAGE}/${d}/"
done
mkdir -p "${STAGE}/worklists"
python3 - "${HERE}" "${STAGE}" <<'EOF'
import json, subprocess, sys, time
here, stage = sys.argv[1], sys.argv[2]
def git(*a):
    r = subprocess.run(["git", "-C", here, *a], capture_output=True, text=True); return r.stdout.strip() if r.returncode == 0 else None
json.dump({"sha": git("rev-parse", "HEAD"), "dirty": bool(git("status", "--porcelain")), "staged": time.strftime("%Y-%m-%dT%H:%M:%S%z")},
          open(f"{stage}/ndp_version.json", "w"), indent=1)
EOF
du -sh "${STAGE}" "${STAGE}/env"
echo "sanity check with a scrubbed environment:"
env -i PATH="${STAGE}/env/bin:/usr/bin:/bin" LD_LIBRARY_PATH="${STAGE}/env/lib" PYTHONPATH="${STAGE}" PYTHONNOUSERSITE=1 NDP_CONFIG="${STAGE}/grid/ndp_grid.yaml" \
  "${STAGE}/env/bin/python3" -c "import numpy, uproot, awkward, XRootD, fsspec_xrootd, yaml, ndp; from ndp.channels import load_channel; load_channel('minerva_me_ccqelike_1mu1p'); print('payload imports ok')"
echo "next: python3 <skill>/scripts/jobsub.py tarball build --build-dir ${STAGE} --include env ndp channels measurements resources grid worklists ndp_version.json --exclude-component __pycache__ --name-prefix ndp-stream --background"
