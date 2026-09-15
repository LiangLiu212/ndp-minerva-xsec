#!/bin/bash
# stage_gibuu_payload.sh — assemble the RCDS payload for a GiBUU generation campaign:
#   GiBUU.x            the in-repo build (external/gibuu/release2025/objects/GiBUU.x)
#   lib/               the shared libraries it needs from the pixi env (per ldd: libgfortran, libquadmath, libgcc_s, libbz2)
#   buuinput/          external/gibuu/buuinput2025 (65 MB of tables)
#   cards/             <model>.job.tmpl with __PATH_TO_INPUT__/__FLUX_FILE__/__SEED__ left for the worker, + flux_gibuu.dat
#   gibuu_payload.json provenance (binary sha256, model fingerprint, git state)
# Usage: grid/stage_gibuu_payload.sh <gibuu cache dir (runs/_generator_cache/gibuu_<fp>)> [staging_dir]
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
CACHE=${1:?"usage: stage_gibuu_payload.sh <runs/_generator_cache/gibuu_<fingerprint>> [staging_dir]"}
STAGE=${2:-/exp/dune/data/users/${USER}/ndp-gibuu-payload}
GX="${HERE}/external/gibuu/release2025/objects/GiBUU.x"
INP="${HERE}/external/gibuu/buuinput2025"
[ -x "${GX}" ] || { echo "GiBUU.x missing: pixi run build-gibuu" >&2; exit 1; }
[ -f "${CACHE}/card.job.tmpl" ] || { echo "no card template in ${CACHE} (run: ndp gibuu plan <name> <model.yaml>)" >&2; exit 1; }
mkdir -p "${STAGE}/lib" "${STAGE}/cards"
cp "${GX}" "${STAGE}/GiBUU.x"
for lib in $(ldd "${GX}" | awk '$3 ~ /\.pixi\/envs/ {print $3}'); do cp -L "${lib}" "${STAGE}/lib/"; done
rsync -a --delete "${INP}/" "${STAGE}/buuinput/"
cp "${CACHE}/card.job.tmpl" "${STAGE}/cards/$(basename "$(dirname "${CACHE}/card.job.tmpl")")_card.job.tmpl" 2>/dev/null || true
cp "${CACHE}/card.job.tmpl" "${STAGE}/cards/card.job.tmpl"
cp "${CACHE}/flux_gibuu.dat" "${STAGE}/cards/flux_gibuu.dat"
[ -f "${CACHE}/energies.txt" ] && cp "${CACHE}/energies.txt" "${STAGE}/cards/energies.txt"
python3 - "${HERE}" "${STAGE}" "${CACHE}" <<'EOF'
import hashlib, json, subprocess, sys, time
here, stage, cache = sys.argv[1], sys.argv[2], sys.argv[3]
def git(*a):
    r = subprocess.run(["git", "-C", here, *a], capture_output=True, text=True); return r.stdout.strip() if r.returncode == 0 else None
prep = json.load(open(f"{cache}/gibuu_prepare.json"))
json.dump({"gibuu_x_sha256": hashlib.sha256(open(f"{stage}/GiBUU.x", "rb").read()).hexdigest(), "fingerprint": prep["fingerprint"],
           "gibuu_version": prep["gibuu_version"], "flux_file_sha256": prep["flux_file_sha256"], "template_sha256": prep["template_sha256"],
           "git": {"sha": git("rev-parse", "HEAD"), "dirty": bool(git("status", "--porcelain"))}, "staged": time.strftime("%Y-%m-%dT%H:%M:%S%z")},
          open(f"{stage}/gibuu_payload.json", "w"), indent=1)
EOF
du -sh "${STAGE}"
echo "sanity check with a scrubbed environment (the binary must resolve its libraries from lib/ only):"
env -i PATH=/usr/bin:/bin LD_LIBRARY_PATH="${STAGE}/lib" ldd "${STAGE}/GiBUU.x" | grep -v "^\s*linux-vdso\|ld-linux" | awk '{print "  " $1, "->", $3}'
env -i PATH=/usr/bin:/bin LD_LIBRARY_PATH="${STAGE}/lib" ldd "${STAGE}/GiBUU.x" | grep -qi "not found" && { echo "unresolved libraries" >&2; exit 1; } || echo "  all resolved"
echo "next: python3 .claude/skills/jobsub-lite/scripts/jobsub.py tarball build --build-dir ${STAGE} --include GiBUU.x lib buuinput cards gibuu_payload.json --name-prefix ndp-gibuu --background"
