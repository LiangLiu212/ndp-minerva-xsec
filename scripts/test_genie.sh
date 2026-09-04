#!/usr/bin/env bash
# Smoke test of the in-repo GENIE: 50 numu CC events on C12 at 3 GeV, then gst conversion.
set -euo pipefail
: "${GENIE:?run via pixi}"
W="$NDP_EXTERNAL/genie-smoke"; mkdir -p "$W"; cd "$W"
XS=()
if [ -n "${GENIEXSECFILE:-}" ] && [ -f "$GENIEXSECFILE" ]; then XS=(--cross-sections "$GENIEXSECFILE"); echo "splines: $GENIEXSECFILE"; else echo "no spline file: gevgen will integrate cross sections on the fly (slow)"; fi
gevgen -n 50 -p 14 -t 1000060120 -e 3.0 --tune G18_02a_00_000 --event-generator-list CC --seed 7 -r 1 -o smoke.ghep.root \
       --message-thresholds Messenger_laconic.xml "${XS[@]}" > gevgen.log 2>&1
gntpc -i smoke.ghep.root -f gst -o smoke.gst.root --tune G18_02a_00_000 > gntpc.log 2>&1
python - <<'PY'
import uproot; t = uproot.open("smoke.gst.root")["gst"]; a = t.arrays(["Ev","El","cc"], library="np")
print(f"gst entries: {t.num_entries}, all CC: {bool(a['cc'].all())}, <El> = {a['El'].mean():.2f} GeV")
PY
echo "GENIE smoke test OK ($W)"
