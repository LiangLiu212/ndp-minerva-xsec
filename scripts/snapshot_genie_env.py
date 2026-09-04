#!/usr/bin/env python3
"""Write the activated pixi/GENIE environment to a JSON the platform's generator runner loads.

Run inside pixi:  pixi run snapshot-genie-env   (-> external/genie_env.json)
The platform (ndp/theory/generator.py::load_genie_env) passes this dict verbatim as the
environment of gevgen/gntpc, so GENIE never sees whichever Python launched the pipeline.
"""
import json, os, subprocess, sys, time
from pathlib import Path

out = Path(sys.argv[1] if len(sys.argv) > 1 else "external/genie_env.json")
env = {k: v for k, v in os.environ.items() if not k.startswith(("PIXI_", "_", "PS1", "OLDPWD", "PWD", "SHLVL", "TERM"))}
genie = env.get("GENIE", "")
version = subprocess.run(["git", "-C", genie, "describe", "--tags", "--always"], capture_output=True, text=True).stdout.strip() if genie else ""
try:
    banner = subprocess.run(["gevgen", "-h"], env=env, capture_output=True, text=True, timeout=120)
    text = banner.stdout + banner.stderr
    ok = "Syntax" in text and "error while loading shared libraries" not in text
except Exception:
    ok = False
doc = {"created": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "genie_version": version, "genie_dir": genie,
       "root_version": subprocess.run(["root-config", "--version"], capture_output=True, text=True).stdout.strip(),
       "gevgen_responds": ok, "env": env}
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(doc, indent=1, sort_keys=True))
print(f"wrote {out} (GENIE {version}, gevgen responds: {ok})")
