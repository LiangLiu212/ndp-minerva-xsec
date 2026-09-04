#!/usr/bin/env python3
"""Fallback test runner for environments without pytest: discovers test_* functions."""
import importlib.util
import sys
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
from _helpers import Skip  # noqa: E402

passed = failed = skipped = 0
for f in sorted(HERE.glob("test_*.py")):
    spec = importlib.util.spec_from_file_location(f.stem, f)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception:
        print(f"IMPORT FAIL {f.name}"); traceback.print_exc(); failed += 1; continue
    for name in sorted(n for n in dir(mod) if n.startswith("test_")):
        fn = getattr(mod, name)
        try:
            fn(); passed += 1; print(f"PASS  {f.stem}::{name}")
        except Skip as e:
            skipped += 1; print(f"SKIP  {f.stem}::{name} ({e})")
        except Exception:
            failed += 1; print(f"FAIL  {f.stem}::{name}"); traceback.print_exc()
print(f"\n{passed} passed, {failed} failed, {skipped} skipped")
sys.exit(1 if failed else 0)
