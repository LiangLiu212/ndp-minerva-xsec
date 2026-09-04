"""Loader for the analysis branch catalog (config/branches.json).

The catalog is the single source of truth for which AnaTuple branches the
analysis uses and why (see its _meta.description). Code never lists branch
names directly — it asks the catalog by role, so documentation and code
cannot drift.
"""
import json
from functools import lru_cache
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CATALOG_PATH = REPO_ROOT / "config" / "branches.json"


@lru_cache(maxsize=1)
def load(path=None):
    return json.loads(Path(path or CATALOG_PATH).read_text())


def branches_for_role(role, path=None):
    """Names of all catalogued branches carrying `role`, in catalog order."""
    cat = load(path)
    known = set(cat["_meta"]["roles"])
    if role not in known:
        raise KeyError(f"unknown role {role!r}; catalog declares {sorted(known)}")
    return tuple(b["name"] for b in cat["branches"] if role in b["roles"])


def conventions(path=None):
    """The _meta.conventions block (e.g. mc_current CC/NC encoding)."""
    return load(path)["_meta"]["conventions"]
