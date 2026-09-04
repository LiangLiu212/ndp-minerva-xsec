"""Small I/O helpers shared across the package: manifests, hashing, run-dir writing."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

try:  # PyYAML is optional: manifests may be written as JSON where it is absent
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


def load_yaml_or_json(path: str | Path):
    path = Path(path)
    text = path.read_text()
    if path.suffix.lower() in (".yaml", ".yml"):
        if yaml is None:
            raise ImportError(f"PyYAML is needed to read {path}; install it or use a .json manifest")
        return yaml.safe_load(text)
    return json.loads(text)


def dump_yaml_or_json(obj, path: str | Path) -> None:
    path = Path(path)
    if path.suffix.lower() in (".yaml", ".yml") and yaml is not None:
        path.write_text(yaml.safe_dump(obj, sort_keys=False))
    else:
        path.write_text(json.dumps(obj, indent=2, default=_json_default))


def _json_default(o):
    import numpy as np
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    if isinstance(o, (set, tuple)):
        return list(o)
    raise TypeError(f"not JSON serialisable: {type(o)}")


def dump_json(obj, path: str | Path) -> None:
    Path(path).write_text(json.dumps(obj, indent=2, default=_json_default))


def sha256_file(path: str | Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def cheap_fingerprint(path: str | Path, head: int = 1 << 20) -> dict:
    """Size + mtime + sha256 of the first and last MiB.

    A full sha256 of a 20 GB AnaTuple costs minutes; this fingerprint costs
    milliseconds and still detects truncation, replacement and most corruption.
    The manifest records which kind was used, so nobody mistakes one for the other.
    """
    path = Path(path)
    st = path.stat()
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read(head))
        if st.st_size > 2 * head:
            fh.seek(st.st_size - head)
            h.update(fh.read(head))
    return {"path": str(path), "size_bytes": st.st_size, "mtime": int(st.st_mtime),
            "fingerprint": "sha256(head1MiB+tail1MiB)", "sha256_partial": h.hexdigest()}


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def git_state(repo: str | Path) -> dict:
    """{'sha': ..., 'dirty': bool} for a git checkout, or {'sha': None} if not a repo."""
    try:
        sha = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True,
                             text=True, timeout=20)
        if sha.returncode != 0:
            return {"sha": None, "dirty": None}
        dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain"], capture_output=True,
                               text=True, timeout=20)
        return {"sha": sha.stdout.strip(), "dirty": bool(dirty.stdout.strip())}
    except Exception:  # pragma: no cover
        return {"sha": None, "dirty": None}


def versions() -> dict:
    import numpy
    out = {"python": sys.version.split()[0], "numpy": numpy.__version__}
    for mod in ("scipy", "uproot", "awkward", "matplotlib", "yaml"):
        try:
            m = __import__(mod)
            out[mod] = getattr(m, "__version__", "present")
        except ImportError:
            out[mod] = None
    return out


def timestamp() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def ensure_dir(p: str | Path) -> Path:
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


def unique_run_dir(root: str | Path, slug: str) -> Path:
    """runs/<YYYY-MM-DD>_<slug>[_N] — never overwrites an existing run."""
    root = ensure_dir(root)
    base = f"{time.strftime('%Y-%m-%d')}_{slug}"
    cand = root / base
    n = 1
    while cand.exists():
        n += 1
        cand = root / f"{base}_{n}"
    cand.mkdir()
    return cand
