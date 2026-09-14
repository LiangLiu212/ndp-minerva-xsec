"""Build a clean environment for invoking jobsub_lite / ifdh.

jobsub_lite is its own python venv (invoked by absolute path), but pixi/conda
`PYTHONHOME` / `PYTHONPATH` / `PIXI_*` / `CONDA_*` leak into the child and break
it. We **scrub** those poisoning vars from a copy of `os.environ` while passing
*everything else* through, so kerberos/token auth (`KRB5CCNAME`,
`BEARER_TOKEN_FILE`, `XDG_RUNTIME_DIR`, `X509_USER_PROXY`, …) survives
untouched. Auth is live and per-session, so we scrub at call time. Every
`subprocess.run([cfg["jobsub_*"]…])` passes `env=build_submit_env()`.
"""
from __future__ import annotations

import os

# Drop any var whose name starts with one of these (pixi/conda machinery that
# poisons jobsub_lite's own python).
_DROP_PREFIXES = ("PIXI_", "CONDA_", "MAMBA_", "_CE_", "BASH_FUNC_")

# Drop these exact names (interpreter overrides + the conda venv marker).
_DROP_EXACT = frozenset({"PYTHONHOME", "PYTHONPATH", "VIRTUAL_ENV"})


def _is_poison(key: str) -> bool:
    return key in _DROP_EXACT or key.startswith(_DROP_PREFIXES)


def build_submit_env(base: dict[str, str] | None = None) -> dict[str, str]:
    """Return a copy of the environment with pixi/conda/python-override vars
    removed. Everything else (notably auth/runtime vars) passes through."""
    src = os.environ if base is None else base
    return {k: v for k, v in src.items() if not _is_poison(k)}


def leaked_vars(env: dict[str, str]) -> list[str]:
    """Names in `env` that would poison jobsub_lite — empty after scrubbing."""
    return sorted(k for k in env if _is_poison(k))
