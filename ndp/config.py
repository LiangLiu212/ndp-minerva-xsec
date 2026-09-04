"""Site configuration: where the platform finds external resources.

Nothing personal is hard-coded in the package. `ndp.yaml` at the repository root (or
the file named by `$NDP_CONFIG`) declares the site paths; environment variables
override individual entries:

    NDP_MINERVA_REPO   the ndp-minerva-data-release-exploration checkout
                       (paper releases, benchmark engine, selection tool, run artifacts)
    NDP_DATA_DIR       directory holding the MINERvA open-data AnaTuples
    NDP_GENIE_ENV      a genie-agent environment snapshot JSON (config/env/<install>.json)
    NDP_GENIE_SPLINES  a GENIE cross-section spline XML

Every path is optional: features that need a missing resource report it and skip,
they never guess a path.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field, asdict
from pathlib import Path

from .io import load_yaml_or_json

REPO_ROOT = Path(__file__).resolve().parents[1]

_ENV_KEYS = {
    "minerva_repo": "NDP_MINERVA_REPO",
    "data_dir": "NDP_DATA_DIR",
    "genie_env_json": "NDP_GENIE_ENV",
    "genie_splines": "NDP_GENIE_SPLINES",
}


@dataclass
class SiteConfig:
    minerva_repo: Path | None = None
    data_dir: Path | None = None
    genie_env_json: Path | None = None
    genie_splines: Path | None = None
    repo_root: Path = REPO_ROOT
    extra: dict = field(default_factory=dict)

    @property
    def resources(self) -> Path:
        return self.repo_root / "resources"

    @property
    def surrogates(self) -> Path:
        return self.repo_root / "surrogates"

    @property
    def runs(self) -> Path:
        return self.repo_root / "runs"

    @property
    def channels_dir(self) -> Path:
        return self.repo_root / "channels"

    def as_dict(self) -> dict:
        d = asdict(self)
        return {k: (str(v) if isinstance(v, Path) else v) for k, v in d.items()}

    def require(self, key: str) -> Path:
        """Return a configured path or raise a clear error naming the config key."""
        p = getattr(self, key)
        if p is None:
            raise FileNotFoundError(
                f"site config has no '{key}' (set it in ndp.yaml or ${_ENV_KEYS.get(key, key.upper())})")
        p = Path(p)
        if not p.exists():
            raise FileNotFoundError(f"site config '{key}' points to a missing path: {p}")
        return p


def load_site_config(path: str | Path | None = None) -> SiteConfig:
    cfg_path = Path(path) if path else Path(os.environ.get("NDP_CONFIG", REPO_ROOT / "ndp.yaml"))
    raw: dict = {}
    if cfg_path.exists():
        raw = load_yaml_or_json(cfg_path) or {}
    sc = SiteConfig()
    for key, env in _ENV_KEYS.items():
        val = os.environ.get(env) or raw.get(key)
        if val:
            p = Path(os.path.expandvars(os.path.expanduser(str(val))))
            if not p.is_absolute():          # relative site paths are relative to the repository
                p = REPO_ROOT / p
            setattr(sc, key, p)
    sc.extra = {k: v for k, v in raw.items() if k not in _ENV_KEYS}
    return sc
