"""Shared fixtures for the platform tests (pytest-free so the fallback runner can use them)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ndp.events import TruthTable  # noqa: E402
from ndp.config import load_site_config  # noqa: E402


class Skip(Exception):
    """Raised by a test that cannot run in this environment (mirrors pytest.skip)."""


def skip(msg):
    try:
        import pytest
        pytest.skip(msg)
    except ImportError:
        raise Skip(msg)


def site():
    return load_site_config()


def have_minerva_repo():
    cfg = site()
    return cfg.minerva_repo is not None and (cfg.minerva_repo / "benchmark" / "papers" / "2106.16210.yaml").exists()


def have_mc_cache():
    cfg = site()
    return cfg.data_dir is not None and (cfg.data_dir / "cache" / "truth_mc110040.npz").exists() \
        and (cfg.data_dir / "cache" / "reco_mc110040.npz").exists()


def have_data_file():
    cfg = site()
    return cfg.data_dir is not None and (cfg.data_dir / "MasterAnaDev_data_AnaTuple_run00010066_Playlist.root").exists()


def toy_truth(n=20000, seed=0, with_fs=False) -> TruthTable:
    """A synthetic numu CC sample with a crude but sensible kinematic shape (GeV)."""
    rng = np.random.default_rng(seed)
    E = rng.gamma(4.0, 1.6, n)                       # peaks ~5 GeV
    y = rng.beta(1.2, 2.5, n)
    El = np.maximum(E * (1 - y), 0.12)
    th = np.abs(rng.normal(0, 0.12, n))
    p = np.sqrt(np.maximum(El ** 2 - 0.1057 ** 2, 1e-6))
    phi = rng.uniform(0, 2 * np.pi, n)
    it = rng.choice([1, 2, 3, 5], size=n, p=[0.2, 0.3, 0.45, 0.05])
    cols = {
        "nu_pdg": np.full(n, 14), "E_nu": E, "lep_pdg": np.full(n, 13),
        "lep_px": p * np.sin(th) * np.cos(phi), "lep_py": p * np.sin(th) * np.sin(phi), "lep_pz": p * np.cos(th), "lep_E": El,
        "current": np.ones(n, int), "int_type": it, "target_Z": np.full(n, 6), "target_A": np.full(n, 12),
        "Q2": 2 * E * El * (1 - np.cos(th)), "W": rng.uniform(0.94, 3.0, n), "weight": np.ones(n),
        "vtx_x": rng.uniform(-800, 800, n), "vtx_y": rng.uniform(-800, 800, n), "vtx_z": rng.uniform(6000, 8400, n),
    }
    if with_fs:
        counts = rng.integers(1, 5, n)
        tot = int(counts.sum())
        pdg = rng.choice([2212, 2112, 211, -211, 111], size=tot, p=[0.4, 0.3, 0.1, 0.1, 0.1])
        KE = rng.exponential(0.3, tot)
        mass = np.where(np.abs(pdg) == 211, 0.13957, np.where(pdg == 111, 0.13498, 0.9383))
        cols.update({"fs_offsets": np.concatenate([[0], np.cumsum(counts)]), "fs_pdg": pdg, "fs_E": KE + mass,
                     "fs_px": np.zeros(tot), "fs_py": np.zeros(tot), "fs_pz": np.sqrt((KE + mass) ** 2 - mass ** 2)})
    return TruthTable(cols, {"source": "toy", "norm": {"kind": "xsec_per_nucleon", "xsec_per_unit_weight": 4.0e-38 / n}})
