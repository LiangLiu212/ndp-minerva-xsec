"""Bridge to the MINERvA exploration repo's certified `benchmark/` engine.

The platform does not re-implement the covariance-aware chi2 or the release ingestion;
it imports them from the audited repo (`tools/xsec_benchmark.md` is the contract). This
module isolates the sys.path handling and caches the ingested release as an .npz.
"""
from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

import numpy as np


def _import_benchmark(minerva_repo: Path):
    repo = str(Path(minerva_repo))
    if repo not in sys.path:
        sys.path.insert(0, repo)
    import benchmark.schema as schema      # noqa
    import benchmark.ingest as ingest      # noqa
    import benchmark.chi2 as chi2          # noqa
    import benchmark.predictions as pred   # noqa
    return schema, ingest, chi2, pred


class PaperRelease:
    """One published result: manifest + ingested arrays in the paper's GlobalID basis."""

    def __init__(self, minerva_repo: str | Path, arxiv: str, cache_dir: str | Path | None = None):
        self.repo = Path(minerva_repo)
        self.arxiv = arxiv
        self.schema, self.ingest_mod, self.chi2, self.pred = _import_benchmark(self.repo)
        self.manifest_path = self.repo / "benchmark" / "papers" / f"{arxiv}.yaml"
        if not self.manifest_path.exists():
            raise FileNotFoundError(self.manifest_path)
        self.manifest = self.schema.load_manifest(self.manifest_path)
        self.schema.validate_manifest(self.manifest)
        self.cache_dir = Path(cache_dir) if cache_dir else self.repo / "bench" / "_ndp_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._npz = self.cache_dir / f"ingest_{arxiv}.npz"
        self.ingest_notes: list[str] = []
        if not self._npz.exists():
            try:
                res = self.ingest_mod.ingest(self.manifest_path)
            except ModuleNotFoundError as e:
                if "ROOT" not in str(e):
                    raise
                # The release's optional .root precision cross-check needs PyROOT. Without it,
                # ingest the same manifest minus that block (the engine supports releases that
                # ship no .root copy) and say so in the provenance.
                import copy, yaml
                m2 = copy.deepcopy(self.manifest)
                m2.get("covariances", {}).pop("root_precision", None)
                alt = self.cache_dir / f"{arxiv}.no_root_precision.yaml"
                alt.write_text(yaml.safe_dump(m2, sort_keys=False))
                res = self.ingest_mod.ingest(alt)
                self.ingest_notes.append("root_precision cross-check skipped: PyROOT not importable in this environment")
            self.ingest_mod.save_npz(res, self._npz)
            (self._npz.with_suffix(".notes.json")).write_text(__import__("json").dumps(self.ingest_notes))
        else:
            notes = self._npz.with_suffix(".notes.json")
            if notes.exists():
                self.ingest_notes = __import__("json").loads(notes.read_text())
        z = np.load(self._npz)
        self.data = z["data"]
        self.mask = z["mask"]
        self.cov_total = z["cov_total"]
        self.cov = {k[4:]: z[k] for k in z.files if k.startswith("cov_")}
        self.pt_edges = z["pt_edges"] if "pt_edges" in z.files else None
        self.pz_edges = z["pz_edges"] if "pz_edges" in z.files else None
        b = self.manifest["basis"]
        self.n_cells, self.n_pt, self.n_pz = b["n_cells"], b["n_pt"], b["n_pz"]
        self.formula = b.get("global_bin_formula", "ipt*n_pz + ipz")
        self.areas = (self.chi2.cell_areas_from_edges(self.pt_edges, self.pz_edges, self.n_pz, self.n_cells,
                                                      n_pt=self.n_pt, formula=self.formula)
                      if self.pt_edges is not None else None)

    def shipped_models(self) -> list[str]:
        return [m["name"] for m in self.manifest.get("models", [])]

    def shipped_curve(self, name: str) -> np.ndarray:
        return self.pred.load_model_txt(self.manifest, name)

    def compare(self, model_vec: np.ndarray, cov: np.ndarray | None = None, mask: np.ndarray | None = None) -> dict:
        cov = self.cov_total if cov is None else cov
        mask = self.mask if mask is None else mask
        return self.chi2.compare(np.asarray(model_vec, float), self.data, cov, mask, areas=self.areas)

    def project_1d(self, vec: np.ndarray, cov: np.ndarray | None, axis: str):
        return self.chi2.project_1d(vec, cov, axis, self.pt_edges, self.pz_edges, formula=self.formula)
