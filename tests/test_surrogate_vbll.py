"""The vbll_event surrogate: smearing machinery, the fold against a binned response built from the model's own
pairs, every grid evaluating on synthetic reco columns, and persistence. Needs torch + vbll:
    pixi install -e ml && pixi run -e ml python -m pytest tests/test_surrogate_vbll.py -q
Skipped in the torch-free default environment."""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

from _helpers import skip

try:
    import torch  # noqa: F401
    import vbll   # noqa: F401
    HAVE_TORCH = True
except ImportError:
    HAVE_TORCH = False

from ndp.channels import load_channel, load_measurement, list_measurements
from ndp.events import TruthTable
from ndp.surrogate.binned import BinnedResponse
from ndp.surrogate.base import load_surrogate

CHANNEL = "minerva_me_ccqelike_1mu1p"
# realistic MeV normaliser statistics (mean, std) per particle and component, of the order of the x60 port's
_STATS = {"muon": ((5300.0, 2700.0), (30.0, 530.0), (-300.0, 560.0), (5250.0, 2670.0)),
          "proton": ((1250.0, 130.0), (-30.0, 400.0), (-50.0, 400.0), (500.0, 250.0))}
NORM = {p: {c: [m, s] for c, (m, s) in zip(("E", "px", "py", "pz"), st)} for p, st in _STATS.items()}


def _need():
    if not HAVE_TORCH:
        skip("torch/vbll not installed (pixi install -e ml)")


def _toy_model(directory: Path, head: str = "standard"):
    """A randomly initialised (untrained) model of the ported architecture, saved and loaded back through the
    platform's own persistence — a deterministic, if physically meaningless, p(reco | truth)."""
    from ndp.surrogate.vbll_model import VBLLSpec, VBLLSurrogateModel, build_module, load_vbll_model
    spec = VBLLSpec(head_type=head, d_embed=4, hidden=16, n_layers=2, n_train_per_particle=1000, frame="detector",
                    normaliser={"input": NORM, "output": NORM}, meta={"source": "toy (random init)"})
    torch.manual_seed(1)
    mod = build_module(spec)
    VBLLSurrogateModel.save_dir(directory, spec, mod.state_dict())
    return load_vbll_model(directory)


def _toy_truth(n: int = 4000, seed: int = 0) -> TruthTable:
    """Generator-like (beam-frame, no geometry) 1mu1p signal events: one muon and one proton inside the windows."""
    rng = np.random.default_rng(seed)
    p_mu = rng.uniform(2.5, 12.0, n); th = rng.uniform(0.0, np.deg2rad(15.0), n); ph = rng.uniform(0, 2 * np.pi, n)
    p_p = rng.uniform(0.55, 1.05, n); thp = rng.uniform(np.deg2rad(10.0), np.deg2rad(60.0), n); php = rng.uniform(0, 2 * np.pi, n)
    E_mu = np.sqrt(p_mu ** 2 + 0.1057 ** 2); E_p = np.sqrt(p_p ** 2 + 0.93827 ** 2)
    cols = {"nu_pdg": np.full(n, 14), "E_nu": E_mu + E_p, "lep_pdg": np.full(n, 13),
            "lep_px": p_mu * np.sin(th) * np.cos(ph), "lep_py": p_mu * np.sin(th) * np.sin(ph), "lep_pz": p_mu * np.cos(th), "lep_E": E_mu,
            "current": np.ones(n, int), "int_type": np.ones(n, int), "target_Z": np.full(n, 6), "target_A": np.full(n, 12),
            "Q2": np.full(n, 0.3), "W": np.full(n, 0.94), "weight": np.ones(n),
            "fs_offsets": np.arange(n + 1), "fs_pdg": np.full(n, 2212), "fs_E": E_p,
            "fs_px": p_p * np.sin(thp) * np.cos(php), "fs_py": p_p * np.sin(thp) * np.sin(php), "fs_pz": p_p * np.cos(thp)}
    return TruthTable(cols, {"source": "toy", "generator": "toy", "has_geometry": False, "frame": "beam",
                             "norm": {"kind": "xsec_per_nucleon", "xsec_per_unit_weight": 1e-38}})


def test_smear_reproduces_the_models_own_predictive():
    _need()
    with tempfile.TemporaryDirectory() as d:
        model = _toy_model(Path(d))
    rng = np.random.default_rng(2)
    truth = np.stack([rng.uniform(3000, 9000, 500), rng.normal(0, 300, 500), rng.normal(-300, 300, 500), rng.uniform(3000, 9000, 500)], 1)
    mean, sig = model.predict("muon", truth.astype(np.float32))
    assert mean.shape == (500, 4) and np.all(sig > 0)
    s = model.smear("muon", truth, n_samples=400, seed=3)
    z = (s - mean[None]) / sig[None]                       # standard head: sigma is deterministic, so z ~ N(0, 1)
    assert abs(z.mean()) < 0.03 and abs(z.std() - 1.0) < 0.03
    assert np.all(np.abs(s.mean(0) - mean) < 5 * sig / np.sqrt(400))


def _binned_from_pairs(measurement, channel, model, t: TruthTable, seed: int = 11) -> BinnedResponse:
    """A binned response learned from one smeared copy of `t` (eff = 1): the empirical response of the model."""
    from ndp.surrogate.vbll import truth_4vectors, synthetic_reco
    mask = channel.in_phase_space(t) & channel.is_signal(t)
    assert mask.all(), "every toy event must be signal in phase space"
    mu, pr = truth_4vectors(channel, t, mask, model.spec.frame)
    r = synthetic_reco(model.smear("muon", mu, 1, seed)[0], model.smear("proton", pr, 1, seed + 1)[0], model.spec.frame)
    xt, yt = measurement.truth_observables(channel, t)
    xr, yr = measurement.reco_observables(r, params=channel.observable_params)
    return BinnedResponse.fit(measurement.binning, x_true_den=xt, y_true_den=yt, x_true_num=xt, y_true_num=yt,
                              x_reco_num=xr, y_reco_num=yr, pot_mc=1.0)


def test_fold_agrees_with_a_binned_response_built_from_the_models_own_pairs():
    _need()
    from ndp.surrogate.vbll import VBLLEventSurrogate
    ch = load_channel(CHANNEL); m = load_measurement(ch, "muon_p"); t = _toy_truth(6000)
    with tempfile.TemporaryDirectory() as d:
        model = _toy_model(Path(d))
    b = _binned_from_pairs(m, ch, model, t)
    w = VBLLEventSurrogate(b, model, m, ch, n_samples=60, seed=5, truncate_to_reco_windows=False)
    pred = w.fold_table(ch, t)
    x, y = m.truth_observables(ch, t)
    true_cells, _, _ = m.binning.histogram(x, y)
    ref = b.fold(true_cells)                                  # the reco histogram of the pairs (eff = 1)
    var = np.maximum(ref, 1.0) * (1.0 + 1.0 / 60)             # Poisson on the single realisation + the fold's own K-sample noise
    chi2 = float(np.sum((pred - ref) ** 2 / var)); ndf = int((ref > 0).sum())
    assert chi2 / ndf < 3.0, (chi2, ndf, pred.round(1), ref)
    assert w.last_fold_info["n_events"] == t.n and w.last_fold_info["sample_pass_fraction"] == 1.0


def test_every_grid_evaluates_on_synthetic_reco_columns():
    _need()
    from ndp.surrogate.vbll import VBLLEventSurrogate, truth_4vectors, synthetic_reco, reco_window_pass, RECO_COLUMNS
    ch = load_channel(CHANNEL); t = _toy_truth(1500, seed=4)
    with tempfile.TemporaryDirectory() as d:
        model = _toy_model(Path(d))
    mu, pr = truth_4vectors(ch, t, np.ones(t.n, bool), model.spec.frame)
    r = synthetic_reco(model.smear("muon", mu, 1, 0)[0], model.smear("proton", pr, 1, 1)[0], model.spec.frame)
    assert set(RECO_COLUMNS) <= set(r) and all(np.isfinite(r[k]).all() for k in RECO_COLUMNS)
    ok = reco_window_pass(ch, r)
    assert ok.dtype == bool and ok.shape == (t.n,)
    names = [n for n in list_measurements(ch) if n != "published"]
    assert len(names) >= 18
    for n in names:
        m = load_measurement(ch, n)
        x, y = m.reco_observables(r, params=ch.observable_params)
        assert np.isfinite(x).all() and np.isfinite(y).all(), n
        b = BinnedResponse(m.binning, eff=np.ones(m.binning.n_cells), P=np.eye(m.binning.n_cells))
        w = VBLLEventSurrogate(b, model, m, ch, n_samples=2, seed=1)
        cells = w.fold_table(ch, t)
        assert cells.shape == (m.binning.n_cells,) and np.all(cells >= 0) and cells.sum() <= t.n + 1e-9, n


def test_wrapper_persists_and_reloads_through_load_surrogate():
    _need()
    from ndp.surrogate.vbll import VBLLEventSurrogate
    from ndp.config import load_site_config
    ch = load_channel(CHANNEL); m = load_measurement(ch, "proton_p"); t = _toy_truth(800, seed=7)
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        model = _toy_model(d / "model")   # noqa: F841  (from_parts reloads it from the directory)
        b = BinnedResponse(m.binning, eff=np.full(m.binning.n_cells, 0.5), P=np.eye(m.binning.n_cells))
        b.save(d / "binned")
        w = VBLLEventSurrogate.from_parts(load_site_config(), ch, m, b, d / "binned", model_dir=d / "model", n_samples=3, seed=9)
        c1 = w.fold_table(ch, t)
        w.save(d / "wrapper")
        w2 = load_surrogate(d / "wrapper")
        assert w2.kind == "vbll_event" and w2.n_samples == 3 and w2.truncate and w2.meta["measurement"] == "proton_p"
        assert np.allclose(w2.fold_table(ch, t), c1)
        assert np.allclose(w2.background(1.0), 0.0)
        assert np.allclose(w2.fold_eff_only(np.ones(m.binning.n_cells)), 0.5)
        assert abs(c1.sum() - 0.5 * t.n * w.last_fold_info["sample_pass_fraction"]) <= 0.5 * t.n * 0.5   # eff 0.5, windows
