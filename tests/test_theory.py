import numpy as np
from _helpers import toy_truth, site, skip
from ndp.theory import flux
from ndp.theory.models import ModelSpec, _weights
from ndp.channels import load_channel
from ndp.events import INT_CODE


def test_flux_table_parses_and_integrates():
    p = site().repo_root / "resources/flux/arXiv2110.13372_supplemental.txt"
    t = flux.load_flux_table(p)
    d = flux.per_cm2_per_pot_per_gev(t["values"], t["units"])
    phi = flux.integrated_flux(t["edges"], d, 0, 100)
    assert t["n_rows"] == 128 and t["edges"][0] == 0.0 and t["edges"][-1] == 100.0
    # sanity band vs the published integrated flux 6.32e-8 (2106.16210); table is 2-decimal truncated
    assert 0.95 < phi / 6.32e-8 < 1.02
    # flux_average of a constant is the constant
    assert abs(flux.flux_average(t["edges"], d, lambda e: np.full_like(e, 3.0)) - 3.0) < 1e-9


def test_reweight_expressions():
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    t = toy_truth(2000)
    w = _weights({"weight_expr": "1.0"}, t, ch)
    assert w.shape == (2000,) and np.all(w == 1.0)
    w = _weights({"weight_expr": "where(int_type == MEC, 1.5, 1.0)"}, t, ch)
    assert np.all(w[t["int_type"] == INT_CODE["MEC"]] == 1.5) and np.all(w[t["int_type"] != INT_CODE["MEC"]] == 1.0)
    w = _weights({"weight_expr": "1 + 0.1*lep_pT"}, t, ch)
    assert np.all(w >= 1.0)


def test_model_spec_validation():
    assert ModelSpec.from_dict({"name": "a", "kind": "shipped_curve", "curve": "x"}).validate() == []
    assert ModelSpec.from_dict({"name": "b", "kind": "reweight", "base": "reference_mc"}).validate() != []
    try:
        ModelSpec.from_dict({"name": "c", "kind": "nonsense"})
        assert False, "unknown kind accepted"
    except ValueError:
        pass
