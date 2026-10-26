import json
from shiftcheck import diagnose, make_shift_pair

FAST = dict(n_boot=60, n_perm=60, max_n_decomp=300)


def run(kinds, seed=1, n=1500, diag=None, **kw):
    S, T, truth = make_shift_pair(n_source=n, n_target=n, d=3, kinds=kinds, seed=seed, **kw)
    return diagnose(S, T, seed=seed, **{**FAST, **(diag or {})}), truth


def test_no_shift_not_flagged_and_certified_equivalent():
    r, _ = run([], n=4000)
    p = r.profile()
    assert p["marginal_difference"] == "equivalent_within_tolerance"
    assert p["dependence_difference"] == "equivalent_within_tolerance"
    assert p["shift_attribution"] == "not_detected"


def test_small_sample_no_shift_is_inconclusive_not_equivalent():
    r, _ = run([], n=80)
    p = r.profile()
    assert p["marginal_difference"] != "detected"
    assert p["marginal_difference"] == "inconclusive"


def test_affine_shift_marginal_only():
    r, _ = run(["affine"], strength=0.4)
    p = r.profile()
    assert p["marginal_difference"] == "detected"
    v = [f["verdict"] for f in r.transform_compatibility["features"]]
    assert "nonlinear_monotone" not in v          # never mislabeled nonlinear
    assert v.count("affine_compatible") >= 2      # individual features may be 'inconclusive' at n=1500
    assert p["dependence_difference"] != "detected"
    assert p["shift_attribution"] == "marginal_dominated"


def test_monotone_shift_flagged_nonlinear():
    r, _ = run(["monotone"], strength=1.0)
    p = r.profile()
    assert p["marginal_difference"] == "detected"
    assert p["marginal_form"] == "nonlinear_monotone_detected"
    assert p["dependence_difference"] != "detected"


def test_copula_shift_without_marginal_shift():
    r, _ = run(["copula"], diag={"max_n_decomp": 800})
    p = r.profile()
    assert p["marginal_difference"] != "detected"
    assert p["dependence_difference"] == "detected"
    assert p["shift_attribution"] == "dependence_dominated"


def test_combined_shifts_reported_independently():
    r, _ = run(["monotone", "copula"], strength=1.0)
    p = r.profile()
    assert p["marginal_form"] == "nonlinear_monotone_detected"
    assert p["dependence_difference"] == "detected"
    d = r.decomposition
    assert d["marginal"]["value"] > 0 and d["dependence"]["value"] > 0


def test_decomposition_contributions_track_truth():
    r_m, _ = run(["monotone"], strength=0.6)
    r_d, _ = run(["copula"])
    assert r_m.decomposition["marginal_share"]["value"] > 0.9
    assert r_d.decomposition["marginal_share"]["value"] < 0.2


def test_json_roundtrip():
    r, _ = run(["affine"])
    assert json.loads(r.to_json())["profile"]["outcome_mechanism"] == "not_assessable_from_features"


def test_power_module_mde_shrinks_with_sample_size():
    S_small, T_small, _ = make_shift_pair(n_source=300, n_target=300, d=3, kinds=[], seed=5)
    S_big, T_big, _ = make_shift_pair(n_source=3000, n_target=3000, d=3, kinds=[], seed=5)
    r_s = diagnose(S_small, T_small, seed=5, power=True, power_reps=10, **FAST)
    r_b = diagnose(S_big, T_big, seed=5, power=True, power_reps=10, **FAST)
    mde = lambda r: r.power["marginal"]["features"][0]["mde_location_shift_sd"]
    assert mde(r_s) is not None and mde(r_b) is not None and mde(r_b) < mde(r_s)
    # theory sanity: KS at alpha/3 with n=3000/sample resolves a few hundredths of an SD up to ~0.1 SD
    assert mde(r_b) < 0.12
