import numpy as np
import pytest
from sklearn.base import clone
from sklearn.exceptions import NotFittedError

from shiftcheck.report.report import DiagnosticReport
from shiftcheck.sklearn import ShiftCheck


def test_shiftcheck_obeys_estimator_parameter_and_fit_conventions():
    options = {
        "n_boot": 20,
        "n_perm": 20,
        "max_n_decomp": 100,
        "n_boot_equivalence": 20,
    }
    estimator = ShiftCheck(diagnose_kwargs=options)
    cloned = clone(estimator)
    assert cloned.get_params() == estimator.get_params()
    assert cloned.set_params(diagnose_kwargs=options) is cloned

    with pytest.raises(NotFittedError):
        estimator.compare(np.zeros((30, 2)))

    rng = np.random.default_rng(5)
    reference = rng.normal(size=(40, 2))
    target = rng.normal(loc=0.5, size=(40, 2))
    assert estimator.fit(reference) is estimator
    assert estimator.n_features_in_ == 2

    report = estimator.compare(target)
    assert isinstance(report, DiagnosticReport)
    assert report.n_source == 40
    assert report.n_target == 40
    assert estimator.reference_ is not reference


def test_shiftcheck_records_dataframe_feature_names():
    pd = pytest.importorskip("pandas")
    data = pd.DataFrame({"left": np.arange(25.0), "right": np.arange(25.0)})
    estimator = ShiftCheck(diagnose_kwargs={"n_boot": 10, "n_perm": 10})
    estimator.fit(data)
    assert estimator.n_features_in_ == 2
    assert estimator.feature_names_in_.tolist() == ["left", "right"]
