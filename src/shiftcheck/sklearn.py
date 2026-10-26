"""Optional scikit-learn estimator-style interface for shiftcheck."""
from copy import deepcopy

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.utils.validation import check_is_fitted

from shiftcheck.core.validate import prepare
from shiftcheck.report.report import DiagnosticReport, diagnose as _diagnose


class ShiftCheck(BaseEstimator):
    """Compare new samples against a reference dataset using sklearn conventions.

    This is a fitted monitor, not a transformer or predictive estimator, so it
    is not intended to be placed in a sklearn ``Pipeline``.

    Parameters
    ----------
    diagnose_kwargs : dict, optional
        Keyword arguments forwarded to :func:`shiftcheck.diagnose`.
    """

    def __init__(self, diagnose_kwargs=None):
        self.diagnose_kwargs = diagnose_kwargs

    def fit(self, X, y=None):
        """Store and validate the reference sample; return this monitor."""
        options = {} if self.diagnose_kwargs is None else self.diagnose_kwargs
        prepared = prepare(
            X,
            X,
            feature_names=options.get("feature_names"),
            min_n=1,
            tie_warn=1.0,
            categorical_features=options.get("categorical_features"),
        )
        self.reference_ = deepcopy(X)
        self.n_features_in_ = len(prepared.names) + len(prepared.categorical_names)
        if hasattr(X, "columns") and all(isinstance(name, str) for name in X.columns):
            self.feature_names_in_ = np.asarray(X.columns, dtype=object)
        return self

    def compare(self, X) -> DiagnosticReport:
        """Diagnose distribution shift from the fitted reference to ``X``."""
        check_is_fitted(self, "reference_")
        options = {} if self.diagnose_kwargs is None else self.diagnose_kwargs
        return _diagnose(self.reference_, X, **options)
