import numpy as np

from shiftcheck import diagnose


def test_notebook_html_report_is_readable_and_escapes_feature_names():
    rng = np.random.default_rng(19)
    source = rng.normal(size=(35, 2))
    target = rng.normal(loc=0.3, size=(35, 2))
    report = diagnose(
        source,
        target,
        feature_names=["<img src=x>", "ordinary"],
        n_boot=8,
        n_perm=8,
        max_n_decomp=80,
        n_boot_equivalence=8,
        seed=19,
    )

    rendered = report._repr_html_()

    assert "shiftcheck diagnostic report" in rendered
    assert "Result profile" in rendered
    assert "Continuous marginal shifts" in rendered
    assert "Shift attribution" in rendered
    assert "Full diagnostic data (JSON)" in rendered
    assert "&lt;img src=x&gt;" in rendered
    assert "<img src=x>" not in rendered
