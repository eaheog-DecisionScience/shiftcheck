"""shiftcheck: characterize distribution shift between a source and a target dataset.

Dependency layers (enforced by tests/test_layering.py):
  L0 core         -> (nothing)
  L1 uncertainty, simulate, viz -> core
  L2 marginal, transform, dependence, decomposition, power -> core, uncertainty
  L3 report       -> core, uncertainty, all L2 engines
"""
from shiftcheck.report.report import diagnose, DiagnosticReport
from shiftcheck.simulate.generator import make_shift_pair

__all__ = ["diagnose", "DiagnosticReport", "make_shift_pair"]
__version__ = "0.4.1"
