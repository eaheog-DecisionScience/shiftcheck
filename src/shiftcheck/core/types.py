"""Shared result types. Layer 0: no intra-package imports."""
from dataclasses import dataclass, asdict


@dataclass(frozen=True)
class Estimate:
    """A scalar with its uncertainty and provenance."""
    value: float
    lower: float
    upper: float
    method: str
    n_resamples: int

    def to_dict(self):
        return asdict(self)
