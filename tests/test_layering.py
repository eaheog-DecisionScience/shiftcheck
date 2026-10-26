"""The package's own dependency structure is a tested contract (acyclic layers)."""
import ast
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1] / "src" / "shiftcheck"
ALLOWED = {
    "core": set(),
    "uncertainty": {"core"},
    "simulate": {"core"},
    "viz": {"core"},
    "marginal": {"core", "uncertainty"},
    "transform": {"core", "uncertainty"},
    "dependence": {"core", "uncertainty"},
    "decomposition": {"core", "uncertainty"},
    "power": {"core", "uncertainty"},
    "benchmarks": {"core", "uncertainty", "simulate", "marginal", "transform", "dependence",
                   "decomposition", "power", "report"},
    "report": {"core", "uncertainty", "marginal", "transform", "dependence", "decomposition", "power"},
}


def _imports(path):
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("shiftcheck."):
            yield node.module.split(".")[1]
        elif isinstance(node, ast.Import):
            for a in node.names:
                if a.name.startswith("shiftcheck."):
                    yield a.name.split(".")[1]


def test_layers_respect_allowed_dependencies():
    for pkg, allowed in ALLOWED.items():
        for f in (ROOT / pkg).rglob("*.py"):
            bad = set(_imports(f)) - allowed - {pkg}
            assert not bad, f"{f.relative_to(ROOT)} imports forbidden layer(s): {bad}"
