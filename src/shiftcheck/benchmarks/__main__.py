import argparse
from shiftcheck.benchmarks.runner import run

ap = argparse.ArgumentParser(description="shiftcheck benchmark suite")
ap.add_argument("--level", default="quick", choices=["smoke", "quick", "full"])
ap.add_argument("--jobs", type=int, default=1)
ap.add_argument("--out", default="benchmark_results")
ap.add_argument("--only", nargs="*", default=None, help="experiment names to run")
a = ap.parse_args()
run(a.level, a.jobs, a.out, only=a.only)
