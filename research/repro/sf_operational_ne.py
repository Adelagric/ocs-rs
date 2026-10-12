"""Time verified dense optima at the exact caps calibrated by the R benchmark.

Usage: python3 research/repro/sf_operational_ne.py wheat mouse
REPRO_OUTPUT_DIR and REPRO_INPUT_DIR select the run's artifact/input directories.
Ne_c=1/sum(c^2), proxy_delta_f_pct=50/Ne_c; neither is realized inbreeding.
"""
import argparse
import csv
import os
from pathlib import Path
import sys
import time

import numpy as np
from dense_support import load_panel, support_first
from operational_data import validate_native


def run(panel):
    outdir = Path(os.environ.get("REPRO_OUTPUT_DIR", "artifacts"))
    source = outdir / ("operational_ne_pig.csv" if panel == "pig" else "operational_ne_realpanels.csv")
    with source.open() as stream:
        caps = [r for r in csv.DictReader(stream) if r["panel"] == panel]
    if len(caps) != 3 or {float(r["target_ne"]) for r in caps} != {25., 50., 100.}:
        raise ValueError(f"missing calibrated operating points for {panel}")
    for row in caps:
        validate_native(row)
    K, bv, male, hashes = load_panel(panel)
    result = []
    for row in caps:
        if any(row[key] != value for key, value in hashes.items()):
            raise ValueError(f"{panel}: exports do not match R benchmark fingerprints")
        ub = float(row["ub"])
        ts = []
        for _ in range(3):
            start = time.perf_counter()
            c, S, iterations = support_first(K, bv, male, ub)
            ts.append(time.perf_counter() - start)
        ne = 1 / np.sum(c ** 2)
        if not np.isfinite(ne) or abs(ne - float(row["target_ne"])) / float(row["target_ne"]) > .02:
            raise ValueError(f"{panel}: dense result did not reach the contribution-Ne target")
        if abs(ne - float(row["achieved_ne"])) > .02 * float(row["target_ne"]):
            raise ValueError(f"{panel}: dense and matrix-free optima disagree on contribution Ne")
        if abs(float(bv @ c) - float(row["gain_ocsrs"])) > 1e-5 * max(1., abs(float(bv @ c))):
            raise ValueError(f"{panel}: dense and matrix-free objective values disagree")
        result.append(dict(panel=panel, target_ne=row["target_ne"], achieved_ne=ne,
                           proxy_delta_f_pct=50 / ne, ub=ub, support=len(S),
                           algorithm_s=float(np.median(ts)), iterations=iterations,
                           status="Solved", verified=True, **hashes))
        print(f"{panel}: target Ne_c={row['target_ne']} actual={ne:.2f} support={len(S)} "
              f"given-G median={np.median(ts):.4f}s verified=True")
    dest = outdir / f"algorithm_operational_{panel}.csv"
    with dest.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panels", nargs="+", choices=("wheat", "mouse", "pig"))
    for panel in parser.parse_args().panels:
        run(panel)
