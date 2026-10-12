"""Measured build times and structured contribution-Ne scaling.

Requires both scaling_matrixfree.csv (build timing, independent-locus panels) and
scaling_structured.csv (AlphaSimR solve timing at target contribution Ne_c=50).
Memory curves are theoretical float64 array footprints, not measured peak RAM.
Dashed/open points project an n^2 build trend beyond the benchmark's allocation
policy; a policy omission is not evidence of mathematical infeasibility.
"""
import argparse
import csv
import os
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def read_rows(path):
    if not path.exists():
        raise ValueError(f"{path} not found; run the corresponding reproduction stage")
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{path} contains no measurements")
    return rows


def load_data(directory, reported=False):
    rows = read_rows(directory / "scaling_matrixfree.csv")
    structured = read_rows(directory / "scaling_structured.csv")
    for row in rows + structured:
        for key in ("n", "m", "support"):
            if not np.isfinite(float(row[key])) or float(row[key]) <= 0:
                raise ValueError("scaling coordinates and support must be finite and positive")
    marker_counts = {int(r["m"]) for r in rows + structured}
    if len(marker_counts) != 1 or min(marker_counts) <= 0:
        raise ValueError("build and structured solve curves require the same positive marker count m")
    for row in structured:
        if not np.isfinite(float(row["ne"])) or not np.isfinite(float(row["solve_s"])):
            raise ValueError("structured measurements must be finite")
        target = 50 if reported else float(row["target_ne"])
        if target != 50 or abs(float(row["ne"]) - 50) > .75:
            raise ValueError("structured scaling requires verified target/achieved contribution Ne_c=50")
        if (not reported and row["status"] != "Solved") or float(row["solve_s"]) <= 0:
            raise ValueError("structured scaling contains an unverified or invalid solve")
    for row in rows:
        if row["status"] != "Solved":
            raise ValueError("matrix-free scaling contains an unsuccessful solve")
    built, omitted = [], []
    for row in rows:
        if row["g_build_s"] in ("not_measured", "infeasible"):
            # The older CSV's 'infeasible' meant the same allocation-policy cutoff.
            omitted.append(row)
        elif not np.isfinite(float(row["g_build_s"])):
            raise ValueError("dense build timing must be finite")
        elif float(row["g_build_s"]) > 0:
            built.append(row)
        else:
            raise ValueError("dense build timings must be positive")
    if not built:
        raise ValueError("at least one measured dense build is required for the figure")
    return sorted(built, key=lambda r: int(r["n"])), omitted, structured, marker_counts.pop()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reported", action="store_true", help="plot hash-checked archived measurements")
    args = parser.parse_args()
    directory = Path(os.environ.get("REPRO_OUTPUT_DIR", ROOT / "artifacts"))
    try:
        if args.reported:
            from repro.reported_data import reported_directory
            directory, _ = reported_directory()
        built, omitted, rows, m = load_data(directory, args.reported)
    except (ValueError, KeyError) as error:
        sys.exit(str(error))
    n = np.array([float(r["n"]) for r in rows])
    support = np.array([float(r["support"]) for r in rows])
    solve_s = np.array([float(r["solve_s"]) for r in rows])
    build_n = np.array([float(r["n"]) for r in built])
    build_s = np.array([float(r["g_build_s"]) for r in built])
    omitted_n = sorted(float(r["n"]) for r in omitted if float(r["n"]) > build_n[-1])
    gib = 2. ** 30
    nn = np.logspace(np.log10(min(n.min(), build_n.min()) * .8),
                     np.log10(max(n.max(), build_n.max()) * 1.5), 300)
    plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans",
                        "axes.spines.top": False, "axes.spines.right": False,
                        "savefig.dpi": 200})
    green, red = "#2e7d32", "#c0392b"
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.8, 4.4))
    a.set_title("A   Measured build and structured solve times", loc="left", fontweight="bold")
    a.set_xscale("log")
    a.set_yscale("log")
    a.plot(build_n, build_s, "o-", color=red, label="dense G build: measured (independent loci)")
    if omitted_n:
        projection_n = np.array([build_n[-1]] + omitted_n)
        projection_s = build_s[-1] * (projection_n / build_n[-1]) ** 2
        a.plot(projection_n, projection_s, "--", color=red, label="n² build projection; omitted by policy")
        a.plot(projection_n[1:], projection_s[1:], "o", mfc="white", mec=red)
    a.plot(n, solve_s, "s-", color=green, label=("historical structured solve: Ne_c≈50" if args.reported else "verified structured solve: Ne_c=50"))
    a.set_xlabel("candidates n (log scale)")
    a.set_ylabel("time (s, log scale)")
    a.legend(frameon=False, fontsize=8.3)
    a.text(.03, .97, f"m={m}; support {int(support.min())}–{int(support.max())}\n"
           "ΔF proxy = 1%; realized inbreeding not measured",
           transform=a.transAxes, va="top", fontsize=8.3, color=green)
    b.set_title("B   Theoretical array footprints", loc="left", fontweight="bold")
    b.set_xscale("log")
    b.set_yscale("log")
    b.plot(nn, nn ** 2 * 8 / gib, color=red, label="dense float64 G: 8n² bytes")
    b.plot(nn, nn * m * 8 / gib, color=green, label="dense float64 Z: 8nm bytes")
    b.plot(n, n ** 2 * 8 / gib, "o", color=red, ms=4)
    b.plot(n, n * m * 8 / gib, "o", color=green, ms=4)
    b.axhline(16, ls=":", color="#777", label="16 GiB reference capacity")
    b.set_xlabel("candidates n (log scale)")
    b.set_ylabel("array footprint (GiB, log scale)")
    b.legend(frameon=False, fontsize=8.3)
    b.text(.03, .97, "Array storage only; excludes solver workspace.\n"
           "Packed genotypes use a different footprint.",
           transform=b.transAxes, va="top", fontsize=8.3, bbox=dict(facecolor="white", edgecolor="none", alpha=.85))
    if args.reported:
        fig.suptitle("Reported v0.4.0 summaries; historical solver status not recorded", fontsize=9)
    fig.tight_layout()
    figure_dir = Path(os.environ.get("REPRO_FIGURE_DIR", ROOT / "research"))
    figure_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "png"):
        fig.savefig(figure_dir / f"fig_scaling.{suffix}", bbox_inches="tight")
    print(f"wrote {figure_dir}/fig_scaling.{{pdf,png}}")


if __name__ == "__main__":
    main()
