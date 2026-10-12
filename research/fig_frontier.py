"""Support and matrix-free time along a contribution-Ne proxy frontier.

The proxy 50/Ne_c percent is not realized genomic inbreeding. Fresh data require
all three panels and all three targets. --reported uses hash-checked historical
summaries without claiming a new timing run or certified native solver status.
"""
import argparse
import os
import csv
import pathlib
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--reported", action="store_true", help="plot hash-checked archived measurements")
args = parser.parse_args()
directory = pathlib.Path(os.environ.get("REPRO_OUTPUT_DIR", ROOT / "artifacts"))
if args.reported:
    from repro.reported_data import reported_directory
    directory, _ = reported_directory()
CSV = directory / "table2_operational.csv"
if not CSV.exists():
    sys.exit(f"{CSV} not found")

with CSV.open() as stream:
    rows = list(csv.DictReader(stream))
panels = ["wheat", "mouse", "pig"]
for panel in panels:
    points = [r for r in rows if r["panel"] == panel]
    if len(points) != 3 or {float(r["target_ne"]) for r in points} != {25., 50., 100.}:
        sys.exit(f"complete wheat/mouse/pig frontier required: missing points for {panel}")
    for row in points:
        proxy = row.get("proxy_delta_f_pct", row.get("delta_f_pct"))
        if not args.reported and "proxy_delta_f_pct" not in row:
            sys.exit("fresh frontier requires explicit contribution-proxy fields")
        if any(not np.isfinite(float(value)) or float(value) <= 0
               for value in (proxy, row["matrixfree_s"], row["support"])):
            sys.exit("frontier measurements must be finite and positive")
        if not args.reported and (row.get("native_status") != "Solved" or
                                  row.get("algorithm_status") != "Solved" or
                                  row.get("algorithm_verified") != "True"):
            sys.exit("fresh frontier requires verified native/dense statuses")
        row["proxy_delta_f_pct"] = proxy
labels = {"wheat": "wheat (m=1.3k)", "mouse": "mouse (m=10k)", "pig": "pig (m=53k)"}
colours = {"wheat": "#2e7d32", "mouse": "#1a4e8a", "pig": "#c0392b"}

by = {p: sorted([r for r in rows if r["panel"] == p], key=lambda r: float(r["proxy_delta_f_pct"])) for p in panels}

plt.rcParams.update(
    {
        "font.size": 11,
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "savefig.dpi": 200,
    }
)
fig, (axA, axB) = plt.subplots(1, 2, figsize=(10.8, 4.4))

# ---- Panel A: support |S| vs ΔF ----
axA.set_title("A   Support grows as the cap tightens", loc="left", fontweight="bold", fontsize=11.5)
for p in panels:
    df = [float(r["proxy_delta_f_pct"]) for r in by[p]]
    supp = [float(r["support"]) for r in by[p]]
    axA.plot(df, supp, "o-", color=colours[p], lw=2, ms=6, label=labels[p])
axA.set_xlabel("contribution proxy  50/Ne_c  (%)")
axA.set_ylabel("optimal support  |S|")
axA.invert_xaxis()  # tighter caps (smaller ΔF, larger Nₑ) to the right
axA.legend(loc="upper left", frameon=False, fontsize=9.5)

# ---- Panel B: matrix-free solve TIME vs ΔF (the m-law, cleanly) ----
# Time ∝ iters · n · m, so it orders by marker count (pig on top) and rises as the cap
# tightens — unlike the speed-up vs optiSel, which also carries optiSel's n-dependence.
axB.set_title("B   Matrix-free time rises with markers and tighter caps", loc="left", fontweight="bold", fontsize=11.5)
for p in panels:
    df = [float(r["proxy_delta_f_pct"]) for r in by[p]]
    tm = [float(r["matrixfree_s"]) for r in by[p]]
    axB.plot(df, tm, "s-", color=colours[p], lw=2, ms=6, label=labels[p])
axB.set_yscale("log")
axB.set_xlabel("contribution proxy  50/Ne_c  (%)")
axB.set_ylabel("matrix-free solve time  (s, log)")
axB.invert_xaxis()
axB.legend(loc="upper left", frameon=False, fontsize=9.5)

if args.reported:
    fig.suptitle("Reported v0.4.0 summaries; no new simulation or benchmark", fontsize=9)
fig.tight_layout()
figure_dir = pathlib.Path(os.environ.get("REPRO_FIGURE_DIR", ROOT / "research"))
figure_dir.mkdir(parents=True, exist_ok=True)
fig.savefig(figure_dir / "fig_frontier.pdf", bbox_inches="tight")
fig.savefig(figure_dir / "fig_frontier.png", bbox_inches="tight")
print("wrote research/fig_frontier.pdf and research/fig_frontier.png")
