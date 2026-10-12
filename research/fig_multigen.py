"""Recorded recurrent-selection gain and observed heterozygosity.

OCS targets contribution Ne_c=50; this is not a measured genomic inbreeding rate.
Fresh trajectories require three replicates, both methods and twelve generations.
--reported redraws hash-checked historical summaries without a new simulation.
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
CSV = directory / "multigen_validation.csv"
if not CSV.exists():
    sys.exit(f"{CSV} not found")

with CSV.open() as stream:
    rows = list(csv.DictReader(stream))
if not rows or {r["method"] for r in rows} != {"ocs", "truncation"}:
    sys.exit("both multigeneration methods are required")
for row in rows:
    if not np.isfinite(float(row["gain"])) or not np.isfinite(float(row["het"])):
        sys.exit("multigeneration gain and heterozygosity must be finite")
    if not 0 <= float(row["het"]) <= 1:
        sys.exit("heterozygosity must be between zero and one")
keys = [(r["rep"], r["method"], r["gen"]) for r in rows]
if len(keys) != len(set(keys)):
    sys.exit("duplicate multigeneration observations")
reps = {r["rep"] for r in rows}
if len(reps) != 3 or any({int(r["gen"]) for r in rows if r["rep"] == rep and r["method"] == method} != set(range(13))
                         for rep in reps for method in ("ocs", "truncation")):
    sys.exit("complete three-replicate, twelve-generation trajectories required")
gens = sorted({int(r["gen"]) for r in rows})
methods = ["ocs", "truncation"]
style = {"ocs": ("#1a4e8a", "-", "support-first OCS (Ne_c=50 proxy target)"),
         "truncation": ("#c0392b", "--", "truncation")}


def mean_by_gen(method, field):
    return [np.mean([float(r[field]) for r in rows if r["method"] == method and int(r["gen"]) == g]) for g in gens]


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

axA.set_title("A   Genetic gain", loc="left", fontweight="bold", fontsize=11.5)
for m in methods:
    col, ls, lab = style[m]
    axA.plot(gens, mean_by_gen(m, "gain"), ls, color=col, lw=2.2, label=lab)
axA.set_xlabel("generation")
axA.set_ylabel("mean genetic value")
axA.legend(loc="upper left", frameon=False, fontsize=9.5)

axB.set_title("B   Diversity retained (heterozygosity)", loc="left", fontweight="bold", fontsize=11.5)
for m in methods:
    col, ls, lab = style[m]
    axB.plot(gens, mean_by_gen(m, "het"), ls, color=col, lw=2.2, label=lab)
axB.set_xlabel("generation")
axB.set_ylabel("observed heterozygosity")
axB.legend(loc="lower left", frameon=False, fontsize=9.5)

if args.reported:
    fig.suptitle("Reported v0.4.0 summaries; no new simulation or benchmark", fontsize=9)
fig.tight_layout()
figure_dir = pathlib.Path(os.environ.get("REPRO_FIGURE_DIR", ROOT / "research"))
figure_dir.mkdir(parents=True, exist_ok=True)
fig.savefig(figure_dir / "fig_multigen.pdf", bbox_inches="tight")
fig.savefig(figure_dir / "fig_multigen.png", bbox_inches="tight")
print("wrote research/fig_multigen.pdf and research/fig_multigen.png")
