"""Figure 3: multi-generation validation. Recurrent selection over 12 generations,
support-first OCS at ΔF=1% (Nₑ=50) against truncation, on AlphaSimR populations
(3 replicates). OCS preserves heterozygosity and sustains genetic gain — controlling
the rate of inbreeding keeps the variance that long-term response feeds on.

Reads artifacts/multigen_validation.csv. Outputs research/fig_multigen.{pdf,png}.
"""
import csv
import pathlib
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = pathlib.Path(__file__).resolve().parents[1] / "artifacts" / "multigen_validation.csv"
if not CSV.exists():
    sys.exit(f"{CSV} not found")

rows = list(csv.DictReader(CSV.open()))
gens = sorted({int(r["gen"]) for r in rows})
methods = ["ocs", "truncation"]
style = {"ocs": ("#1a4e8a", "-", "support-first OCS (ΔF=1%)"),
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

fig.tight_layout()
fig.savefig("research/fig_multigen.pdf", bbox_inches="tight")
fig.savefig("research/fig_multigen.png", bbox_inches="tight")
print("wrote research/fig_multigen.pdf and research/fig_multigen.png")
