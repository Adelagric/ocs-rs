"""Figure 2: the operating-point frontier — support size and speed-up along the rate
of inbreeding ΔF, on the three real panels. Makes §3.6 (the most panel-specific result)
visual: the support grows as the cap tightens, and the matrix-free speed-up follows the
marker count m (winning on the low-m panels, dropping below 1× on the 52k-SNP pig),
while the given-G algorithm stays ahead throughout.

Reads artifacts/table2_operational.csv. Outputs research/fig_frontier.{pdf,png}.
"""
import csv
import pathlib
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = pathlib.Path(__file__).resolve().parents[1] / "artifacts" / "table2_operational.csv"
if not CSV.exists():
    sys.exit(f"{CSV} not found")

rows = list(csv.DictReader(CSV.open()))
panels = ["wheat", "mouse", "pig"]
labels = {"wheat": "wheat (m=1.3k)", "mouse": "mouse (m=10k)", "pig": "pig (m=53k)"}
colours = {"wheat": "#2e7d32", "mouse": "#1a4e8a", "pig": "#c0392b"}

by = {p: sorted([r for r in rows if r["panel"] == p], key=lambda r: float(r["delta_f_pct"])) for p in panels}

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
    df = [float(r["delta_f_pct"]) for r in by[p]]
    supp = [float(r["support"]) for r in by[p]]
    axA.plot(df, supp, "o-", color=colours[p], lw=2, ms=6, label=labels[p])
axA.set_xlabel("rate of inbreeding  ΔF  (%)")
axA.set_ylabel("optimal support  |S|")
axA.invert_xaxis()  # tighter caps (smaller ΔF, larger Nₑ) to the right
axA.legend(loc="upper left", frameon=False, fontsize=9.5)

# ---- Panel B: matrix-free solve TIME vs ΔF (the m-law, cleanly) ----
# Time ∝ iters · n · m, so it orders by marker count (pig on top) and rises as the cap
# tightens — unlike the speed-up vs optiSel, which also carries optiSel's n-dependence.
axB.set_title("B   Matrix-free time rises with markers and tighter caps", loc="left", fontweight="bold", fontsize=11.5)
for p in panels:
    df = [float(r["delta_f_pct"]) for r in by[p]]
    tm = [float(r["matrixfree_s"]) for r in by[p]]
    axB.plot(df, tm, "s-", color=colours[p], lw=2, ms=6, label=labels[p])
axB.set_yscale("log")
axB.set_xlabel("rate of inbreeding  ΔF  (%)")
axB.set_ylabel("matrix-free solve time  (s, log)")
axB.invert_xaxis()
axB.legend(loc="upper left", frameon=False, fontsize=9.5)

fig.tight_layout()
fig.savefig("research/fig_frontier.pdf", bbox_inches="tight")
fig.savefig("research/fig_frontier.png", bbox_inches="tight")
print("wrote research/fig_frontier.pdf and research/fig_frontier.png")
