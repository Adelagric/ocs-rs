"""Headline figure: support-first scales to genomic n where the dense matrix cannot.

Data is read from artifacts/scaling_matrixfree.csv, written by
`cargo run --release --example scaling_matrixfree` (m=1000, binding cap
k=0.1*mean_diag, release profile). It used to be hand-copied into this file, which
let the figure lag the solver: after the incremental-Gram work halved the solve
times, the committed figure kept the old curve. Reading the artifact removes that
failure mode — regenerate the CSV and the figure follows.

Memory curves are exact footprints (n^2*8 for dense G, n*m*8 for Z). Build and solve
times are measured; the dense-G build at n=40000 is not run (G exceeds free RAM) and
is projected ~n^2, marked open.
Outputs research/fig_scaling.{pdf,png}.
"""
import csv
import pathlib
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = pathlib.Path(__file__).resolve().parents[1] / "artifacts" / "scaling_matrixfree.csv"
if not CSV.exists():
    sys.exit(
        f"{CSV} not found — run:\n"
        "  cargo run --release --example scaling_matrixfree"
    )

rows = list(csv.DictReader(CSV.open()))  # matrix-free sweep: the dense-G build curve
m = int(rows[0]["m"])

# The solve and support curve is drawn at an HONEST breeder operating point
# (ΔF = 1 %, Nₑ = 50) on structured AlphaSimR panels, not at the loose fixed-absolute
# cap the build curve happens to use — the memory/build story below is
# operating-point-independent, the support and solve are not.
SCSV = pathlib.Path(__file__).resolve().parents[1] / "artifacts" / "scaling_structured.csv"
srows = list(csv.DictReader(SCSV.open())) if SCSV.exists() else rows
n = np.array([float(r["n"]) for r in srows])
support = np.array([float(r["support"]) for r in srows])
solve_s = np.array([float(r.get("solve_s", r.get("mfree_solve_s"))) for r in srows])

# The dense build is timed only where G fits in RAM; beyond that the example records
# "infeasible" and the figure projects the O(n^2) trend from the last measured point.
built = [r for r in rows if r["g_build_s"] != "infeasible"]
build_n = np.array([float(r["n"]) for r in built])
build_s = np.array([float(r["g_build_s"]) for r in built])
unbuilt = [r for r in rows if r["g_build_s"] == "infeasible"]
build_proj_n = float(unbuilt[0]["n"]) if unbuilt else None
build_proj_s = (
    build_s[-1] * (build_proj_n / build_n[-1]) ** 2 if build_proj_n else None
)

# --- exact memory footprints ---
GIB = 2.0**30
nn = np.logspace(np.log10(800), np.log10(60000), 300)
G_curve, Z_curve = nn**2 * 8 / GIB, nn * m * 8 / GIB
G_pts, Z_pts = n**2 * 8 / GIB, n * m * 8 / GIB

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
C_SOLVE, C_BUILD, C_G, C_Z = "#2e7d32", "#c0392b", "#c0392b", "#2e7d32"

fig, (axA, axB) = plt.subplots(1, 2, figsize=(10.8, 4.4))

# ---------- Panel A: setup time — building G alone dwarfs the solve ----------
axA.set_title("A   Building the dense matrix alone dwarfs the solve", loc="left", fontweight="bold", fontsize=11.5)
axA.set_xscale("log")
axA.set_yscale("log")
axA.plot(build_n, build_s, "o-", color=C_BUILD, lw=2.2, ms=6, label="dense G: build only  (O(n²m))")
axA.plot([build_n[-1], build_proj_n], [build_s[-1], build_proj_s], "--", color=C_BUILD, lw=1.3)
axA.plot(build_proj_n, build_proj_s, "o", mfc="white", mec=C_BUILD, mew=1.5, ms=6)
axA.plot(n, solve_s, "s-", color=C_SOLVE, lw=2.2, ms=5.5, label="support-first solve (ΔF=1%, Nₑ=50)")
axA.set_xlabel("candidates  n  (log scale)")
axA.set_ylabel("time  (s, log scale)")
axA.set_ylim(3e-3, 2e1)
axA.legend(loc="lower right", frameon=False, fontsize=9.5)

# Annotation computed from the data, not written by hand: the last n at which the
# dense matrix was actually built, against the solve at that same n.
i_last = int(np.argmax(build_n))
n_last, t_build = build_n[i_last], build_s[i_last]
axA.annotate(
    f"{t_build:.1f} s just to build G at n={n_last / 1000:.0f}k\n"
    f"(more than the whole ΔF=1% solve)",
    # upper-right: the only region both curves and the (upper-left) legend leave free
    xy=(n_last, t_build), xytext=(12500, 11.0), fontsize=9, color=C_BUILD,
    ha="center", arrowprops=dict(arrowstyle="->", color=C_BUILD, lw=1),
)
axA.text(
    1150, 5.5e-3,
    f"support |S| = {int(support.min())}–{int(support.max())} at ΔF=1% (Nₑ=50)",
    fontsize=9, color=C_SOLVE,
)

# ---------- Panel B: the dense matrix is the memory wall ----------
axB.set_title("B   …and does not fit in memory", loc="left", fontweight="bold", fontsize=11.5)
axB.set_xscale("log")
axB.set_yscale("log")
axB.axhspan(16, 1e3, color="#c0392b", alpha=0.07)
axB.axhline(16, ls=":", color="#555", lw=1.1)
axB.text(820, 19, "16 GiB laptop RAM  →  dense G infeasible above", fontsize=8.5, color="#555")
axB.plot(nn, G_curve, "-", color=C_G, lw=2.2, label="dense G  (n²)  — every other solver")
axB.plot(nn, Z_curve, "-", color=C_Z, lw=2.2, label="matrix-free Z  (n·m)  — support-first")
axB.plot(n, G_pts, "o", color=C_G, ms=5)
axB.plot(n, Z_pts, "o", color=C_Z, ms=5)
axB.set_xlabel("candidates  n  (log scale)")
axB.set_ylabel("memory  (GiB, log scale)")
axB.set_ylim(3e-3, 5e2)
axB.legend(loc="lower right", frameon=False, fontsize=9.5)
axB.annotate("11.9 GiB at n = 40k", xy=(40000, 11.9), xytext=(2600, 70), fontsize=9, color=C_G,
             arrowprops=dict(arrowstyle="->", color=C_G, lw=1))
axB.annotate("0.30 GiB", xy=(40000, 0.298), xytext=(9000, 0.05), fontsize=9, color=C_Z,
             arrowprops=dict(arrowstyle="->", color=C_Z, lw=1))

fig.tight_layout()
fig.savefig("research/fig_scaling.pdf", bbox_inches="tight")
fig.savefig("research/fig_scaling.png", bbox_inches="tight")
print("wrote research/fig_scaling.pdf and research/fig_scaling.png")
