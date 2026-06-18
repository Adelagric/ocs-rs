"""Headline figure: support-first scales to genomic n where the dense matrix cannot.

Data from examples/scaling_matrixfree (m=1000, binding cap k=0.1*mean_diag,
release profile). Memory curves are exact footprints (n^2*8 for dense G,
n*m*8 for Z). Build time and solve time are measured; the dense-G build at
n=40000 is projected (G too large to build inside the free RAM) and marked open.
Outputs research/fig_scaling.{pdf,png}.
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- measured (binding cap), examples/scaling_matrixfree ---
n = np.array([1000, 2000, 5000, 10000, 20000, 30000, 40000], float)
support = np.array([14, 15, 15, 19, 16, 17, 16], float)
solve_s = np.array([0.0079, 0.0085, 0.0158, 0.0307, 0.0413, 0.0573, 0.0785])
# dense-G build time: measured to n=30000; n=40000 not built (RAM), projected ~n^2.
build_n = np.array([1000, 2000, 5000, 10000, 20000, 30000], float)
build_s = np.array([0.006, 0.020, 0.102, 0.402, 1.568, 3.599])
build_proj_n, build_proj_s = 40000.0, 3.599 * (40000 / 30000) ** 2
m = 1000

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
axA.plot(n, solve_s, "s-", color=C_SOLVE, lw=2.2, ms=5.5, label="support-first: full solve")
axA.set_xlabel("candidates  n  (log scale)")
axA.set_ylabel("time  (s, log scale)")
axA.set_ylim(3e-3, 2e1)
axA.legend(loc="upper left", frameon=False, fontsize=9.5)
axA.annotate(
    "63× the solve at n=30k\n(3.6 s vs 57 ms)",
    xy=(30000, 3.6), xytext=(1500, 6.0), fontsize=9, color=C_BUILD,
    arrowprops=dict(arrowstyle="->", color=C_BUILD, lw=1),
)
axA.text(1150, 5.5e-3, "support |S| = 14–19 throughout", fontsize=9, color=C_SOLVE)

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
