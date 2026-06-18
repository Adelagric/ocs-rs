"""Evaluate an AlphaMate contribution vector in the OCS spike's own metric.

AlphaMate reports its solutions in its internal coancestry convention; for a fair
head-to-head we ignore those numbers and re-score AlphaMate's *contribution
vector* with the same b (selection criterion) and K (= sKin, the mean-kinship
matrix) used by the support-first / optiSel benchmark, so every method is judged
by identical definitions:

    gain        = b . c
    coancestry  = c' K c          (K = bench_K_{n}.csv = sKin = G/2 + ridge)

A ContributorsMode*.txt file has columns:
    Id  Gender  SelCriterion  AvgCoancestryA  AvgCoancestryC  Contribution  nContribution
Ids are 1..n here (we emitted them that way), so c[id-1] = Contribution.

Usage:  python3 eval_contrib.py <n> <ContributorsFile> [<ContributorsFile> ...]
Prints, for each file and for the optiSel reference (column oc of bench_bc),
the support size, sum(c), gain and coancestry.
"""

import sys
import numpy as np


def load_instance(n):
    K = np.loadtxt(f"/tmp/bench_K_{n}.csv", delimiter=",")
    rows = [r.split(",") for r in open(f"/tmp/bench_bc_{n}.csv").read().splitlines()[1:]]
    bv = np.array([float(r[0]) for r in rows])
    oc = np.array([float(r[1]) for r in rows])
    sex = np.array([r[2].strip().strip('"') for r in rows])
    return K, bv, oc, sex


def read_contributions(path, n):
    c = np.zeros(n)
    for line in open(path).read().splitlines()[1:]:
        t = line.split()
        if len(t) < 6:
            continue
        idx = int(t[0]) - 1
        c[idx] = float(t[5])
    return c


def score(name, c, bv, K):
    tol = 1e-9
    gain = float(bv @ c)
    coan = float(c @ K @ c)
    supp = int((c > tol).sum())
    print(f"  {name:<34} support={supp:>4}  sum(c)={c.sum():.6f}  "
          f"gain={gain:+.5f}  coancestry(c'Kc)={coan:.6f}")
    return gain, coan


def main():
    n = int(sys.argv[1])
    K, bv, oc, sex = load_instance(n)
    print(f"instance n={n}  (males={int((sex=='male').sum())}, females={int((sex=='female').sum())})")
    score("optiSel (reference, oc)", oc, bv, K)
    for path in sys.argv[2:]:
        name = "AlphaMate:" + path.split("/")[-1].replace("Contributors", "").replace(".txt", "")
        c = read_contributions(path, n)
        score(name, c, bv, K)


if __name__ == "__main__":
    main()
