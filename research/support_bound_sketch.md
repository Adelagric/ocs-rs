# Toward a bound on the OCS support size

> Working note — a *route* to turning the empirical observation (|S| small and
> bounded in n) into a theorem, with the subtleties stated honestly. Not a finished
> proof. Intended as the seed for the theory follow-on and the INRAE / optimisation
> collaboration. References at the end were verified in the prior-art search.

## The object

Sexed OCS is the second-order cone program

  maximise **bᵀc**  s.t.  **A c = d**  (p rows: p = 1 simplex, p = 2 sexed),
  **c ≥ 0**,  **cᵀGc ≤ k**,   with **G = ZZᵀ/s + εI**, **Z** ∈ ℝ^{n×m}.

Let S = supp(c\*) = {i : c\*ᵢ > 0} at the optimum. We observe (Results, and
independently Waldmann 2025): |S| is small (≈15 at a working cap), **bounded as n
grows** (14–19 up to n = 40000 at a fixed cap), and grows monotonically as the cap
k is tightened toward zero coancestry.

## What is clean

**Loose cap (kinship slack at the optimum).** If cᵀGc\* < k, the active set near
c\* is just {A c = d, c ≥ 0}: a polytope, and c\* is a basic feasible solution — an
LP vertex of a system with p equality rows. Hence **|S| ≤ p**: the unconstrained
sexed optimum puts all mass on the best male and the best female (|S| = 2). This is
exactly the loose end of the frontier.

## Where the subtlety is (and where a naive argument fails)

**Binding cap (cᵀGc\* = k).** Now c\* lies on the *curved* boundary of the
ellipsoid. A naive constraint count — |S| free variables, minus p affine equalities,
minus 1 active quadratic ⇒ |S| = p + 1 — is **wrong**: it is the count for *linear*
active constraints (LP vertices). The quadratic boundary is smooth and strictly
curved, so every point on it (intersected with a coordinate face) is already an
extreme point of the feasible convex set; the curvature, not a 0-dimensional vertex,
makes it extreme. So the support of an extreme point of a set with one quadratic
constraint is **not** bounded by p + 1 — which is consistent with observing |S| ≈ 15,
not 3. Any honest statement must avoid this trap.

## The route to n-independence (the part that does hold)

The leverage is the **low rank** of the genomic matrix. Write the stationarity on
the support, with multipliers μ ∈ ℝ^p (equalities) and λ > 0 (kinship):

  for i ∈ S:  bᵢ − (Aᵀμ)ᵢ − 2λ(Gc\*)ᵢ = 0,   and  (Gc\*)ᵢ = ε c\*ᵢ + (Z g)ᵢ / s,
  where  **g := Zᵀc\*** ∈ ℝ^m.

So c\*ᵢ = [ bᵢ − (Aᵀμ)ᵢ − 2λ (Zᵢ·g)/s ] / (2λε) on S. The whole optimum is therefore
parameterised by **(μ, λ, g) ∈ ℝ^{p+1+m}** — a description whose size depends on the
**number of markers m, not on n**. Equivalently, the kinship constraint acts only
through the m-dimensional image g = Zᵀc; an extreme optimum lives on a face of the
feasible set whose dimension is controlled by rank(G) ≤ m (Barvinok–Pataki
extreme-point counting on the lifted cone, where the SOC has dimension m+1). In the
regime **m < n** — many candidates, a modest marker panel, precisely the scaling
experiment (m = 1000, n up to 40000) — this gives a genuine **n-independent** bound:
|S| is controlled by m + p, not by n. That is the rigorous content behind "the
support does not grow with n," and it is exactly the regime where the matrix-free
solver and the empirical plateau live.

Two honest caveats, both to be discharged in the full proof:
1. **The bound is loose.** rank-counting gives |S| ≲ m + p; we *observe* |S| ≈ 15 ≪ m.
   The tight, typical-case size is a separate question — likely tied to the effective
   number of contributors and the genetics below — and is the real prize.
2. **It is a statement about an *extreme* optimum.** Active-set / vertex-returning
   solvers (support-first, the critical line algorithm) return such a point; conic
   interior-point and ADMM solvers return non-sparse interior points and threshold
   post hoc (Waldmann 2025 truncates at 1e-4). So the theorem is about *precisely the
   object support-first computes* — which is the natural home for it.

## The growth half (cap → 0)

As k → 0 the optimum approaches the minimum-coancestry point, which must spread mass
over much of the population to drive cᵀGc down. The genetics quantifies this exactly:
the rate of inbreeding obeys ΔF ∝ E[Σ cᵢ²] (Wray & Thompson 1990; Woolliams & Bijma
2000), and Ne ≈ 1/(2ΔF). Tightening the coancestry cap forces Σcᵢ² down, i.e. forces
contributions to spread — so |S| grows, monotonically, as observed. The optimisation
side bounds |S| from above (n-independently); the genetics side explains its growth
along the frontier. A clean theorem should marry the two: |S| as a function of the
target ΔF / Ne, capped by the rank term.

## Status

- **Reachable, not yet proved.** The n-independence via rank is a short, standard
  extreme-point argument once the lifted cone is set up carefully; the tight
  typical-case |S| is open and is the publishable core.
- **Most promising single lemma:** an extreme optimum of the OCS SOCP has support
  bounded by p + (active SOC face dimension) ≤ p + rank(G), independent of n; sharpen
  the face term toward the observed constant using the structure of G_SS on the
  support.
- **Validation in hand:** the support-vs-frontier table (Results) and Waldmann (2025)
  give the numbers any bound must reproduce.

## References (verified)

- Pataki, G. (1998). On the rank of extreme matrices in semidefinite programs…
  *Math. Oper. Res.* 23(2):339–358. DOI 10.1287/moor.23.2.339. — extreme-point rank counting.
- Markowitz, H. (1956). *Naval Res. Logist. Q.* 3:111–133. DOI 10.1002/nav.3800030110. — the critical-line / corner-portfolio structure (the "grows along the frontier" analogue).
- Wray, N.R. & Thompson, R. (1990). *Genet. Res.* 55(1):41–54. DOI 10.1017/S0016672300025180. — ΔF ∝ Σ(contribution²).
- Woolliams, J.A. & Bijma, P. (2000). *Genetics* 154(4):1851–1864. DOI 10.1093/genetics/154.4.1851. — contributions ↔ ΔF ↔ Ne.
- Yamashita, M., Mullin, T.J. & Safarina, S. (2018). *Optim. Lett.* 12(7):1683–1697. DOI 10.1007/s11590-018-1229-y. — OCS as a single-SOC program (the cone setup).
- Waldmann, P. (2025). *Bioinform. Adv.* 5(1):vbaf259. DOI 10.1093/bioadv/vbaf259. — independent empirical support counts; post-hoc truncation.
