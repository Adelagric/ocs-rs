# Preprint skeleton — support-first OCS

Target: **bioRxiv** preprint → journal in animal/plant breeding genetics (GSE or G3).
Framing: a **genetics methods** paper. The optimisation engine is a rigorous Methods
section, not the headline; the headline is *exact OCS at genomic scale, orders of
magnitude faster than the tools breeders actually use, validated on real data*.
This framing is robust to the prior-art outcome (see §Novelty): its value does not
depend on the abstract algorithm being new, but on the genomic-scale matrix-free
realisation and the empirical demonstration.

Author: Adel Kaleche. Affiliation: TBD (lycée agricole; INRAE collaboration being explored).
Single-author draft.

---

## Working title (pick one)

1. **Exact optimum contribution selection in milliseconds: a matrix-free, support-first
   solver validated against optiSel and AlphaMate on real genomic data**
2. Support-first: exact optimum contribution selection at genomic scale by active-set
   optimisation on a tiny solution support
3. Orders-of-magnitude faster exact optimum contribution selection by exploiting the
   sparsity of the optimal contribution vector

(Lead with #1 — it states the contribution and the evidence in one line.)

---

## Abstract (draft — numbers locked from the benchmark tables)

Optimum contribution selection (OCS) — maximise genetic gain subject to a cap on the
coancestry of the next generation — is the central tool for balancing gain against the loss
of diversity in breeding programmes. Its computational bottleneck is the relationship matrix
it constrains: at genomic scale that matrix is dense and n×n, costing O(n²m) to build and
O(n²) to store, so the established exact (optiSel) and heuristic (AlphaMate) tools run into a
memory wall exactly as candidate numbers grow. We present **support-first**, an exact solver
that never forms it: on a laptop it solves an instance (n = 40000) whose dense matrix
(11.9 GiB) the established tools cannot allocate, the optimal support holding at about fifteen
individuals — doing what is otherwise infeasible, not merely doing it faster. It exploits two
structural facts of OCS: the optimal
contribution vector is supported on a tiny subset of candidates (typically 2–50 of n),
and the only difficulty is combinatorial — which candidates are free versus pinned at a
bound (zero, or a per-candidate cap). Support-first grows the
support by an active-set / column-generation rule and solves each fixed-support subproblem
— maximise a linear objective over an ellipsoid intersected with the affine sum (and sex)
constraints — in closed form, with no inner iterative solver; the kinship products G·c are
formed matrix-free from the genotype matrix Z, the enabler that lets it never store the dense
n×n G (the order-of-magnitude speed is the tiny active set, not the matrix-free product). On real data — a CIMMYT
wheat panel (n=599), a PIC pig (n=3534, 52k SNP), and a heterogeneous-stock mouse panel
(n=1814, with real sex) — support-first returns the exact optimum on the kinship
boundary — agreeing with a conic interior-point solver to 1e-8, which itself stops just short
of that boundary — while running 90×–2280× faster (0.008 s vs 6.96 s on the sexed mouse
instance) and ~37000× faster than a general conic solver (Clarabel) at n=10000. Against
AlphaMate — a heuristic for the distinct problem of discrete mate allocation — the exact
optimum is no worse at matched coancestry on the continuous relaxation the two share, at a
fraction of the run time; per-candidate contribution caps
(0 ≤ c ≤ u) are supported. Support-first makes exact, reproducible OCS practical at genomic
scale on a laptop.

---

## Section structure (IMRaD)

### 1. Introduction
- OCS: definition, why coancestry control matters (inbreeding depression, long-term gain). Refs: Meuwissen 1997; Woolliams et al. review.
- The genomic-scale bottleneck: dense G, QCQP cost; what breeders run today.
- The existing tools and what each is: GENCONT (Lagrangian), optiSel/Wellmann (cone QP via cccp), AlphaMate/Gorjanc & Hickey (evolutionary heuristic), generic conic solvers.
- The gap: no exact solver that exploits the structure of OCS to reach genomic scale cheaply.
- Contribution bullets (state plainly; defer novelty wording to §Novelty verdict).

### 2. Methods
- **2.1 OCS as a QCQP/SOCP.** maximise bᵀc s.t. 𝟙ᵀc=1, c≥0, cᵀGc≤k. Sexed variant: Σ_males c = Σ_females c = ½. G = VanRaden ZZᵀ/s + εI.
- **2.2 Structure.** (i) optimal support is tiny and bounded; (ii) with the support fixed and c≥0 inactive, the KKT system is a max-linear-over-ellipsoid-∩-affine problem.
- **2.3 Closed form per support.** Null-space reduction of the sum/sex constraints → scalar quadratic in the multiplier μ (Vieta); recovery of (μ, λ) by least squares. No inner solver.
- **2.4 Active set / column generation.** Reduced-cost entry rule; drop on negativity; chunked feasibility phase; anti-cycling (taboo list, relative degeneracy threshold). Exactness/termination argument.
- **2.5 Matrix-free G·c.** G·c = ε·c + Z(Zᵀc)/s; per-product cost O(n·m), memory O(n·m) for Z vs O(n²) for a dense G. Honest tradeoff (measured, n=2000/m=10000): when m > n the Z(Zᵀc) product moves *more* data than a resident dense G·c, so matrix-free is **not** an inner-loop speedup there — it is the memory enabler (never allocate n×n) and wins on time only when n ≫ m. The order-of-magnitude speedup over optiSel/Clarabel is **algorithmic** (tiny active set), independent of this choice.
- **2.6 Implementation.** Rust; pure-Rust dependency stack (faer dense linalg, clarabel as an independent cross-check oracle); zero unsafe; seeded, reproducible. Single-binary CLI.
- **2.7 Data & baselines.** Datasets (provenance, n, m, trait used as b, sex availability). Baselines: optiSel (exact reference), Clarabel (conic cross-check), AlphaMate (heuristic). Hardware + exact commands (BENCHES.md discipline).

### 3. Results
- **3.1 Exactness.** support-first vs optiSel and vs Clarabel: max |Δc|, Δgain, Δcoancestry across k. → Δgain ≈ 1e-8 vs the conic optimum; cross-language vs the numpy prototype 1.5e-14; equal gain at matched coancestry, with support-first on the boundary where optiSel's IPM stops interior.
- **3.2 Speed & scaling.** vs Clarabel to n=10000 (37090×: 1579 s → 0.043 s); vs optiSel on real data (90×–2280×).
- **3.3 Real-data benchmark table** (the headline — see SUPPORT_FIRST.md for current numbers).
- **3.4 vs AlphaMate.** Equal-coancestry comparison: take AlphaMate's contribution vector, score gain and cᵀKc in our metric, compare support-first at the same coancestry; report time (AlphaMate self-reported CPU vs support-first wall). Result (mouse): on the continuous relaxation the two share, these Δgain are the optimality gap the differential-evolution heuristic leaves on the table — +0.004 at the 45° tradeoff, +0.018 / +0.080 at the higher / lower coancestry corners — read as an exact-vs-stochastic consistency check, not a win at AlphaMate's own discrete-mating task (out of scope here). AlphaMate's frontier cost 882 s CPU vs ≤1.1 s per point exact. Also report robustness: AlphaMate required 6 configurations and 3 distinct work-arounds (matings<n; full parent set to avoid a setup segfault; positive-shifted criterion to avoid a value/max sign inversion) to run on a real genomic instance, whereas support-first and optiSel ran unmodified.
- **3.5 Support behaviour.** |support| vs k: tiny at the operating point (19 at the mouse working coancestry), grows only as coancestry → 0 (~1163 to force group coancestry ≈ 0); bounded 14–19 as n→40000 at a fixed binding cap. Interpretable, exact frontier traced in ms–s.

### 4. Discussion
- Practical significance: OCS on a laptop at genomic scale; exact and reproducible vs stochastic.
- **Limitations (state up front):** b is a heritable-phenotype *proxy*, not a true GEBV, on these public panels; sex is real only for the mouse panel (arbitrary 50/50 elsewhere); cross-language timing (Python prototype for the sexed benchmark vs R optiSel — the Rust binary is faster still, but the head-to-head language differs and is labelled); single-trait, single quadratic constraint (no multiple kinship/own-relationship constraints yet).
- Novelty positioning → §Novelty verdict.
- Future work: multiple constraints, integer/mate-allocation layer, true-GEBV datasets, larger n.

### 5. Code & data availability
- Source repository (Rust crate + benchmark scripts + this research note), with the exact public datasets and the commands to reproduce every number.

---

## Novelty positioning (prior-art verdict)

Adversarial prior-art search done (OCS + portfolio-optimization literature). The
method bundles four ingredients that **do not share a fate** — be precise, or a
reviewer separates them for us:

**Concede openly (each is standard, in isolation):**
- Active-set on non-negativity exploiting a small support → **Markowitz Critical
  Line Algorithm (1956)**; Stein–Branke–Schmeck (2008) for large-scale long-only.
- Closed-form "max linear over ellipsoid ∩ affine via a scalar secular equation"
  → **Gander, Golub & von Matt (1989)**; Moré–Sorensen; Merton (1972).
- Matrix-free `G·v` via `Z(Zᵀv)` → **VanRaden 2008 / Legarra–Misztal 2008** (GBLUP).

**Defensible novelty (survives, stated narrowly):** the *combination, specialised
to the OCS quadratic-kinship QCQP* of (i) **reduced-cost / column-generation
support-growing toward a single fixed coancestry cap k** and (ii) **genomic
matrix-free `G·c` evaluated only over the sparse support**. No verified OCS work
does either: in-domain solvers form/ingest a dense relationship matrix and either
truncate post-hoc (optiSel, AlphaMate, Clark, **Waldmann 2025** — newest, forms
dense `G`, truncates ρ<1e-4) or exploit *pedigree-matrix* sparsity inside a
full-population interior-point solve (**Yamashita et al. 2018**).

**Must pre-empt by name (or look unaware of the closest neighbours):**
- **Yamashita, Mullin & Safarina 2018** (*Optimization Letters*, arXiv:1506.04487)
  — OCS + sparsity + SOCP, but it is **A⁻¹ pedigree-matrix** sparsity (a data
  property) inside a full-candidate ECOS solve, **not** solution-support, **not**
  genotype-matrix-free. State our distinction explicitly and up front.
- **Markowitz CLA 1956** — frame our technique as *adapted from* portfolio
  active-set methods. The delta: CLA does **parametric continuation that sweeps λ
  to trace the whole frontier**; support-first does **reduced-cost pricing toward a
  single target k**. (Verified against our code: the loop adds the best positive
  reduced-cost candidate at fixed k — it is pricing, not a λ-sweep. The novelty
  framing must rest on this, and it holds.)
- **Gander–Golub–von Matt 1989** — cite the per-support closed form as known.

**Biggest risk:** a reviewer saying *"this is Markowitz CLA on OCS with a known
matrix-free trick bolted on."* Defense = the pricing-toward-fixed-k mechanism +
that it is exactly what makes the genomic matrix-free product pay (cost scales
with the sparse support × markers, never n²). Source caveats: a few method bodies
(Pong-Wong & Woolliams 2007, Schierenbeck 2011) were paywalled — SDP nature from
abstracts; "not found in OCS" = thoroughly searched English-indexed literature.

## Evidence ledger — verified / measured / extrapolated / missing

| Claim | Status | Source / what's needed |
|---|---|---|
| support-first returns the OCS optimum exactly | **Verified** | matches optiSel & Clarabel to 1e-9–1e-12 across k (correctness tests + sf_at_ub reproduces optiSel gain/coancestry) |
| Optimal support is tiny (2–50) and grows only near zero coancestry | **Measured** | support 19 at mouse operating point; ~1163 forcing coancestry≈0 (sf_at_ub run) |
| Optimal support is **bounded in n** | **Measured** | scaling_matrixfree: support 14–19 across n = 1000→40000 at a fixed binding cap; the basis for cost scaling with the support, not n |
| Matrix-free G·c, never forms G | **Verified + measured** | per-product O(n·m); measured slower than dense G·c when m>n (n=2000/m=10000) — it is the memory/large-n enabler, not an inner-loop speedup. Scaling plot done (scaling_matrixfree): dense G → 11.9 GiB infeasible at n=40k while matrix-free solves in 78 ms |
| Speedup source is algorithmic, not matrix-free | **Measured** | both the dense-K Python proto and the matrix-free Rust beat optiSel by the same order; the win is the tiny active set, not the G·c route. State this explicitly so the matrix-free claim is not oversold |
| ~37000× vs Clarabel | **Measured** | spike benchmark; re-confirm command + hardware in BENCHES.md |
| 22×–2280× vs optiSel on real data | **Measured** | wheat / pig / mouse; pig 0.024 s vs 54.8 s; mouse 0.008 s vs 6.96 s. Caveat: Python proto vs R |
| vs AlphaMate at equal coancestry | **Measured** | mouse run6: support-first gain strictly > AlphaMate at all 3 frontier points (Δ +0.004 at the 45° opt, +0.018 / +0.080 elsewhere); AlphaMate 882 s CPU for the frontier vs ≤1.1 s/point exact |
| AlphaMate fragility on genomic data | **Measured (qualitative)** | 6 configs, 3 work-arounds (matings<n; full parent set vs a setup segfault; positive-shifted criterion vs a value/max sign inversion); tighten into one paragraph |
| Novelty of the method | **Resolved** | see §Novelty positioning — the OCS combination (reduced-cost support-growing + genomic matrix-free) is defensible; 3 ingredients conceded; pre-empt Yamashita 2018 + Markowitz CLA |
| Sexed solver in Rust | **Verified + measured** | `solve_sexed` ported; KKT-certificate + across-k invariant tests; gate green (fmt/clippy -D warnings/9 tests). Cross-language exactness vs the Python prototype on an identical binding instance: **Δgain = 1.5e-14** (machine precision). Release timing (n=2000/m=10000, `examples/bench_sexed`): 1.8 ms (support 2) → 294 ms (support 17) as the cap tightens |
| Reproducibility package | **Partly** | scripts exist under research/; needs a one-command reproduce + pinned dataset fetch |

## Hard gaps to close before submission (priority order)
1. ~~Novelty verdict~~ — **DONE**, see §Novelty positioning. Framing locked: genetics paper, combination-novelty, pre-empt Yamashita + CLA.
2. ~~AlphaMate equal-coancestry point~~ — **DONE** (mouse run6): exact dominates at every frontier point; 882 s vs ≤1.1 s.
3. ~~Rust sexed solver~~ — **DONE**: ported, tested, gate green, cross-language exact (Δgain 1.5e-14), release timing measured. Surfaced an honest reframe (below) now baked into the abstract/§2.5.
4. ~~Matrix-free scaling figure~~ — **DONE** (`examples/scaling_matrixfree`, m=1000, binding k=0.1·diag): support stays **14–19 as n goes 1000→40000** (the "support bounded in n" claim, confirmed); matrix-free solve scales ~linearly (8→78 ms); dense G grows O(n²) to **11.9 GiB at n=40000 (infeasible to allocate)** with O(n²·m) build. Turn into the headline figure (two panels: support & solve-time bounded vs G-memory exploding). **Reframe locked:** speedup over optiSel/Clarabel is algorithmic (tiny active set); matrix-free is the memory/large-n enabler, never sold as inner-loop speed.
5. **Reproducibility**: one command that fetches the public datasets and regenerates every table/figure.
6. **AlphaMate cosmetics for the paper**: it was run via a Linux x86-64 binary under Rosetta (no macOS build; MKL-locked source) — report its self-timed CPU, label the emulation, don't quote wall-time as native.
