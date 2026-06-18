# References (verified bibliography)

> Every field below was confirmed against CrossRef, the publisher DOI landing
> page, or Europe PMC. One software entry (faer) has no published DOI and is cited
> as software — see its note.

## OCS methods and tools

1. **Meuwissen, T.H.E. (1997).** Maximizing the response of selection with a
   predefined rate of inbreeding. *Journal of Animal Science* 75(4):934–940.
   DOI 10.2527/1997.754934x.

2. **Dagnachew, B.S. & Meuwissen, T.H.E. (2016).** A fast Newton–Raphson based
   iterative algorithm for large scale optimal contribution selection. *Genetics
   Selection Evolution* 48(1):70. DOI 10.1186/s12711-016-0249-2.

3. **Pong-Wong, R. & Woolliams, J.A. (2007).** Optimisation of contribution of
   candidate parents to maximise genetic gain and restricting inbreeding using
   semidefinite programming. *Genetics Selection Evolution* 39(1):3–25.
   DOI 10.1186/1297-9686-39-1-3.

4. **Wellmann, R. (2019).** Optimum contribution selection for animal breeding and
   conservation: the R package optiSel. *BMC Bioinformatics* 20(1):25.
   DOI 10.1186/s12859-018-2450-5.

5. **Gorjanc, G. & Hickey, J.M. (2018).** AlphaMate: a program for optimizing
   selection, maintenance of diversity and mate allocation in breeding programs.
   *Bioinformatics* 34(19):3408–3411. DOI 10.1093/bioinformatics/bty375.

6. **Waldmann, P. (2025).** Genomic optimum contribution selection and mate
   allocation using JuMP. *Bioinformatics Advances* 5(1):vbaf259.
   DOI 10.1093/bioadv/vbaf259.

7. **Yamashita, M., Mullin, T.J. & Safarina, S. (2018).** An efficient
   second-order cone programming approach for optimal selection in tree breeding.
   *Optimization Letters* 12(7):1683–1697. DOI 10.1007/s11590-018-1229-y;
   arXiv:1506.04487. — Exploits *pedigree-inverse* sparsity, not solution support.

## Genomic relationships and matrix-free products

8. **VanRaden, P.M. (2008).** Efficient methods to compute genomic predictions.
   *Journal of Dairy Science* 91(11):4414–4423. DOI 10.3168/jds.2007-0980.

9. **Legarra, A. & Misztal, I. (2008).** Technical note: Computing strategies in
   genome-wide selection. *Journal of Dairy Science* 91(1):360–366.
   DOI 10.3168/jds.2007-0403.

## Active-set and closed-form lineage (conceded prior art)

10. **Markowitz, H. (1956).** The optimization of a quadratic function subject to
    linear constraints. *Naval Research Logistics Quarterly* 3(1–2):111–133.
    DOI 10.1002/nav.3800030110. — The critical line algorithm.

11. **Gander, W., Golub, G.H. & von Matt, U. (1989).** A constrained eigenvalue
    problem. *Linear Algebra and its Applications* 114–115:815–839.
    DOI 10.1016/0024-3795(89)90494-1. — The per-support closed form (secular equation).

## Software and data

12. **Goulart, P.J. & Chen, Y. (2024).** Clarabel: an interior-point solver for
    conic programs with quadratic objectives. arXiv:2405.12762. (Peer-reviewed
    version: *Mathematical Programming Computation*, 2026,
    DOI 10.1007/s12532-026-00320-7.) — Used only as an independent cross-check oracle.

13. **Quiñones, S.** faer: a linear-algebra library for the Rust programming
    language. Software, version 0.24.0 (pinned in `Cargo.toml`);
    <https://github.com/sarah-quinones/faer-rs>. (JOSS submission under review; no
    published DOI at time of writing.)

14. **Pérez, P. & de los Campos, G. (2014).** Genome-wide regression and
    prediction with the BGLR statistical package. *Genetics* 198(2):483–495.
    DOI 10.1534/genetics.114.164442. — Source of the wheat (CIMMYT) and mouse panels.

15. **Cleveland, M.A., Hickey, J.M. & Forni, S. (2012).** A common dataset for
    genomic analysis of livestock populations. *G3: Genes|Genomes|Genetics*
    2(4):429–435. DOI 10.1534/g3.111.001453. — The PIC pig panel.

---

### Not cited in the prose but held in reserve (from the prior-art search)

- **Stein, M., Branke, J. & Schmeck, H. (2008).** Efficient implementation of an
  active set algorithm for large-scale portfolio selection. *Computers &
  Operations Research* 35(12):3945–3961. DOI 10.1016/j.cor.2007.05.004. — For the
  portfolio-active-set lineage, if a reviewer presses it.
- **Clarkson, K.L. (2010).** Coresets, sparse greedy approximation, and the
  Frank–Wolfe algorithm. *ACM Transactions on Algorithms* 6(4):article 63, 1–30.
  DOI 10.1145/1824777.1824783. — Pre-empts a "merely Frank–Wolfe" attack on
  incremental support growth (rebuttal: FW is approximate, not exact).
