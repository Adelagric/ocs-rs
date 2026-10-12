# Repository and manuscript corrections

## Scope and authorization

Correct the reviewed source, documentation, English and French manuscripts locally. Preserve the continuous OCS model and solver behavior unless a reproduced bug warrants a targeted change. No push, release, archive modification or message to a third party is authorized until the user reviews the concrete result and approves that destination.

## Findings and acceptance criteria

1. Qualify the zero-ridge support bound: the uncapped theorem does not bound total positive support under finite individual upper bounds. State the capped extension in terms of free coordinates and account separately for upper-bound coordinates. Include a counterexample and correct KKT conditions and equality counts in both languages.
2. Align performance claims with the operational Delta-F=1% measurements (structured 40,000 candidates / 1,000 markers: approximately 0.65 seconds). Separate measurements, projections and resource-policy omissions; distinguish given-G prototype comparisons from shipped matrix-free comparisons.
3. Reproduction must identify which scripts produce each current table and figure, run dependencies in order and report failures honestly. Preserve historical measurements with explicit operating-point labels. A scaling figure must never substitute an incompatible data set silently.
4. Provide runnable R/Python examples with status checks, centered-Z inputs and aligned vector requirements. State packed Rust versus dense wrapper memory behavior accurately. Align R package version with the current repository version, declare required runtime versions and add meaningful combined sex/cap wrapper checks and CI test gates.
5. Keep published DOI identifiers and release metadata truthful; document the difference between editing archive metadata and publishing corrected files as a new version. Do not invent a new version DOI.

## Ordered execution

- Inspect current source and script dependency/output paths; obtain an independent critical review of this plan before implementation.
- Correct theory and both manuscripts as one owned edit area.
- Correct reproduction, scaling labels and provenance as a separate owned edit area.
- Correct bindings, examples, metadata and CI without overlapping ownership.
- Integrate root documentation; search for inconsistent claims and obsolete paths.
- Run Rust tests and relevant formatting/lint checks, freshly build R and Python bindings, run their tests and quickstarts, run reproduction regression checks and selected actual reproduction stages.
- Compile both manuscript sources, inspect compiler diagnostics and render the PDFs for visual review.
- Obtain independent specification and correctness review, resolve findings and repeat affected checks.
- Deliver corrected PDFs, a reviewable patch and a concise validation/publication note. Seek approval only for a concrete external publication action.

## Verification limits

A complete large-panel performance rerun may require external data and long computations. Distinguish historical reported measurements from measurements obtained during this correction. Do not claim that changing documentation revalidates every benchmark or that a successful core test suite verifies wrapper packaging.

## Critical analysis

Independent review completed and accepted below. Main risks: mistranslated bound assumptions; a script that produces the wrong operating-point data; unsupported solver superiority claims; wrapper tests using an installed stale package; PDF figures lagging source edits; premature release metadata changes.

### Independent critical review, accepted before implementation

- Use an existential LP-slice argument from any attained optimum; no Slater or duality assumption is needed. Count the rank of the stacked linear constraints, at most q+r+1. Do not claim that the implemented active-set solver necessarily returns the slice vertex.
- With finite caps, bound strictly free coordinates and count positive coordinates at their upper bounds separately. Upper-bound constraints invalidate a bound on total support alone.
- State both upper and lower KKT multipliers. The affine identity requires positive ridge, positive quadratic multiplier and valid KKT multipliers, and only applies to free coordinates. The parameter count is m+q+1.
- The structured scaling cross-check currently uses different ridge values in the two solver formulations; align them and recompute affected results before making exactness claims.
- The displayed Delta-F derived from 1/sum(c_i^2) is a contribution-based proxy, not a measured realized genomic inbreeding rate. Label this distinction in the manuscript and reproduction guide.
- R source metadata is already version 0.4.0; a locally installed 0.1.0 package is stale and must be replaced by an isolated fresh source build for verification. Keep published version/DOI metadata unchanged until the concrete next release is reviewed.

### Evidence recovered during execution

The original local experiment checkout at e6790c99 contains small historical summary
CSVs for scaling, operational comparisons and multigeneration validation. Recover
these unchanged with checksums as reported measurements; keep fresh corrected-run
outputs distinct. Historical given-G and R columns were independently calibrated
to target contribution proxies, not exported from a common cap. Wheat has material
gain differences at requested caps. The structured uniform coancestry approaches
zero through panel centering; it is not a positive structure floor.

The loose-cap support bound also needs the quadratic constraint to be redundant
on the budget polytope (or a feasible optimal LP vertex), not merely slack at one
optimum. A participation-ratio bound does not imply a support-cardinality bound.

Verification so far: core test suite passes all 45 tests on the correction branch.
