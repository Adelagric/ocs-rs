#!/usr/bin/env python3
"""Build AlphaMate input files from an OCS instance.

AlphaMate (Gorjanc & Hickey 2018) is driven by a specification file plus three
plain-text data files, all keyed by an individual identifier in the first
column. This script converts the three native arrays of an OCS instance into
exactly those files:

    relationship matrix  (n x n)  ->  Nrm.txt        (id + n values per row)
    EBV / merit vector   (n)      ->  Criterion.txt  (id + 1 value per row)
    sex vector           (n)      ->  Gender.txt     (id + {1=male, 2=female})

and emits a ready-to-run AlphaMateSpec.txt.

Mapping note (relationship vs coancestry)
-----------------------------------------
The OCS spike enforces  c' G c <= k  with G a VanRaden genomic relationship
matrix and contributions c summing to 1. AlphaMate's "group coancestry" of a
contribution vector is  C = 1/2 c' A c  for the matrix A it is given. Feeding G
through `NrmMatrixFile` (A := G) therefore makes AlphaMate's group coancestry
equal to  1/2 c' G c , so the solver's bound  c' G c <= k  is the AlphaMate
target  TargetCoancestry = k / 2 . Pass --target-coancestry with the HALVED
value, or use --target-degree for the angle-based frontier target the manual's
example uses.

Input file formats accepted by this script
------------------------------------------
--matrix : whitespace- or comma-delimited square numeric matrix, n rows x n
           cols, NO id column and NO header (a bare Gram/relationship matrix).
           If your matrix already has an id column, pass --matrix-has-ids.
--ebv    : one value per line (length n), OR two columns "id value".
--sex    : one token per line (length n), OR two columns "id sex". Sex tokens
           may be 1/2, m/f, male/female, M/F (case-insensitive).

If no ids are supplied anywhere, ids default to 1..n (as plain integers), which
AlphaMate accepts.

Usage
-----
    python3 make_alphamate_inputs.py \
        --matrix G.txt --ebv b.txt --sex sex.txt \
        --outdir run1 \
        --matings 50 --male-parents 25 --female-parents 50 \
        --target-degree 30

This writes run1/{Nrm.txt,Criterion.txt,Gender.txt,AlphaMateSpec.txt}. Then:

    cd run1 && /path/to/AlphaMate        # reads AlphaMateSpec.txt by default
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _split(line: str) -> list[str]:
    line = line.strip()
    if not line:
        return []
    if "," in line:
        return [t.strip() for t in line.split(",") if t.strip() != ""]
    return line.split()


def read_matrix(path: Path, has_ids: bool) -> tuple[list[str], list[list[float]]]:
    rows: list[list[float]] = []
    ids: list[str] = []
    for raw in path.read_text().splitlines():
        toks = _split(raw)
        if not toks:
            continue
        if has_ids:
            ids.append(toks[0])
            vals = toks[1:]
        else:
            vals = toks
        rows.append([float(v) for v in vals])
    n = len(rows)
    if n == 0:
        sys.exit(f"error: matrix {path} is empty")
    for i, r in enumerate(rows):
        if len(r) != n:
            sys.exit(
                f"error: matrix {path} is not square: row {i} has {len(r)} "
                f"values but there are {n} rows"
            )
    if not ids:
        ids = [str(i + 1) for i in range(n)]
    return ids, rows


def read_vector(path: Path) -> tuple[list[str] | None, list[str]]:
    ids: list[str] = []
    vals: list[str] = []
    have_ids = False
    for raw in path.read_text().splitlines():
        toks = _split(raw)
        if not toks:
            continue
        if len(toks) == 1:
            vals.append(toks[0])
        elif len(toks) == 2:
            have_ids = True
            ids.append(toks[0])
            vals.append(toks[1])
        else:
            sys.exit(f"error: {path} has a row with {len(toks)} columns (expected 1 or 2)")
    return (ids if have_ids else None), vals


_SEX_MAP = {
    "1": "1", "m": "1", "male": "1",
    "2": "2", "f": "2", "female": "2",
}


def normalise_sex(tokens: list[str]) -> list[str]:
    out = []
    for t in tokens:
        key = t.strip().lower()
        if key not in _SEX_MAP:
            sys.exit(
                f"error: unrecognised sex code {t!r}; use 1/2, m/f, or male/female"
            )
        out.append(_SEX_MAP[key])
    return out


def fmt(x: float) -> str:
    # AlphaMate reads free-format reals; keep full double precision.
    return repr(float(x))


def main() -> None:
    ap = argparse.ArgumentParser(description="Build AlphaMate inputs from an OCS instance.")
    ap.add_argument("--matrix", required=True, type=Path,
                    help="square relationship/coancestry matrix (n x n)")
    ap.add_argument("--matrix-has-ids", action="store_true",
                    help="first column of --matrix is an id column")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="multiply every matrix entry by this factor before writing. "
                         "AlphaMate expects a relationship matrix (self ~ 1+F); pass "
                         "--scale 2 to turn a coancestry/kinship matrix sKin into the "
                         "relationship matrix G = 2*sKin (diagonal ~1), otherwise its "
                         "criterion-maximisation mode reads self-coancestry ~0.5 as an "
                         "inbreeding of -0.5 and faults.")
    ap.add_argument("--as-coancestry", action="store_true",
                    help="write CoancestryMatrixFile instead of NrmMatrixFile "
                         "(use when the matrix is already a coancestry/kinship matrix)")
    ap.add_argument("--ebv", required=True, type=Path,
                    help="EBV / selection-criterion vector (length n, or 'id value')")
    ap.add_argument("--sex", required=True, type=Path,
                    help="sex vector (length n, or 'id sex'); 1=male 2=female")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=15)
    ap.add_argument("--matings", type=int, default=None)
    ap.add_argument("--male-parents", type=int, default=None)
    ap.add_argument("--female-parents", type=int, default=None)
    ap.add_argument("--target-degree", type=float, default=None,
                    help="angle target (deg) on the gain-vs-coancestry frontier; "
                         "0=max coancestry-control, 90=max gain")
    ap.add_argument("--target-coancestry", type=float, default=None,
                    help="absolute group-coancestry target (= k/2 for a c'Gc<=k bound)")
    ap.add_argument("--target-coancestry-rate", type=float, default=None,
                    help="target rate of group coancestry (delta F) vs current population")
    ap.add_argument("--equalize", action="store_true",
                    help="add EqualizeMaleContributions/EqualizeFemaleContributions = Yes")
    args = ap.parse_args()

    ids, mat = read_matrix(args.matrix, args.matrix_has_ids)
    n = len(ids)

    ebv_ids, ebv = read_vector(args.ebv)
    sex_ids, sex_raw = read_vector(args.sex)
    sex = normalise_sex(sex_raw)

    for name, v in (("ebv", ebv), ("sex", sex)):
        if len(v) != n:
            sys.exit(f"error: {name} has length {len(v)} but matrix has {n} individuals")

    # If the vector files carried their own ids, they must line up with the matrix ids.
    if ebv_ids is not None and ebv_ids != ids:
        sys.exit("error: --ebv ids do not match --matrix ids (same order required)")
    if sex_ids is not None and sex_ids != ids:
        sys.exit("error: --sex ids do not match --matrix ids (same order required)")

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    matrix_name = "Coancestry.txt" if args.as_coancestry else "Nrm.txt"
    sc = args.scale
    with (out / matrix_name).open("w") as fh:
        for i in range(n):
            fh.write(ids[i] + " " + " ".join(fmt(x * sc) for x in mat[i]) + "\n")

    with (out / "Criterion.txt").open("w") as fh:
        for i in range(n):
            fh.write(f"{ids[i]} {fmt(float(ebv[i]))}\n")

    with (out / "Gender.txt").open("w") as fh:
        for i in range(n):
            fh.write(f"{ids[i]} {sex[i]}\n")

    n_male = sum(1 for s in sex if s == "1")
    n_female = sum(1 for s in sex if s == "2")

    spec: list[str] = []
    spec.append(f"Seed                        , {args.seed}")
    spec.append("GenderFile                  , Gender.txt")
    key = "CoancestryMatrixFile" if args.as_coancestry else "NrmMatrixFile"
    spec.append(f"{key:<28}, {matrix_name}")
    spec.append("SelCriterionFile            , Criterion.txt")
    matings = args.matings if args.matings is not None else max(1, n_female)
    spec.append(f"NumberOfMatings             , {matings}")
    mp = args.male_parents if args.male_parents is not None else max(1, n_male)
    spec.append(f"NumberOfMaleParents         , {mp}")
    if args.equalize:
        spec.append("EqualizeMaleContributions   , Yes")
    fp = args.female_parents if args.female_parents is not None else max(1, n_female)
    spec.append(f"NumberOfFemaleParents       , {fp}")
    if args.equalize:
        spec.append("EqualizeFemaleContributions , Yes")
    if args.target_degree is not None:
        spec.append(f"TargetDegree                , {args.target_degree:g}")
    if args.target_coancestry is not None:
        spec.append(f"TargetCoancestry            , {args.target_coancestry:g}")
    if args.target_coancestry_rate is not None:
        spec.append(f"TargetCoancestryRate        , {args.target_coancestry_rate:g}")
    spec.append("Stop")

    (out / "AlphaMateSpec.txt").write_text("\n".join(spec) + "\n")

    print(f"wrote {out}/{matrix_name}  ({n} x {n})")
    print(f"wrote {out}/Criterion.txt ({n} rows)")
    print(f"wrote {out}/Gender.txt     ({n} rows; {n_male} male, {n_female} female)")
    print(f"wrote {out}/AlphaMateSpec.txt")
    print()
    print("spec:")
    print("\n".join("  " + line for line in spec))
    print()
    print(f"run with:  cd {out} && /path/to/AlphaMate")


if __name__ == "__main__":
    main()
