"""The PLINK `.bed` reader, checked against an encoder written from the format spec.

The Rust test suite encodes its own fixtures; this one re-implements the packing in a
different language, so both decoders would have to be wrong in the same way for the
two suites to agree. The panel then goes straight into a solve — the path a user with
real data actually takes.
"""

from pathlib import Path

import numpy as np

import ocs_rs

# PLINK 1 two-bit codes, keyed by A1 dosage; None is a missing call.
CODE = {2: 0b00, 1: 0b10, 0: 0b11, None: 0b01}


def write_trio(prefix: Path, calls) -> None:
    """Write `<prefix>.bed/.bim/.fam` for `calls[j][i]` = genotype of i at marker j."""
    n = len(calls[0])
    data = bytearray([0x6C, 0x1B, 0x01])  # magic + SNP-major
    for column in calls:
        row = bytearray((n + 3) // 4)
        for i, call in enumerate(column):
            row[i // 4] |= CODE[call] << (2 * (i % 4))
        data += row
    Path(f"{prefix}.bed").write_bytes(bytes(data))
    Path(f"{prefix}.fam").write_text("".join(f"FAM{i} IND{i} 0 0 1 -9\n" for i in range(n)))
    Path(f"{prefix}.bim").write_text(
        "".join(f"1 rs{j} 0 {(j + 1) * 1000} A G\n" for j in range(len(calls)))
    )


def random_calls(n, m, seed=0, missing_rate=0.05):
    rng = np.random.default_rng(seed)
    p = rng.uniform(0.05, 0.5, m)
    dosages = rng.binomial(2, p[:, None], size=(m, n))
    missing = rng.random((m, n)) < missing_rate
    return [
        [None if missing[j, i] else int(dosages[j, i]) for i in range(n)] for j in range(m)
    ]


def expected_panel(calls):
    """Centred genotypes and VanRaden scale, computed independently with numpy."""
    raw = np.array(
        [[np.nan if g is None else float(g) for g in column] for column in calls]
    ).T  # (n, m)
    p = np.nanmean(raw, axis=0) / 2.0
    z = np.nan_to_num(raw - 2.0 * p, nan=0.0)
    return z, p, float(2.0 * np.sum(p * (1.0 - p)))


def test_reads_genotypes_frequencies_and_scale(tmp_path):
    # n = 37 is not a multiple of 4, so every row ends in padding bits the reader
    # must ignore; some calls are missing.
    calls = random_calls(n=37, m=11, seed=7)
    prefix = tmp_path / "panel"
    write_trio(prefix, calls)

    panel = ocs_rs.read_plink(prefix)
    z, p, s = expected_panel(calls)

    assert (panel.n, panel.m) == (37, 11)
    assert panel.ids[:2] == ["IND0", "IND1"]
    assert np.allclose(panel.p, p, atol=1e-12)
    assert np.allclose(panel.Z, z, atol=1e-12)
    assert abs(panel.s - s) < 1e-12


def test_panel_feeds_a_solve(tmp_path):
    calls = random_calls(n=120, m=200, seed=3, missing_rate=0.02)
    prefix = tmp_path / "panel"
    write_trio(prefix, calls)

    panel = ocs_rs.read_plink(prefix)
    rng = np.random.default_rng(0)
    b = panel.Z @ rng.standard_normal(panel.m)

    ridge = 1e-5
    g = panel.Z @ panel.Z.T / panel.s + ridge * np.eye(panel.n)
    one = np.ones(panel.n)
    kmin = float(one @ g @ one) / panel.n**2
    k = kmin + 0.3 * (float(np.diag(g).mean()) - kmin)

    res = ocs_rs.solve(panel.Z, b, k=k, s=panel.s, ridge=ridge)

    assert res.status == "Solved"
    assert abs(res.c.sum() - 1.0) < 1e-9
    assert res.quad <= k * (1 + 1e-9)
    assert 0 < len(res.support) < panel.n


def test_missing_trio_is_reported(tmp_path):
    try:
        ocs_rs.read_plink(tmp_path / "absent")
    except ValueError as exc:
        assert "absent" in str(exc)
    else:
        raise AssertionError("reading an absent trio should raise")
