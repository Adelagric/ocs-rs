//! PLINK 1 binary genotype input (`.bed` / `.bim` / `.fam`).
//!
//! Reads a marker panel from the format breeding programmes actually store it in,
//! straight into the centred genotype matrix `Z` the matrix-free solver consumes —
//! so a real dataset reaches an OCS optimum without the dense `n×n` relationship
//! matrix ever existing.
//!
//! Scope and conventions:
//!
//! - Only SNP-major `.bed` (mode byte `0x01`) is read — what PLINK has written by
//!   default since 1.9. Individual-major files are rejected with a message rather
//!   than silently transposed.
//! - Missing genotypes are imputed to the marker mean, i.e. they contribute `0` to
//!   the centred column, the usual convention for VanRaden `G`.
//! - Which allele the `.bim` calls A1 does not matter: flipping a marker's coding
//!   negates its column of `Z`, and `Z Zᵀ = Σⱼ zⱼ zⱼᵀ` is invariant under that.
//! - [`read_panel`] returns a dense `f64` `Z`: a genotype stored in 2 bits on disk
//!   occupies 64 in memory (10 MiB → 320 MiB at n=40000, m=1000). That is the cost of
//!   that representation, not an artefact of the reader — and it is avoidable.
//! - [`read_packed`] returns the panel in the solver's 2-bit store instead. The `.bed`
//!   layout *is* that store under a different code assignment (4 genotypes per byte,
//!   column-major, `ceil(n/4)` bytes per marker), so the file becomes the working
//!   representation through a byte remap: memory is the size of the file, and no dense
//!   matrix is ever formed. Both routes return the same optimum, missing calls
//!   included.

use crate::error::OcsError;
use crate::packed::{MISSING, PackedGeno};
use faer::Mat;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

/// A marker panel read from a PLINK trio, in the form the solver takes.
pub struct Panel {
    /// Number of individuals (rows of `z`, lines of the `.fam`).
    pub n: usize,
    /// Number of markers (columns of `z`, lines of the `.bim`).
    pub m: usize,
    /// Column-centred genotypes `Z = M - 2p`, shape `n x m`.
    pub z: Mat<f64>,
    /// Allele frequency per marker, estimated from the non-missing genotypes.
    pub p: Vec<f64>,
    /// VanRaden scaling `s = 2 * sum_j p_j (1 - p_j)`.
    pub s: f64,
    /// Individual IDs (`.fam` column 2), in row order of `z`.
    pub ids: Vec<String>,
}

/// `prefix` + `ext`, appended rather than substituted: a prefix such as
/// `data/panel.qc` must yield `data/panel.qc.bed`, which `Path::with_extension`
/// would mangle into `data/panel.bed`.
fn with_suffix(prefix: &Path, ext: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(ext);
    PathBuf::from(s)
}

fn io_err(path: &Path, e: std::io::Error) -> OcsError {
    OcsError::Io(format!("{}: {e}", path.display()))
}

/// Individual IDs from the `.fam` (FID IID PID MID SEX PHENO); its line count is `n`.
fn read_fam(path: &Path) -> Result<Vec<String>, OcsError> {
    let file = File::open(path).map_err(|e| io_err(path, e))?;
    let mut ids = Vec::new();
    for (idx, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|e| io_err(path, e))?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split_whitespace();
        match (fields.next(), fields.next()) {
            (Some(_fid), Some(iid)) => ids.push(iid.to_string()),
            _ => {
                return Err(OcsError::Format(format!(
                    "{}: line {} has fewer than the 2 leading fields (FID IID)",
                    path.display(),
                    idx + 1
                )));
            }
        }
    }
    if ids.is_empty() {
        return Err(OcsError::Format(format!(
            "{}: no individuals",
            path.display()
        )));
    }
    Ok(ids)
}

/// Marker count from the `.bim` (CHR SNP CM BP A1 A2).
fn read_bim(path: &Path) -> Result<usize, OcsError> {
    let file = File::open(path).map_err(|e| io_err(path, e))?;
    let mut m = 0usize;
    for (idx, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|e| io_err(path, e))?;
        if line.trim().is_empty() {
            continue;
        }
        if line.split_whitespace().nth(1).is_none() {
            return Err(OcsError::Format(format!(
                "{}: line {} has fewer than 2 fields",
                path.display(),
                idx + 1
            )));
        }
        m += 1;
    }
    if m == 0 {
        return Err(OcsError::Format(format!("{}: no markers", path.display())));
    }
    Ok(m)
}

/// Open `<prefix>.bed` for decoding: check its length against the `n` and `m` the
/// `.fam` and `.bim` imply, check the magic, and leave the reader on the first
/// genotype byte. Returns the reader and the bytes per marker.
fn open_bed(bed_path: &Path, n: usize, m: usize) -> Result<(BufReader<File>, usize), OcsError> {
    // Rows are padded to whole bytes: 4 genotypes per byte, one row per marker.
    let row_bytes = n.div_ceil(4);
    let expected = 3 + (m as u64) * (row_bytes as u64);
    let actual = std::fs::metadata(bed_path)
        .map_err(|e| io_err(bed_path, e))?
        .len();
    if actual != expected {
        return Err(OcsError::Format(format!(
            "{}: is {actual} bytes, but n={n} (.fam) and m={m} (.bim) imply {expected}",
            bed_path.display()
        )));
    }

    let file = File::open(bed_path).map_err(|e| io_err(bed_path, e))?;
    let mut reader = BufReader::new(file);
    let mut magic = [0u8; 3];
    reader
        .read_exact(&mut magic)
        .map_err(|e| io_err(bed_path, e))?;
    if magic[0] != 0x6c || magic[1] != 0x1b {
        return Err(OcsError::Format(format!(
            "{}: not a PLINK .bed file (magic {:#04x} {:#04x})",
            bed_path.display(),
            magic[0],
            magic[1]
        )));
    }
    if magic[2] != 0x01 {
        return Err(OcsError::Format(format!(
            "{}: individual-major .bed is not supported; rewrite it with PLINK >= 1.9",
            bed_path.display()
        )));
    }
    Ok((reader, row_bytes))
}

/// PLINK's four 2-bit codes remapped to the packed store's, four at a time:
/// `00` (hom A1) → dosage 2, `10` (het) → 1, `11` (hom A2) → 0, `01` → [`MISSING`].
/// Nothing else differs between the two layouts, which is why [`read_packed`] is a
/// byte remap rather than a decode.
const fn bed_to_packed() -> [u8; 256] {
    let mut lut = [0u8; 256];
    let mut byte = 0usize;
    while byte < 256 {
        let mut out = 0u8;
        let mut slot = 0;
        while slot < 4 {
            let code = ((byte >> (2 * slot)) & 0b11) as u8;
            let packed = match code {
                0b00 => 2,
                0b10 => 1,
                0b11 => 0,
                _ => MISSING,
            };
            out |= packed << (2 * slot);
            slot += 1;
        }
        lut[byte] = out;
        byte += 1;
    }
    lut
}

/// A marker panel read straight into the solver's 2-bit store.
pub struct PackedPanel {
    /// Genotypes and the centring derived from them; hand it to
    /// [`crate::support_first::solve_sexed_packed`].
    pub geno: PackedGeno,
    /// Individual IDs, in the `.fam` order, so a support can be reported by name.
    pub ids: Vec<String>,
}

/// Read a PLINK trio into the packed store, forming no dense matrix.
///
/// The genotype block is read once and remapped in place, so peak memory is the size
/// of the `.bed` file — 32× below what [`read_panel`] must hold for the same panel.
/// Allele frequencies come from the non-missing calls, as they do there, and a missing
/// call is imputed to its marker mean exactly (it decodes to `2p`, contributing `0` to
/// the centred column), so both routes reach the same optimum.
pub fn read_packed(prefix: &Path) -> Result<PackedPanel, OcsError> {
    let ids = read_fam(&with_suffix(prefix, ".fam"))?;
    let n = ids.len();
    let m = read_bim(&with_suffix(prefix, ".bim"))?;
    let bed_path = with_suffix(prefix, ".bed");
    let (mut reader, row_bytes) = open_bed(&bed_path, n, m)?;

    let mut data = vec![0u8; row_bytes * m];
    reader
        .read_exact(&mut data)
        .map_err(|e| io_err(&bed_path, e))?;

    let lut = bed_to_packed();
    for byte in data.iter_mut() {
        *byte = lut[*byte as usize];
    }
    // The padding slots past individual n-1 decode to a dosage of 2 under that remap;
    // zero them, so the store holds exactly what packing the same panel would produce.
    let tail = n % 4;
    if tail != 0 {
        let mask = (1u8 << (2 * tail)) - 1;
        for j in 0..m {
            data[j * row_bytes + row_bytes - 1] &= mask;
        }
    }

    let geno = PackedGeno::from_raw_2bit(n, m, data).ok_or_else(|| {
        OcsError::Format(format!(
            "{}: packed genotypes have the wrong length",
            bed_path.display()
        ))
    })?;
    Ok(PackedPanel { geno, ids })
}

/// Read `<prefix>.bed`, `<prefix>.bim` and `<prefix>.fam` into a centred panel.
///
/// The `.bed` length is checked against `n` and `m` before decoding, so a trio whose
/// three files do not belong together fails with a shape message instead of reading
/// shifted genotypes. For large panels prefer [`read_packed`], which holds the same
/// genotypes in 2 bits each instead of 64.
pub fn read_panel(prefix: &Path) -> Result<Panel, OcsError> {
    let fam_path = with_suffix(prefix, ".fam");
    let bim_path = with_suffix(prefix, ".bim");
    let bed_path = with_suffix(prefix, ".bed");

    let ids = read_fam(&fam_path)?;
    let n = ids.len();
    let m = read_bim(&bim_path)?;
    let (mut reader, row_bytes) = open_bed(&bed_path, n, m)?;

    let mut z = Mat::<f64>::zeros(n, m);
    let mut p = vec![0.0_f64; m];
    let mut s = 0.0_f64;
    let mut row = vec![0u8; row_bytes];
    // Dosage of the A1 allele; NAN marks a missing call until the column mean is known.
    let mut dosage = vec![0.0_f64; n];

    for j in 0..m {
        reader
            .read_exact(&mut row)
            .map_err(|e| io_err(&bed_path, e))?;
        let mut sum = 0.0_f64;
        let mut observed = 0usize;
        for (i, d) in dosage.iter_mut().enumerate() {
            // Two bits per genotype, least-significant pair first. Any bits past
            // individual n-1 in the final byte are padding and never read.
            let code = (row[i / 4] >> (2 * (i % 4))) & 0b11;
            *d = match code {
                0b00 => 2.0, // homozygous A1
                0b10 => 1.0, // heterozygous
                0b11 => 0.0, // homozygous A2
                _ => f64::NAN,
            };
            if !d.is_nan() {
                sum += *d;
                observed += 1;
            }
        }
        // An all-missing marker carries no information: p = 0 leaves its column at 0
        // and contributes nothing to s, which is what dropping it would do.
        let pj = if observed == 0 {
            0.0
        } else {
            sum / (2.0 * observed as f64)
        };
        p[j] = pj;
        s += pj * (1.0 - pj);
        let centre = 2.0 * pj;
        for (i, &d) in dosage.iter().enumerate() {
            z[(i, j)] = if d.is_nan() { 0.0 } else { d - centre };
        }
    }
    s *= 2.0;

    Ok(Panel { n, m, z, p, s, ids })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// PLINK code for an A1 dosage; `None` is a missing call. Written from the format
    /// definition rather than by inverting the reader, so a decoding error cannot be
    /// cancelled out by a matching encoding error.
    fn code_of(dosage: Option<u8>) -> u8 {
        match dosage {
            Some(2) => 0b00,
            Some(1) => 0b10,
            Some(0) => 0b11,
            None => 0b01,
            Some(other) => panic!("dosage {other} is not 0, 1 or 2"),
        }
    }

    /// `calls[j][i]` = genotype of individual `i` at marker `j`.
    fn encode_bed(calls: &[Vec<Option<u8>>], n: usize) -> Vec<u8> {
        let mut out = vec![0x6c, 0x1b, 0x01];
        for column in calls {
            let mut row = vec![0u8; n.div_ceil(4)];
            for (i, &call) in column.iter().enumerate() {
                row[i / 4] |= code_of(call) << (2 * (i % 4));
            }
            out.extend_from_slice(&row);
        }
        out
    }

    struct Trio {
        dir: PathBuf,
        prefix: PathBuf,
    }

    impl Drop for Trio {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.dir);
        }
    }

    fn write_trio(tag: &str, calls: &[Vec<Option<u8>>], n: usize) -> Trio {
        let mut dir = std::env::temp_dir();
        dir.push(format!("ocs_rs_plink_{tag}_{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("panel");

        let fam: String = (0..n)
            .map(|i| format!("FAM{i} IND{i} 0 0 1 -9\n"))
            .collect();
        std::fs::write(with_suffix(&prefix, ".fam"), fam).expect("write fam");
        let bim: String = (0..calls.len())
            .map(|j| format!("1 rs{j} 0 {} A G\n", (j + 1) * 1000))
            .collect();
        std::fs::write(with_suffix(&prefix, ".bim"), bim).expect("write bim");
        std::fs::write(with_suffix(&prefix, ".bed"), encode_bed(calls, n)).expect("write bed");

        Trio { dir, prefix }
    }

    /// `G·c` from a dense centred `Z`, written out from the definition so the packed
    /// kernel is checked against arithmetic that shares no code with it.
    fn dense_g_matvec(z: faer::MatRef<'_, f64>, s: f64, ridge: f64, c: &[f64]) -> Vec<f64> {
        let (n, m) = (z.nrows(), z.ncols());
        let t: Vec<f64> = (0..m)
            .map(|j| (0..n).map(|i| z[(i, j)] * c[i]).sum::<f64>())
            .collect();
        (0..n)
            .map(|i| ridge * c[i] + (0..m).map(|j| z[(i, j)] * t[j]).sum::<f64>() / s)
            .collect()
    }

    /// The two routes out of one `.bed` must agree exactly: same kinship, same
    /// optimum. `n = 6` is not a multiple of 4, so every marker's last byte carries
    /// padding the packed route has to blank rather than read as a dosage, and two
    /// markers carry missing calls, which it must impute to the marker mean the way
    /// the dense route does.
    #[test]
    fn packed_route_matches_dense_route() {
        let n = 6;
        let calls: Vec<Vec<Option<u8>>> = vec![
            vec![Some(0), Some(1), Some(2), Some(1), Some(0), Some(2)],
            vec![Some(2), None, Some(1), Some(0), Some(2), Some(1)],
            vec![Some(1), Some(1), Some(0), Some(2), Some(1), Some(0)],
            vec![Some(0), Some(2), None, Some(1), None, Some(2)],
            vec![Some(2), Some(0), Some(1), Some(1), Some(2), Some(0)],
        ];
        let trio = write_trio("packed_vs_dense", &calls, n);
        let dense = read_panel(&trio.prefix).expect("dense read");
        let packed = read_packed(&trio.prefix).expect("packed read");

        assert_eq!(packed.geno.n(), dense.n);
        assert_eq!(packed.geno.m(), dense.m);
        assert_eq!(packed.ids, dense.ids);
        // The store is the genotype block of the file itself, byte for byte.
        assert_eq!(packed.geno.packed_bytes(), n.div_ceil(4) * calls.len());

        let ridge = 1e-5;
        let c: Vec<f64> = (0..n).map(|i| 0.05 + 0.1 * i as f64).collect();
        let want = dense_g_matvec(dense.z.as_ref(), dense.s, ridge, &c);
        let got = packed.geno.g_matvec(&c, ridge);
        for (i, (&w, &g)) in want.iter().zip(&got).enumerate() {
            assert!((w - g).abs() < 1e-12, "G·c[{i}]: dense {w}, packed {g}");
        }
        for i in 0..n {
            for j in 0..n {
                let w = {
                    let raw: f64 = (0..dense.m)
                        .map(|l| dense.z[(i, l)] * dense.z[(j, l)])
                        .sum();
                    raw / dense.s + if i == j { ridge } else { 0.0 }
                };
                let g = packed.geno.gram_entry(i, j, ridge);
                assert!((w - g).abs() < 1e-12, "G[{i},{j}]: dense {w}, packed {g}");
            }
        }

        let male: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();
        let b = [0.9_f64, 0.4, 0.7, 0.1, 0.5, 0.8];
        let mean_diag: f64 = (0..n)
            .map(|i| (0..dense.m).map(|l| dense.z[(i, l)].powi(2)).sum::<f64>() / dense.s)
            .sum::<f64>()
            / n as f64;
        let k = 0.6 * mean_diag;
        let from_dense = crate::support_first::solve_sexed(
            dense.z.as_ref(),
            dense.s,
            ridge,
            &b,
            &male,
            k,
            10_000,
            1e-10,
        );
        let from_packed = crate::support_first::solve_sexed_packed(
            &packed.geno,
            ridge,
            &b,
            &male,
            k,
            10_000,
            1e-10,
        );
        assert_eq!(from_dense.status, from_packed.status);
        assert_eq!(from_dense.support, from_packed.support);
        assert!(
            (from_dense.gain - from_packed.gain).abs() < 1e-12,
            "gain: dense {}, packed {}",
            from_dense.gain,
            from_packed.gain
        );
    }

    /// A panel whose markers are entirely missing has no frequencies to estimate; the
    /// packed route must survive it the way the dense one does (`p = 0`, zero column).
    #[test]
    fn packed_route_handles_an_all_missing_marker() {
        let n = 4;
        let calls: Vec<Vec<Option<u8>>> = vec![
            vec![Some(0), Some(1), Some(2), Some(1)],
            vec![None, None, None, None],
            vec![Some(2), Some(0), Some(1), Some(1)],
        ];
        let trio = write_trio("packed_all_missing", &calls, n);
        let dense = read_panel(&trio.prefix).expect("dense read");
        let packed = read_packed(&trio.prefix).expect("packed read");
        let ridge = 1e-5;
        let c = [0.3_f64, 0.1, 0.4, 0.2];
        let want = dense_g_matvec(dense.z.as_ref(), dense.s, ridge, &c);
        let got = packed.geno.g_matvec(&c, ridge);
        for (&w, &g) in want.iter().zip(&got) {
            assert!((w - g).abs() < 1e-12, "dense {w}, packed {g}");
        }
    }

    #[test]
    fn round_trips_genotypes_frequencies_and_scale() {
        // n = 6 is deliberately not a multiple of 4, so the last byte of every row
        // carries padding bits the reader must ignore; marker 1 has a missing call.
        let n = 6;
        let calls: Vec<Vec<Option<u8>>> = vec![
            vec![Some(0), Some(1), Some(2), Some(2), Some(1), Some(0)],
            vec![Some(2), None, Some(0), Some(1), Some(1), Some(2)],
            vec![Some(0), Some(0), Some(0), Some(0), Some(0), Some(0)],
        ];
        let trio = write_trio("roundtrip", &calls, n);

        let panel = read_panel(&trio.prefix).expect("read panel");

        assert_eq!(panel.n, n);
        assert_eq!(panel.m, calls.len());
        assert_eq!(panel.ids[0], "IND0");
        assert_eq!(panel.ids[5], "IND5");

        // Frequencies computed from the calls, independently of the reader.
        let mut expect_s = 0.0;
        for (j, column) in calls.iter().enumerate() {
            let observed: Vec<u8> = column.iter().filter_map(|c| *c).collect();
            let pj =
                observed.iter().map(|&d| d as f64).sum::<f64>() / (2.0 * observed.len() as f64);
            assert!(
                (panel.p[j] - pj).abs() < 1e-12,
                "marker {j}: {} vs {pj}",
                panel.p[j]
            );
            expect_s += pj * (1.0 - pj);
            for (i, &call) in column.iter().enumerate() {
                let want = match call {
                    Some(d) => d as f64 - 2.0 * pj,
                    None => 0.0,
                };
                assert!(
                    (panel.z[(i, j)] - want).abs() < 1e-12,
                    "z[{i},{j}] = {} but expected {want}",
                    panel.z[(i, j)]
                );
            }
        }
        assert!((panel.s - 2.0 * expect_s).abs() < 1e-12);
    }

    #[test]
    fn rejects_a_file_that_is_not_a_bed() {
        let calls = vec![vec![Some(0), Some(1), Some(2), Some(1)]];
        let trio = write_trio("magic", &calls, 4);
        let bed = with_suffix(&trio.prefix, ".bed");
        let mut bytes = std::fs::read(&bed).expect("read bed");
        bytes[0] = 0x00;
        std::fs::write(&bed, bytes).expect("rewrite bed");

        match read_panel(&trio.prefix) {
            Err(OcsError::Format(msg)) => assert!(msg.contains("not a PLINK"), "{msg}"),
            other => panic!(
                "expected a format error, got {other:?}",
                other = other.map(|p| p.n)
            ),
        }
    }

    #[test]
    fn rejects_a_mismatched_trio() {
        // A .fam with one individual too many: the .bed can no longer be the right
        // size, and the reader must say so instead of decoding shifted genotypes.
        let calls = vec![vec![Some(0), Some(1), Some(2), Some(1)]];
        let trio = write_trio("mismatch", &calls, 4);
        let fam = with_suffix(&trio.prefix, ".fam");
        let mut text = std::fs::read_to_string(&fam).expect("read fam");
        text.push_str("FAM9 IND9 0 0 2 -9\n");
        std::fs::write(&fam, text).expect("rewrite fam");

        match read_panel(&trio.prefix) {
            Err(OcsError::Format(msg)) => assert!(msg.contains("imply"), "{msg}"),
            other => panic!(
                "expected a format error, got {other:?}",
                other = other.map(|p| p.n)
            ),
        }
    }

    #[test]
    fn rejects_individual_major_layout() {
        let calls = vec![vec![Some(0), Some(1), Some(2), Some(1)]];
        let trio = write_trio("major", &calls, 4);
        let bed = with_suffix(&trio.prefix, ".bed");
        let mut bytes = std::fs::read(&bed).expect("read bed");
        bytes[2] = 0x00;
        std::fs::write(&bed, bytes).expect("rewrite bed");

        match read_panel(&trio.prefix) {
            Err(OcsError::Format(msg)) => assert!(msg.contains("individual-major"), "{msg}"),
            other => panic!(
                "expected a format error, got {other:?}",
                other = other.map(|p| p.n)
            ),
        }
    }
}
