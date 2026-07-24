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
//! - `Z` is dense `f64`: a genotype stored in 2 bits on disk occupies 64 in memory
//!   (10 MiB → 320 MiB at n=40000, m=1000). That is the cost of the representation
//!   the solver consumes, not an artefact of the reader.

use crate::error::OcsError;
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

/// Read `<prefix>.bed`, `<prefix>.bim` and `<prefix>.fam` into a centred panel.
///
/// The `.bed` length is checked against `n` and `m` before decoding, so a trio whose
/// three files do not belong together fails with a shape message instead of reading
/// shifted genotypes.
pub fn read_panel(prefix: &Path) -> Result<Panel, OcsError> {
    let fam_path = with_suffix(prefix, ".fam");
    let bim_path = with_suffix(prefix, ".bim");
    let bed_path = with_suffix(prefix, ".bed");

    let ids = read_fam(&fam_path)?;
    let n = ids.len();
    let m = read_bim(&bim_path)?;

    // Rows are padded to whole bytes: 4 genotypes per byte, one row per marker.
    let row_bytes = n.div_ceil(4);
    let expected = 3 + (m as u64) * (row_bytes as u64);
    let actual = std::fs::metadata(&bed_path)
        .map_err(|e| io_err(&bed_path, e))?
        .len();
    if actual != expected {
        return Err(OcsError::Format(format!(
            "{}: is {actual} bytes, but n={n} (.fam) and m={m} (.bim) imply {expected}",
            bed_path.display()
        )));
    }

    let file = File::open(&bed_path).map_err(|e| io_err(&bed_path, e))?;
    let mut reader = BufReader::new(file);
    let mut magic = [0u8; 3];
    reader
        .read_exact(&mut magic)
        .map_err(|e| io_err(&bed_path, e))?;
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
