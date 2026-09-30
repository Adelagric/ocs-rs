//! The real PIC pig panel through the packed `.bed` route, and what 2 bits cost there.
//!
//! The panel ships as imputed *dosages* (mostly 0/1/2, a small fraction fractional).
//! PLINK `.bed` — and any 2-bit store — holds hard calls only, so the packed route must
//! round. This measures both things at once: the memory the packed route saves, and the
//! price the rounding puts on the optimum, against the raw-dosage solve the manuscript's
//! pig results use.
//!
//!   cargo run --release -- <FileS1 dir> [k_frac]
use std::io::{BufRead, BufWriter, Write};
use std::time::Instant;

use faer::Mat;
use ocs_rs::support_first::{solve_sexed, solve_sexed_packed};
use ocs_rs::{plink, packed::PackedGeno};

const RIDGE: f64 = 1e-5;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::path::PathBuf::from(std::env::args().nth(1).ok_or("usage: pig_bed <FileS1 dir>")?);
    let k_frac: f64 = std::env::args().nth(2).map(|s| s.parse()).transpose()?.unwrap_or(0.05);
    let geno = dir.join("genotypes.txt");

    // --- one pass: raw dosages into a dense matrix, hard calls into a .bed buffer ---
    let t = Instant::now();
    let f = std::io::BufReader::with_capacity(1 << 22, std::fs::File::open(&geno)?);
    let mut lines = f.lines();
    let header = lines.next().ok_or("empty genotypes.txt")??;
    let m = header.split(',').count() - 1;
    let mut ids: Vec<String> = Vec::new();
    let mut dosages: Vec<f64> = Vec::new(); // row-major while reading
    let mut fractional = 0usize;
    let mut total = 0usize;
    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let mut it = line.split(',');
        ids.push(it.next().ok_or("no id")?.to_string());
        let before = dosages.len();
        for tok in it {
            let v: f64 = tok.parse()?;
            if (v - v.round()).abs() > 1e-9 {
                fractional += 1;
            }
            dosages.push(v);
            total += 1;
        }
        if dosages.len() - before != m {
            return Err(format!("row {} has {} markers, header says {m}", ids.len(), dosages.len() - before).into());
        }
    }
    let n = ids.len();
    let t_read = t.elapsed().as_secs_f64();
    println!("PIC pig panel: n={n} m={m}  read in {t_read:.1}s");
    println!(
        "  imputed dosages: {fractional} of {total} values are non-integer ({:.3}%) — a 2-bit store cannot hold them",
        100.0 * fractional as f64 / total as f64
    );

    // --- write a PLINK trio with the rounded hard calls (what `plink --make-bed` would) ---
    let prefix = std::env::temp_dir().join("pig_panel");
    {
        let mut fam = BufWriter::new(std::fs::File::create(prefix.with_extension("fam"))?);
        for (i, id) in ids.iter().enumerate() {
            writeln!(fam, "F{id} {id} 0 0 {} -9", if i % 2 == 0 { 1 } else { 2 })?;
        }
        let mut bim = BufWriter::new(std::fs::File::create(prefix.with_extension("bim"))?);
        for j in 0..m {
            writeln!(bim, "1 snp{j} 0 {} A G", j + 1)?;
        }
        let mut bed = BufWriter::new(std::fs::File::create(prefix.with_extension("bed"))?);
        bed.write_all(&[0x6c, 0x1b, 0x01])?;
        let row_bytes = n.div_ceil(4);
        let mut block = vec![0u8; row_bytes];
        for j in 0..m {
            block.iter_mut().for_each(|b| *b = 0);
            for i in 0..n {
                let code = match dosages[i * m + j].round() as i64 {
                    2 => 0b00u8,
                    1 => 0b10,
                    _ => 0b11,
                };
                block[i >> 2] |= code << ((i & 3) * 2);
            }
            bed.write_all(&block)?;
        }
    }
    let bed_bytes = std::fs::metadata(prefix.with_extension("bed"))?.len();
    println!("  wrote {} ({:.1} MB)", prefix.with_extension("bed").display(), bed_bytes as f64 / 1e6);

    // --- the raw-dosage solve: centre in place, as the R pipeline does ---
    let mut p = vec![0.0f64; m];
    for j in 0..m {
        p[j] = (0..n).map(|i| dosages[i * m + j]).sum::<f64>() / (2.0 * n as f64);
    }
    let s: f64 = 2.0 * p.iter().map(|pj| pj * (1.0 - pj)).sum::<f64>();
    let z_raw = Mat::from_fn(n, m, |i, j| dosages[i * m + j] - 2.0 * p[j]);
    drop(dosages);
    let male: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();
    // A deterministic stand-in criterion: the manuscript uses the panel's real EBV, but
    // the comparison here is representation against representation on one criterion.
    let b: Vec<f64> = (0..n).map(|i| ((i * 2654435761) % 1000) as f64 / 1000.0).collect();
    let mean_diag: f64 = (0..n)
        .map(|i| (0..m).map(|l| z_raw[(i, l)].powi(2)).sum::<f64>() / s)
        .sum::<f64>()
        / n as f64;
    let k = k_frac * mean_diag;
    println!("  cap k = {k_frac} x mean diag G = {k:.6}");

    let t = Instant::now();
    let raw = solve_sexed(z_raw.as_ref(), s, RIDGE, &b, &male, k, 50_000, 1e-9);
    println!(
        "  raw dosages, dense f64 : {:>8.1} MB  solve {:.2}s  |S|={}  gain {:.9}",
        (n * m * 8) as f64 / 1e6, t.elapsed().as_secs_f64(), raw.support.len(), raw.gain
    );
    drop(z_raw);

    // --- the same panel through .bed, both routes ---
    let t = Instant::now();
    let dense = plink::read_panel(&prefix)?;
    let t_dense = t.elapsed().as_secs_f64();
    let hard = solve_sexed(dense.z.as_ref(), dense.s, RIDGE, &b, &male, k, 50_000, 1e-9);
    println!(
        "  rounded, dense f64     : {:>8.1} MB  read {t_dense:.2}s  |S|={}  gain {:.9}",
        (n * m * 8) as f64 / 1e6, hard.support.len(), hard.gain
    );
    let t = Instant::now();
    let pk = plink::read_packed(&prefix)?;
    let t_pk = t.elapsed().as_secs_f64();
    let pkout = solve_sexed_packed(&pk.geno, RIDGE, &b, &male, k, 50_000, 1e-9);
    println!(
        "  rounded, 2-bit packed  : {:>8.1} MB  read {t_pk:.2}s  |S|={}  gain {:.9}",
        pk.geno.packed_bytes() as f64 / 1e6, pkout.support.len(), pkout.gain
    );
    println!(
        "  packed vs dense-from-bed: support match {}  Δgain {:.2e}",
        hard.support == pkout.support,
        (hard.gain - pkout.gain).abs()
    );
    println!(
        "  rounding cost (raw vs rounded): support {} vs {}, shared {}, Δgain {:.2e}",
        raw.support.len(), hard.support.len(),
        raw.support.iter().filter(|i| hard.support.contains(i)).count(),
        (raw.gain - hard.gain).abs()
    );
    println!("  (a dense n x n G for this panel: {:.1} MB)", (n * n * 8) as f64 / 1e6);
    for ext in ["bed", "bim", "fam"] {
        std::fs::remove_file(prefix.with_extension(ext)).ok();
    }
    Ok(())
}
