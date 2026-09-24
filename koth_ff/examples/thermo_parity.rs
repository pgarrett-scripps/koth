//! Compare the native Thermo `.raw` reader against an mzML conversion of the
//! same run (msconvert with vendor peak picking).
//!
//! ```text
//! cargo run --release --example thermo_parity -- ms1  RUN.raw RUN.mzML[.gz]
//! cargo run --release --example thermo_parity -- ms2  RUN.raw RUN.mzML[.gz]
//! cargo run --release --example thermo_parity -- load RUN.raw|RUN.mzML
//! ```
//!
//! `ms1` pairs MS1 spectra in RT order and reports scan counts, RT, per-scan
//! peak counts and m/z / intensity agreement. `ms2` compares the DIA isolation
//! windows each reader recovers. `load` reads MS1 once and reports the time.

use std::collections::BTreeMap;
use std::path::Path;
use std::time::Instant;

use koth_ff::config::FileConfig;
use koth_ff::io;
use koth_ff::models::Spectrum;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("ms1") => ms1(Path::new(&args[2]), Path::new(&args[3])),
        Some("ms2") => ms2(Path::new(&args[2]), Path::new(&args[3])),
        Some("load") => load(Path::new(&args[2])),
        _ => anyhow::bail!("usage: thermo_parity ms1|ms2 RAW MZML | load PATH"),
    }
}

fn load(path: &Path) -> anyhow::Result<()> {
    let t = Instant::now();
    let mut n = 0usize;
    let mut peaks = 0usize;
    for s in io::stream_spectra(path, &FileConfig::default())? {
        let s = s?;
        n += 1;
        peaks += s.peaks.len();
    }
    println!(
        "load {}: {n} MS1 spectra, {peaks} peaks, {:.2} s",
        path.display(),
        t.elapsed().as_secs_f64()
    );
    Ok(())
}

fn read_timed(path: &Path) -> anyhow::Result<Vec<Spectrum>> {
    let t = Instant::now();
    let spectra = io::read_spectra(path, &FileConfig::default())?;
    println!(
        "read {:>6} MS1 spectra in {:6.2} s from {}",
        spectra.len(),
        t.elapsed().as_secs_f64(),
        path.display()
    );
    Ok(spectra)
}

fn ms1(raw: &Path, mzml: &Path) -> anyhow::Result<()> {
    let a = read_timed(raw)?;
    let b = read_timed(mzml)?;
    println!("MS1 spectra: raw {} mzML {}", a.len(), b.len());
    let n = a.len().min(b.len());
    let mut max_rt = 0.0f64;
    let mut same_count = 0usize;
    let (mut peaks_a, mut peaks_b) = (0usize, 0usize);
    let mut max_ppm = 0.0f64;
    let mut max_rel_int = 0.0f64;
    let mut compared = 0usize;
    let mut faims_mismatch = 0usize;
    for (sa, sb) in a.iter().zip(b.iter()).take(n) {
        max_rt = max_rt.max((sa.retention_time - sb.retention_time).abs());
        peaks_a += sa.peaks.len();
        peaks_b += sb.peaks.len();
        if sa.faims_cv != sb.faims_cv {
            faims_mismatch += 1;
        }
        if sa.peaks.len() != sb.peaks.len() {
            continue;
        }
        same_count += 1;
        for (pa, pb) in sa.peaks.iter().zip(&sb.peaks) {
            let ppm = ((pa.mz - pb.mz) as f64).abs() / pb.mz as f64 * 1e6;
            max_ppm = max_ppm.max(ppm);
            let rel =
                ((pa.intensity - pb.intensity) as f64).abs() / (pb.intensity as f64).max(1e-12);
            max_rel_int = max_rel_int.max(rel);
            compared += 1;
        }
    }
    println!("paired {n}; max |dRT| = {:.3e} min", max_rt);
    println!(
        "identical peak count in {same_count}/{n} scans; total peaks raw {peaks_a} mzML {peaks_b}"
    );
    println!(
        "over {compared} paired peaks: max |dm/z| = {max_ppm:.4} ppm, max rel |dI| = {max_rel_int:.3e}"
    );
    println!("FAIMS CV mismatches: {faims_mismatch}");
    Ok(())
}

fn ms2(raw: &Path, mzml: &Path) -> anyhow::Result<()> {
    let t = Instant::now();
    let a = io::thermo::read_thermo_ms2(raw)?;
    println!(
        "raw: {} MS2 spectra in {:.2} s",
        a.len(),
        t.elapsed().as_secs_f64()
    );
    let b: Vec<Spectrum> = io::mzml::stream_mzml_ms2(mzml)?.collect();
    println!("mzML: {} MS2 spectra", b.len());
    let windows = |v: &[Spectrum]| {
        let mut m: BTreeMap<(i64, i64, i64), usize> = BTreeMap::new();
        for s in v {
            *m.entry(s.isolation_window.unwrap().key()).or_default() += 1;
        }
        m
    };
    let (wa, wb) = (windows(&a), windows(&b));
    println!("distinct windows: raw {} mzML {}", wa.len(), wb.len());
    let fmt = |k: &(i64, i64, i64)| {
        format!(
            "{:.3} [{:.3}, {:.3}]",
            k.0 as f64 / 1e4,
            k.1 as f64 / 1e4,
            k.2 as f64 / 1e4
        )
    };
    for (k, c) in wa.iter().take(5) {
        println!("  raw  {} x{c}", fmt(k));
    }
    for (k, c) in wb.iter().take(5) {
        println!("  mzML {} x{c}", fmt(k));
    }
    // Pair scans in RT order and compare windows and peaks.
    let n = a.len().min(b.len());
    let (mut same_w, mut same_n, mut max_rt, mut max_bound) = (0usize, 0usize, 0.0f64, 0.0f64);
    for (sa, sb) in a.iter().zip(&b).take(n) {
        let (wa, wb) = (sa.isolation_window.unwrap(), sb.isolation_window.unwrap());
        max_rt = max_rt.max((sa.retention_time - sb.retention_time).abs());
        let d = (wa.target - wb.target)
            .abs()
            .max((wa.lower - wb.lower).abs())
            .max((wa.upper - wb.upper).abs());
        max_bound = max_bound.max(d);
        if wa.key() == wb.key() {
            same_w += 1;
        }
        if sa.peaks.len() == sb.peaks.len() {
            same_n += 1;
        }
    }
    println!(
        "paired {n}: identical window key {same_w}, identical peak count {same_n}, \
         max |dRT| {max_rt:.3e} min, max window bound diff {max_bound:.4} Th"
    );
    Ok(())
}
