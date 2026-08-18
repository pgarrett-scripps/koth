use std::collections::BTreeMap;
use std::path::Path;

use serde::Serialize;

use crate::error::KothError;
use crate::models::{Hill, ScoredFeature};

// ── helpers ──────────────────────────────────────────────────────────────────

fn median(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = v.len();
    if n % 2 == 0 {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    } else {
        v[n / 2]
    }
}

fn stats(values: &mut [f64]) -> StatSummary {
    StatSummary {
        median: median(values),
        std: crate::stats::std_dev(values),
    }
}

// ── serialisable types ────────────────────────────────────────────────────────

#[derive(Serialize)]
pub struct StatSummary {
    pub median: f64,
    pub std: f64,
}

#[derive(Serialize)]
pub struct HillsReport {
    pub n: usize,
    pub mz_std: StatSummary,
    pub im_std: StatSummary,
    pub rt_width: StatSummary,
    pub n_scans: StatSummary,
    pub intensity_sum: StatSummary,
}

#[derive(Serialize)]
pub struct FeaturesReport {
    pub n: usize,
    /// Charge state → feature count, sorted by charge.
    pub charge_distribution: BTreeMap<u8, usize>,
    /// Chromatographic cosine score (mean adjacent-isotope-hill cosine).
    pub cosine_score: StatSummary,
    /// Bhattacharyya isotope-pattern score.
    pub isotope_score: StatSummary,
    /// `isotope_score × cosine_score`.
    pub combined_score: StatSummary,
    pub ppm_error: StatSummary,
    pub n_isotopes: StatSummary,
    pub rt_width: StatSummary,
}

#[derive(Serialize)]
pub struct RunReport {
    pub hills: HillsReport,
    pub features: FeaturesReport,
}

// ── builders ──────────────────────────────────────────────────────────────────

pub fn build_hills_report(hills: &[Hill]) -> HillsReport {
    let mut mz_stds = hills.iter().map(|h| h.mz_std).collect::<Vec<_>>();
    let mut im_stds = hills.iter().map(|h| h.im_std).collect::<Vec<_>>();
    let mut rt_widths = hills.iter().map(|h| h.rt_width).collect::<Vec<_>>();
    let mut n_scans = hills.iter().map(|h| h.n_scans as f64).collect::<Vec<_>>();
    let mut int_sums = hills.iter().map(|h| h.intensity_sum).collect::<Vec<_>>();

    HillsReport {
        n: hills.len(),
        mz_std: stats(&mut mz_stds),
        im_std: stats(&mut im_stds),
        rt_width: stats(&mut rt_widths),
        n_scans: stats(&mut n_scans),
        intensity_sum: stats(&mut int_sums),
    }
}

pub fn build_features_report(features: &[ScoredFeature]) -> FeaturesReport {
    let scored: Vec<&ScoredFeature> = features.iter().filter(|f| f.feature.charge > 0).collect();

    let mut charge_dist: BTreeMap<u8, usize> = BTreeMap::new();
    for sf in &scored {
        *charge_dist.entry(sf.feature.charge).or_insert(0) += 1;
    }

    let mut cosine = scored.iter().map(|sf| sf.cosine_score).collect::<Vec<_>>();
    let mut isotope = scored.iter().map(|sf| sf.isotope_score).collect::<Vec<_>>();
    let mut combined = scored
        .iter()
        .map(|sf| sf.combined_score)
        .collect::<Vec<_>>();
    let mut ppm = scored
        .iter()
        .map(|sf| sf.feature.ppm_error)
        .collect::<Vec<_>>();
    let mut n_iso = scored
        .iter()
        .map(|sf| sf.feature.hills.len() as f64)
        .collect::<Vec<_>>();
    let mut rt_w = scored
        .iter()
        .map(|sf| sf.feature.rt_end() - sf.feature.rt_start())
        .collect::<Vec<_>>();

    FeaturesReport {
        n: scored.len(),
        charge_distribution: charge_dist,
        cosine_score: stats(&mut cosine),
        isotope_score: stats(&mut isotope),
        combined_score: stats(&mut combined),
        ppm_error: stats(&mut ppm),
        n_isotopes: stats(&mut n_iso),
        rt_width: stats(&mut rt_w),
    }
}

// ── writer ────────────────────────────────────────────────────────────────────

pub fn write_report(report: &RunReport, path: &Path) -> Result<(), KothError> {
    let json = serde_json::to_string_pretty(report)?;
    std::fs::write(path, json)?;
    log::info!("Wrote report to {}", path.display());
    Ok(())
}
