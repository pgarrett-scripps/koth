//! JSON diagnostics report for the `koth_align` stage.
//!
//! Captures per-anchor data + fitted models so downstream plotting scripts can
//! visualise RT warp, mass-ppm drift, IM drift, and per-run quality without
//! re-running anything. Consensus feature data is intentionally NOT duplicated
//! here — it already lives in `consensus_features.tsv` / `intensity_matrix.tsv`.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::time::Instant;

use rand::seq::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;
use serde::Serialize;

/// Max anchor records emitted per run in the JSON report.
/// `residual_stats` still aggregates over every anchor; only the per-anchor
/// detail array gets capped. With 20 runs and ~70k anchors each, the full
/// list would balloon the JSON past 800 MB and turn the writer into the new
/// bottleneck — see lfq/mod.rs for the matching loop-hoisting fix.
const MAX_ANCHORS_PER_RUN: usize = 500;

use crate::alignment::{
    drift::DriftFit, warp::RtWarp, AlignmentResult, RunAlignment, RunInput,
};
use crate::config::AlignConfig;
use crate::error::KothError;
use crate::lfq::IntensityMatrix;

// ── helpers ──────────────────────────────────────────────────────────────────

fn median_sorted(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    let idx = ((p * (n as f64 - 1.0)).round() as usize).min(n - 1);
    sorted[idx]
}

fn summary(values: &[f64]) -> ResidualSummary {
    if values.is_empty() {
        return ResidualSummary {
            n: 0,
            median: f64::NAN,
            mad: f64::NAN,
            std: f64::NAN,
        };
    }
    let mut sorted: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let med = median_sorted(&sorted);
    let mut abs_dev: Vec<f64> = sorted.iter().map(|x| (x - med).abs()).collect();
    abs_dev.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad = median_sorted(&abs_dev);
    ResidualSummary {
        n: sorted.len(),
        median: med,
        mad,
        std: crate::stats::std_dev(&sorted),
    }
}

// ── serialisable types ───────────────────────────────────────────────────────

#[derive(Serialize)]
pub struct ResidualSummary {
    pub n: usize,
    pub median: f64,
    pub mad: f64,
    pub std: f64,
}

#[derive(Serialize)]
pub struct AnchorRecord {
    pub ref_rt_norm: f64,
    pub run_rt_norm: f64,
    pub ref_rt: f64,
    pub run_rt: f64,
    pub ref_mz: f64,
    pub run_mz: f64,
    pub ppm_error: f64,
    pub ref_im: f64,
    pub run_im: f64,
    pub im_delta: f64,
    pub rt_residual_after_warp: f64,
    pub ppm_residual_after_fit: f64,
    pub im_residual_after_fit: f64,
    pub rt_active: bool,
    pub mass_active: bool,
    pub im_active: bool,
}

#[derive(Serialize)]
pub struct ResidualStats {
    pub ppm_before: ResidualSummary,
    pub ppm_after: ResidualSummary,
    pub im_before: ResidualSummary,
    pub im_after: ResidualSummary,
    pub rt_norm_before: ResidualSummary,
    pub rt_norm_after: ResidualSummary,
}

#[derive(Serialize)]
pub struct AlignmentDiagnostics {
    pub n_anchors: usize,
    pub n_anchors_rt_active: usize,
    pub n_anchors_mass_active: usize,
    pub n_anchors_im_active: usize,
    pub rt_warp: RtWarp,
    pub mass_drift: DriftFit,
    pub im_drift: DriftFit,
    pub anchors: Vec<AnchorRecord>,
    pub residual_stats: ResidualStats,
}

#[derive(Serialize)]
pub struct LfqRunSummary {
    pub n_detected: usize,
    pub detection_rate: f64,
    pub qvalue_percentiles: QvaluePercentiles,
}

#[derive(Serialize)]
pub struct QvaluePercentiles {
    pub p25: f64,
    pub p50: f64,
    pub p75: f64,
    pub p95: f64,
}

#[derive(Serialize)]
pub struct RunAlignReport {
    pub name: String,
    pub is_reference: bool,
    pub rt_range: [f64; 2],
    pub n_features: usize,
    pub n_features_high_score: usize,
    pub n_hills: usize,
    /// Absent for the reference run.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub alignment: Option<AlignmentDiagnostics>,
    pub lfq: LfqRunSummary,
}

#[derive(Serialize)]
pub struct Timing {
    pub alignment_sec: f64,
    pub lfq_sec: f64,
    pub total_sec: f64,
}

#[derive(Serialize)]
pub struct AlignReport {
    pub schema_version: u32,
    pub reference_run: String,
    pub n_runs: usize,
    pub n_consensus_features: usize,
    pub config: AlignConfig,
    pub runs: Vec<RunAlignReport>,
    pub timing: Timing,
}

// ── builder ──────────────────────────────────────────────────────────────────

fn build_anchors_block(al: &RunAlignment) -> AlignmentDiagnostics {
    let (run_lo, run_hi) = al.run_rt_range;
    let (ref_lo, ref_hi) = al.ref_rt_range;
    let run_span = run_hi - run_lo;
    let ref_span = ref_hi - ref_lo;

    let mut records: Vec<AnchorRecord> = Vec::with_capacity(al.anchors.len());

    let mut ppm_before: Vec<f64> = Vec::with_capacity(al.anchors.len());
    let mut ppm_after: Vec<f64> = Vec::with_capacity(al.anchors.len());
    let mut im_before: Vec<f64> = Vec::new();
    let mut im_after: Vec<f64> = Vec::new();
    let mut rt_before: Vec<f64> = Vec::with_capacity(al.anchors.len());
    let mut rt_after: Vec<f64> = Vec::with_capacity(al.anchors.len());

    for (i, a) in al.anchors.iter().enumerate() {
        let ppm = a.ppm_error();
        let im_d = a.im_delta();
        let rt_norm_delta = a.ref_rt_norm - a.run_rt_norm;

        // Apply the fitted models to compute residuals.
        let warp_predicted_delta = al.rt_warp.delta_at_run_norm(a.run_rt_norm);
        let rt_resid = rt_norm_delta - warp_predicted_delta;

        let ppm_resid = ppm - al.mass_drift.predict(a.ref_rt_norm);

        // IM is only meaningful when both sides have IM.
        let im_resid = if a.ref_im != 0.0 && a.run_im != 0.0 {
            im_d - al.im_drift.predict(a.ref_rt_norm)
        } else {
            f64::NAN
        };

        // Absolute RTs derived from the per-run RT span.
        let run_rt_abs = run_lo + a.run_rt_norm * run_span;
        let ref_rt_abs = ref_lo + a.ref_rt_norm * ref_span;

        let rt_active = *al.rt_active.get(i).unwrap_or(&false);
        let mass_active = *al.mass_active.get(i).unwrap_or(&false);
        let im_active = *al.im_active.get(i).unwrap_or(&false);

        records.push(AnchorRecord {
            ref_rt_norm: a.ref_rt_norm,
            run_rt_norm: a.run_rt_norm,
            ref_rt: ref_rt_abs,
            run_rt: run_rt_abs,
            ref_mz: a.ref_mz,
            run_mz: a.run_mz,
            ppm_error: ppm,
            ref_im: a.ref_im,
            run_im: a.run_im,
            im_delta: im_d,
            rt_residual_after_warp: rt_resid,
            ppm_residual_after_fit: ppm_resid,
            im_residual_after_fit: im_resid,
            rt_active,
            mass_active,
            im_active,
        });

        ppm_before.push(ppm);
        ppm_after.push(ppm_resid);
        if a.ref_im != 0.0 && a.run_im != 0.0 {
            im_before.push(im_d);
            im_after.push(im_resid);
        }
        rt_before.push(rt_norm_delta);
        rt_after.push(rt_resid);
    }

    // Subsample the per-anchor detail to keep the JSON report tractable.
    // Residual stats above were computed over the full anchor set, so the
    // summary numbers stay exact — only the scatter-plot detail is capped.
    // Deterministic seed so re-running on the same input produces the same
    // sample (helpful for diffing reports).
    let n_total = records.len();
    let sampled_records = if n_total > MAX_ANCHORS_PER_RUN {
        let mut rng = rand::rngs::StdRng::seed_from_u64(0xa11c_0_de);
        records.shuffle(&mut rng);
        records.truncate(MAX_ANCHORS_PER_RUN);
        records
    } else {
        records
    };

    AlignmentDiagnostics {
        n_anchors: al.n_anchors,
        n_anchors_rt_active: al.rt_active.iter().filter(|&&v| v).count(),
        n_anchors_mass_active: al.mass_active.iter().filter(|&&v| v).count(),
        n_anchors_im_active: al.im_active.iter().filter(|&&v| v).count(),
        rt_warp: al.rt_warp.clone(),
        mass_drift: al.mass_drift.clone(),
        im_drift: al.im_drift.clone(),
        anchors: sampled_records,
        residual_stats: ResidualStats {
            ppm_before: summary(&ppm_before),
            ppm_after: summary(&ppm_after),
            im_before: summary(&im_before),
            im_after: summary(&im_after),
            rt_norm_before: summary(&rt_before),
            rt_norm_after: summary(&rt_after),
        },
    }
}

fn build_lfq_summary(matrix: &IntensityMatrix, run_idx: usize) -> LfqRunSummary {
    let mut detected_q: Vec<f64> = Vec::new();
    for feat in 0..matrix.n_features {
        let intensity = matrix.intensity(feat, run_idx);
        if intensity > 0.0 {
            detected_q.push(matrix.q_value(feat, run_idx));
        }
    }
    let n_detected = detected_q.len();
    detected_q.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let qp = QvaluePercentiles {
        p25: percentile_sorted(&detected_q, 0.25),
        p50: percentile_sorted(&detected_q, 0.50),
        p75: percentile_sorted(&detected_q, 0.75),
        p95: percentile_sorted(&detected_q, 0.95),
    };

    let detection_rate = if matrix.n_features == 0 {
        0.0
    } else {
        n_detected as f64 / matrix.n_features as f64
    };

    LfqRunSummary {
        n_detected,
        detection_rate,
        qvalue_percentiles: qp,
    }
}

pub fn build_align_report(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    matrix: &IntensityMatrix,
    config: &AlignConfig,
    timing: Timing,
) -> AlignReport {
    let t_outer = Instant::now();
    let min_combined = config.alignment.min_anchor_combined_score;

    // Each run's block is independent — parallelise across runs. The big
    // per-run cost is `build_anchors_block`, which does six sorts over the
    // ~70k anchors (for residual stats) plus the per-anchor record build.
    // Precompute `rt_range` for every run up front: `RunInput::rt_range`
    // falls back to a full feature sweep when `scan_times` is empty (which
    // it is in the align pipeline), so calling it inside the per-run closure
    // would be the same trap as before — only this time the alignment table
    // already stores it, so reuse that and fall back for the reference run.
    let rt_ranges: Vec<(f64, f64)> = runs
        .par_iter()
        .enumerate()
        .map(|(run_idx, ri)| {
            if run_idx == alignment.reference_idx {
                ri.rt_range()
            } else {
                alignment
                    .alignments
                    .get(&ri.name)
                    .map(|al| al.run_rt_range)
                    .unwrap_or_else(|| ri.rt_range())
            }
        })
        .collect();

    let runs_out: Vec<RunAlignReport> = runs
        .par_iter()
        .enumerate()
        .map(|(run_idx, ri)| {
            let is_reference = run_idx == alignment.reference_idx;
            let rt_range = rt_ranges[run_idx];
            let n_features = ri.features.len();
            let n_features_high_score = ri
                .features
                .iter()
                .filter(|f| f.combined_score >= min_combined)
                .count();
            // Hills are streamed per-run during LFQ (not held on `ri`), so the
            // count comes from the matrix, which records it as each run loads.
            let n_hills = matrix.hills_per_run.get(run_idx).copied().unwrap_or(0);

            let alignment_block = if is_reference {
                None
            } else {
                alignment.alignments.get(&ri.name).map(build_anchors_block)
            };

            let lfq = build_lfq_summary(matrix, run_idx);

            RunAlignReport {
                name: ri.name.clone(),
                is_reference,
                rt_range: [rt_range.0, rt_range.1],
                n_features,
                n_features_high_score,
                n_hills,
                alignment: alignment_block,
                lfq,
            }
        })
        .collect();

    let report = AlignReport {
        schema_version: 1,
        reference_run: alignment.reference_name.clone(),
        n_runs: runs.len(),
        n_consensus_features: matrix.n_features,
        config: config.clone(),
        runs: runs_out,
        timing,
    };
    log::info!("[timing] build_align_report inner: {:.2?}", t_outer.elapsed());
    report
}

// ── writer ───────────────────────────────────────────────────────────────────

pub fn write_align_report(report: &AlignReport, path: &Path) -> Result<(), KothError> {
    // Stream directly to disk via a buffered writer to avoid materialising the
    // full JSON in memory (~hundreds of MB on multi-run datasets) and to batch
    // the file I/O. Compact (non-pretty) output cuts the file ~3× by dropping
    // whitespace — the schema is consumed by `scripts/plot_align_report.py`
    // via `json.load`, which doesn't care about formatting.
    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(1 << 20, file); // 1 MiB buffer
    serde_json::to_writer(&mut writer, report)?;
    use std::io::Write;
    writer.flush()?;
    log::info!("Wrote align report to {}", path.display());
    Ok(())
}

