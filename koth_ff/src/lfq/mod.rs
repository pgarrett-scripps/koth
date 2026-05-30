pub mod consensus;
pub mod grid;
pub mod integrate;
pub mod score;
pub mod tdc;

use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::alignment::{AlignmentResult, RunInput};

use consensus::{build_consensus, ConsensusConfig, ConsensusFeature};
use grid::{build_grid, SortedHills, XicGrid};
use integrate::{integrate, PeakResult};
use score::{score_grid, ColumnScores};

/// Convert peak apex/start/end bins into native-space RT statistics.
/// Returns `(apex_rt, peak_width_rt)`; both NaN/0.0 when no peak was found.
fn peak_rt_stats(grid: &XicGrid, n_cols: usize, peak: &PeakResult) -> (f64, f64) {
    if peak.intensity <= 0.0 || n_cols == 0 {
        return (f64::NAN, 0.0);
    }
    let span = grid.rt_max - grid.rt_min;
    if span <= 0.0 {
        return (f64::NAN, 0.0);
    }
    let bin_w = span / n_cols as f64;
    let apex_rt = grid.rt_min + (peak.apex_bin as f64 + 0.5) * bin_w;
    // Inclusive bin range → +1 bin worth of width.
    let width = (peak.end_bin as f64 - peak.start_bin as f64 + 1.0) * bin_w;
    (apex_rt, width)
}

/// Configuration for LFQ grid extraction and peak integration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LfqConfig {
    /// m/z tolerance for hill lookup (ppm)
    pub mz_ppm: f64,
    /// Half-window size as a fraction of the run's total RT span.
    /// 0.005 = ±0.5% of gradient → 1% total window, centred on the feature RT.
    pub rt_window_pct: f64,
    /// Ion mobility tolerance for hill lookup (absolute 1/K0 units)
    pub im_tolerance: f64,
    /// Number of isotopologue rows: 1 = M only, 2 = M+M1, 3 = M+M1+M2
    pub n_isotopes: usize,
    /// Number of RT bins per grid (default 100)
    pub grid_cols: usize,
    /// Minimum spectral cosine (XIC column vs theoretical pattern) used as a
    /// peak-expansion stopping criterion. Despite the historical "spectral
    /// angle" name in the literature, the actual value computed here is a
    /// cosine similarity in [0, 1].
    pub spectral_cosine_min: f64,
    /// Scoring mode used to find the best integration window
    pub score_mode: ScoreMode,
    /// Whether to run target-decoy competition and compute q-values
    pub run_tdc: bool,
    /// Decoy m/z shift in Da, added to the target m/z and divided by charge
    /// (so the decoy m/z = target_mz + decoy_mz_shift_da / charge). Must be
    /// large enough to fall outside any plausible isotopologue or adduct
    /// envelope of the target.
    #[serde(default = "default_decoy_mz_shift_da")]
    pub decoy_mz_shift_da: f64,
    /// Decoy RT shift as a fraction of the run's RT span, subtracted from the
    /// target RT. With the default 0.01 and a ±`rt_window_pct` half-window,
    /// the decoy window still overlaps the target window — keep this within
    /// `rt_window_pct` if you want some overlap, or larger to fully separate.
    #[serde(default = "default_decoy_rt_shift_pct")]
    pub decoy_rt_shift_pct: f64,
    /// Tolerances for building the multi-run consensus feature list.
    #[serde(default)]
    pub consensus: ConsensusConfig,
}

fn default_decoy_mz_shift_da() -> f64 {
    11.0
}

fn default_decoy_rt_shift_pct() -> f64 {
    0.01
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ScoreMode {
    Hybrid,
    Rt,
    Intensity,
    Spectral,
}

impl Default for LfqConfig {
    fn default() -> Self {
        Self {
            mz_ppm: 10.0,
            rt_window_pct: 0.02,
            im_tolerance: 0.05,
            n_isotopes: 3,
            grid_cols: 100,
            spectral_cosine_min: 0.1,
            score_mode: ScoreMode::Hybrid,
            run_tdc: true,
            decoy_mz_shift_da: default_decoy_mz_shift_da(),
            decoy_rt_shift_pct: default_decoy_rt_shift_pct(),
            consensus: ConsensusConfig::default(),
        }
    }
}

/// Per-feature, per-run quantification result.
///
/// Holds enough information to write the long-format `lfq_details` file:
/// integration outputs (intensity / scores), RT diagnostics, and the
/// observed monoisotopic m/z + ion mobility of the winning M hill.
///
/// "Score" terminology here:
///   - `hybrid_score`    = cbrt(rt_score × intensity_score × spectral_cosine)
///                         at the apex column of the LFQ XIC grid
///   - `spectral_cosine` = cosine similarity (XIC column vs theoretical
///                         pattern) at the apex column. Different from
///                         `Feature.cosine_score` (chromatographic cosine
///                         between isotope hills, feature-finding stage).
#[derive(Debug, Clone)]
pub struct LfqEntry {
    pub feature_idx: usize,
    pub run_idx: usize,
    pub intensity: f64,
    pub hybrid_score: f32,
    /// Spectral cosine (XIC column vs theoretical pattern) at the apex column.
    /// Bounded [0, 1].
    pub spectral_cosine: f32,
    pub n_isotopes_found: u8,
    pub is_decoy: bool,
    /// Expected RT in native run space used to centre the XIC grid. For
    /// targets this is the alignment-predicted RT; for decoys it is the
    /// shifted decoy RT (`expected_rt = target_rt − decoy_rt_shift_pct · rt_span`).
    pub expected_rt: f64,
    /// Observed peak apex RT derived from the grid's apex bin.
    /// NaN when no peak was found.
    pub apex_rt: f64,
    /// Width of the integration window in RT units (end_rt - start_rt).
    /// 0.0 when no peak was found.
    pub peak_width_rt: f64,
    /// Expected monoisotopic m/z used to centre the XIC grid. For targets
    /// this is the alignment-predicted m/z; for decoys it is the shifted
    /// decoy m/z (`expected_mz = target_mz + decoy_mz_shift_da / charge`).
    /// Used as the reference for the per-row ppm residual.
    pub expected_mz: f64,
    /// Monoisotopic (M) m/z of the winning hill. NaN when the M row was empty.
    pub observed_mz: f64,
    /// Expected ion mobility used to centre the XIC grid. Decoys reuse the
    /// target IM today (no IM shift), but the field is stored explicitly so
    /// downstream residuals always have a per-row reference.
    pub expected_im: f64,
    /// Ion mobility of the winning M hill. NaN when the M row was empty
    /// or the run has no IM dimension.
    pub observed_im: f64,
}

/// Final intensity matrix: features (rows) × runs (columns).
pub struct IntensityMatrix {
    pub n_features: usize,
    pub n_runs: usize,
    pub run_names: Vec<String>,
    pub reference_run: String,
    // Per-consensus-feature metadata
    pub feature_mz: Vec<f64>,
    pub feature_mass: Vec<f64>,
    pub feature_charge: Vec<u8>,
    pub feature_rt: Vec<f64>,
    pub feature_im: Vec<f64>,
    /// Seed feature's `combined_score` (isotope × chromato cosine) per consensus row.
    pub feature_combined_score: Vec<f64>,
    /// Name of the run that provided the seed feature for each consensus row.
    pub feature_seed_run: Vec<String>,
    /// Number of runs that contributed at least one detection to each consensus group.
    pub feature_n_contributing_runs: Vec<usize>,
    // Flat [feature * n_runs + run] layout
    pub intensities: Vec<f64>,
    pub q_values: Vec<f64>,
    /// Decoy intensities, same flat layout as `intensities`.
    /// All zeros when `run_tdc = false`.
    pub decoy_intensities: Vec<f64>,
    /// Full per-(feature, run) target + decoy entries with detail stats.
    /// Empty when `run_tdc = false` and details are not requested.
    pub entries: Vec<LfqEntry>,
}

impl IntensityMatrix {
    pub fn intensity(&self, feat: usize, run: usize) -> f64 {
        self.intensities[feat * self.n_runs + run]
    }

    pub fn q_value(&self, feat: usize, run: usize) -> f64 {
        self.q_values[feat * self.n_runs + run]
    }

    pub fn decoy_intensity(&self, feat: usize, run: usize) -> f64 {
        self.decoy_intensities[feat * self.n_runs + run]
    }
}

/// Run LFQ quantification across all runs.
///
/// Consensus features are built from all runs via `build_consensus`, which
/// projects every feature into reference-run coordinate space and groups by
/// (charge, neutral mass ± ppm, aligned RT ± window, IM ± tolerance).
/// The highest-scoring feature in each group seeds the LFQ extraction.
/// For every consensus feature × run combination the function:
///   1. Applies the alignment corrections to the seed's reference-space coordinates.
///   2. Builds a 100-bin × n_isotopes XIC grid from the run's hills.
///   3. Scores each column and integrates the best peak.
///   4. Optionally repeats with a decoy feature for target-decoy q-values.
pub fn quantify(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
) -> IntensityMatrix {
    let t_consensus = Instant::now();
    let consensus: Vec<ConsensusFeature> =
        build_consensus(runs, alignment, config);
    log::info!("[timing] build_consensus: {:.2?}", t_consensus.elapsed());

    let n_features = consensus.len();
    let n_runs = runs.len();

    log::info!(
        "LFQ: {} consensus features × {} runs",
        n_features, n_runs
    );

    // Build a cache-friendly, mz-sorted struct-of-arrays per run.
    // Done in parallel across runs since each run is independent.
    let t_sort = Instant::now();
    let sorted_hills: Vec<SortedHills> =
        runs.par_iter().map(SortedHills::from_run).collect();
    log::info!(
        "[timing] sorted_hills build ({} runs, total {} hills): {:.2?}",
        n_runs,
        sorted_hills.iter().map(|s| s.len()).sum::<usize>(),
        t_sort.elapsed()
    );

    // Precompute each run's RT range once. `RunInput::rt_range` falls back to
    // a full sweep of `run.features` when `scan_times` is empty (which it is
    // for the align pipeline — feature/hill files don't carry scan_times),
    // so calling it inside the hot loop costs ~10ms × 7871 × 20 ≈ 100s.
    let t_rt = Instant::now();
    let run_rt_ranges: Vec<(f64, f64)> =
        runs.par_iter().map(|r| r.rt_range()).collect();
    log::info!(
        "[timing] run_rt_ranges precompute: {:.2?}",
        t_rt.elapsed()
    );

    // Process every (feature, run) pair in parallel.
    // Scratch buffers are allocated once per rayon task (once per feature) and
    // reused across all n_runs iterations and the optional decoy pass.
    let t_quant = Instant::now();
    let progress = AtomicUsize::new(0);
    let progress_step = (n_features / 20).max(1);
    // Per-phase nanosecond accumulators across all rayon tasks.
    let ns_build = AtomicU64::new(0);
    let ns_score = AtomicU64::new(0);
    let ns_integrate = AtomicU64::new(0);
    let all_entries: Vec<LfqEntry> = (0..n_features)
        .into_par_iter()
        .flat_map_iter(|feat_idx| {
            let cf = &consensus[feat_idx];
            let mut entries = Vec::with_capacity(n_runs * 2);

            // Per-task scratch — allocated once, reused for every (run, target/decoy) pair.
            let mut grid = XicGrid::empty(config.n_isotopes, config.grid_cols, 0.0, 1.0);
            let mut scores = ColumnScores::new(config.grid_cols);
            let mut col_totals = vec![0.0f32; config.grid_cols];
            let mut obs = vec![0.0f64; config.n_isotopes];

            for run_idx in 0..n_runs {
                let run = &runs[run_idx];
                let run_rt_range = run_rt_ranges[run_idx];
                let hills = &sorted_hills[run_idx];

                // Seed coordinates are already in reference-run space.
                // For the reference run use them directly; for other runs invert
                // the alignment warp to get expected native-space coordinates.
                let (corr_rt, corr_mz, corr_im) = if run_idx == alignment.reference_idx {
                    (cf.ref_rt, cf.ref_mz, cf.ref_im)
                } else {
                    let al = &alignment.alignments[&run.name];
                    (
                        al.predict_run_rt(cf.ref_rt),
                        al.predict_run_mz(cf.ref_mz, cf.ref_rt),
                        al.predict_run_im(cf.ref_im, cf.ref_rt),
                    )
                };

                let half_window = config.rt_window_pct * (run_rt_range.1 - run_rt_range.0);

                // Target
                let t_b = Instant::now();
                build_grid(
                    &mut grid, run, hills, corr_mz, cf.charge,
                    corr_rt, corr_im, half_window, config,
                );
                ns_build.fetch_add(t_b.elapsed().as_nanos() as u64, Ordering::Relaxed);
                let tgt_slots = grid.n_slots_filled;
                let tgt_obs_mz = grid.winner_mz[0];
                let tgt_obs_im = grid.winner_im[0];
                let t_s = Instant::now();
                score_grid(&grid, &cf.theoretical_pattern, config, &mut scores, &mut col_totals, &mut obs);
                ns_score.fetch_add(t_s.elapsed().as_nanos() as u64, Ordering::Relaxed);
                let t_i = Instant::now();
                let tgt_peak = integrate(&grid, &scores, &col_totals, config);
                ns_integrate.fetch_add(t_i.elapsed().as_nanos() as u64, Ordering::Relaxed);
                let (tgt_apex_rt, tgt_peak_width) =
                    peak_rt_stats(&grid, config.grid_cols, &tgt_peak);

                entries.push(LfqEntry {
                    feature_idx: feat_idx,
                    run_idx,
                    intensity: tgt_peak.intensity,
                    hybrid_score: tgt_peak.hybrid_score,
                    spectral_cosine: tgt_peak.spectral_cosine_at_apex,
                    n_isotopes_found: tgt_slots,
                    is_decoy: false,
                    expected_rt: corr_rt,
                    apex_rt: tgt_apex_rt,
                    peak_width_rt: tgt_peak_width,
                    expected_mz: corr_mz,
                    observed_mz: tgt_obs_mz,
                    expected_im: corr_im,
                    observed_im: tgt_obs_im,
                });

                // Decoy grid: shifted by `decoy_mz_shift_da / charge` in m/z
                // and by `−decoy_rt_shift_pct · rt_span` in RT.
                // Reuses the same scratch buffers sequentially.
                if config.run_tdc {
                    let dec_mz = corr_mz + config.decoy_mz_shift_da / cf.charge as f64;
                    let dec_rt = corr_rt
                        - config.decoy_rt_shift_pct * (run_rt_range.1 - run_rt_range.0);
                    let t_b = Instant::now();
                    build_grid(
                        &mut grid, run, hills, dec_mz, cf.charge,
                        dec_rt, corr_im, half_window, config,
                    );
                    ns_build.fetch_add(t_b.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    let dec_slots = grid.n_slots_filled;
                    let dec_obs_mz = grid.winner_mz[0];
                    let dec_obs_im = grid.winner_im[0];
                    let t_s = Instant::now();
                    score_grid(&grid, &cf.theoretical_pattern, config, &mut scores, &mut col_totals, &mut obs);
                    ns_score.fetch_add(t_s.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    let t_i = Instant::now();
                    let dec_peak = integrate(&grid, &scores, &col_totals, config);
                    ns_integrate.fetch_add(t_i.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    let (dec_apex_rt, dec_peak_width) =
                        peak_rt_stats(&grid, config.grid_cols, &dec_peak);

                    entries.push(LfqEntry {
                        feature_idx: feat_idx,
                        run_idx,
                        intensity: dec_peak.intensity,
                        hybrid_score: dec_peak.hybrid_score,
                        spectral_cosine: dec_peak.spectral_cosine_at_apex,
                        n_isotopes_found: dec_slots,
                        is_decoy: true,
                        expected_rt: dec_rt,
                        apex_rt: dec_apex_rt,
                        peak_width_rt: dec_peak_width,
                        expected_mz: dec_mz,
                        observed_mz: dec_obs_mz,
                        expected_im: corr_im,
                        observed_im: dec_obs_im,
                    });
                }
            }

            let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
            if done % progress_step == 0 || done == n_features {
                log::info!(
                    "[lfq] quantified {}/{} features ({:.0}%) in {:.2?}",
                    done,
                    n_features,
                    100.0 * done as f64 / n_features as f64,
                    t_quant.elapsed(),
                );
            }

            entries
        })
        .collect();
    log::info!(
        "[timing] quantification loop ({} entries): {:.2?}",
        all_entries.len(),
        t_quant.elapsed()
    );
    let total_ns = ns_build.load(Ordering::Relaxed)
        + ns_score.load(Ordering::Relaxed)
        + ns_integrate.load(Ordering::Relaxed);
    let pct = |ns: u64| {
        if total_ns == 0 { 0.0 } else { 100.0 * ns as f64 / total_ns as f64 }
    };
    let ns_b = ns_build.load(Ordering::Relaxed);
    let ns_s = ns_score.load(Ordering::Relaxed);
    let ns_i = ns_integrate.load(Ordering::Relaxed);
    log::info!(
        "[timing] hot-loop breakdown (cpu-ns across all threads): \
         build_grid={:.2}s ({:.0}%), score_grid={:.2}s ({:.0}%), integrate={:.2}s ({:.0}%)",
        ns_b as f64 / 1e9, pct(ns_b),
        ns_s as f64 / 1e9, pct(ns_s),
        ns_i as f64 / 1e9, pct(ns_i),
    );

    // Target-decoy q-values
    let t_tdc = Instant::now();
    let q_map: HashMap<(usize, usize), f64> = if config.run_tdc {
        tdc::compute_qvalues(&all_entries)
    } else {
        all_entries
            .iter()
            .filter(|e| !e.is_decoy)
            .map(|e| ((e.feature_idx, e.run_idx), 1.0))
            .collect()
    };
    log::info!("[timing] tdc q-values: {:.2?}", t_tdc.elapsed());

    // Per-run detection summary
    let t_summary = Instant::now();
    for run_idx in 0..n_runs {
        let n_detected = all_entries
            .iter()
            .filter(|e| !e.is_decoy && e.run_idx == run_idx && e.intensity > 0.0)
            .count();
        let n_decoy = all_entries
            .iter()
            .filter(|e| e.is_decoy && e.run_idx == run_idx && e.intensity > 0.0)
            .count();
        if config.run_tdc {
            log::info!(
                "LFQ: run '{}': {}/{} target detected ({:.1}%), {} decoy",
                runs[run_idx].name,
                n_detected,
                n_features,
                100.0 * n_detected as f64 / n_features as f64,
                n_decoy,
            );
        } else {
            log::info!(
                "LFQ: run '{}': {}/{} features detected ({:.1}%)",
                runs[run_idx].name,
                n_detected,
                n_features,
                100.0 * n_detected as f64 / n_features as f64,
            );
        }
    }

    log::info!("[timing] per-run detection summary: {:.2?}", t_summary.elapsed());

    // Assemble output matrix
    let t_assemble = Instant::now();
    let mut intensities = vec![0.0f64; n_features * n_runs];
    let mut q_values_out = vec![1.0f64; n_features * n_runs];
    let mut decoy_intensities = vec![0.0f64; n_features * n_runs];

    for entry in &all_entries {
        let flat = entry.feature_idx * n_runs + entry.run_idx;
        if entry.is_decoy {
            decoy_intensities[flat] = entry.intensity;
        } else {
            intensities[flat] = entry.intensity;
            if let Some(&qv) = q_map.get(&(entry.feature_idx, entry.run_idx)) {
                q_values_out[flat] = qv;
            }
        }
    }

    log::info!("[timing] intensity matrix assemble: {:.2?}", t_assemble.elapsed());

    IntensityMatrix {
        n_features,
        n_runs,
        run_names: runs.iter().map(|r| r.name.clone()).collect(),
        reference_run: alignment.reference_name.clone(),
        feature_mz: consensus.iter().map(|f| f.ref_mz).collect(),
        feature_mass: consensus.iter().map(|f| f.neutral_mass).collect(),
        feature_charge: consensus.iter().map(|f| f.charge).collect(),
        feature_rt: consensus.iter().map(|f| f.ref_rt).collect(),
        feature_im: consensus.iter().map(|f| f.ref_im).collect(),
        feature_combined_score: consensus.iter().map(|f| f.seed_combined_score).collect(),
        feature_seed_run: consensus
            .iter()
            .map(|f| runs[f.seed_run_idx].name.clone())
            .collect(),
        feature_n_contributing_runs: consensus
            .iter()
            .map(|f| f.n_contributing_runs)
            .collect(),
        intensities,
        q_values: q_values_out,
        decoy_intensities,
        entries: all_entries,
    }
}
