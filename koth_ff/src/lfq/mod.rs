pub mod consensus;
pub mod grid;
pub mod integrate;
pub mod score;
pub mod tdc;

use std::collections::HashMap;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::alignment::{AlignmentResult, RunInput};

use consensus::{build_consensus, ConsensusConfig, ConsensusFeature};
use grid::{build_grid, XicGrid};
use integrate::integrate;
use score::{score_grid, ColumnScores};

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
    /// Minimum spectral angle used as a peak-expansion stopping criterion
    pub spectral_angle_min: f64,
    /// Scoring mode used to find the best integration window
    pub score_mode: ScoreMode,
    /// Whether to run target-decoy competition and compute q-values
    pub run_tdc: bool,
    /// Superseded by `consensus.min_member_score`; kept for TOML backward compatibility.
    #[serde(default)]
    pub min_feature_score: f64,
    /// Tolerances for building the multi-run consensus feature list.
    #[serde(default)]
    pub consensus: ConsensusConfig,
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
            spectral_angle_min: 0.1,
            score_mode: ScoreMode::Hybrid,
            run_tdc: true,
            min_feature_score: 0.0,
            consensus: ConsensusConfig::default(),
        }
    }
}

/// Per-feature, per-run quantification result (internal).
#[derive(Debug, Clone)]
pub struct LfqEntry {
    pub feature_idx: usize,
    pub run_idx: usize,
    pub intensity: f64,
    pub hybrid_score: f32,
    pub spectral_angle: f32,
    pub n_isotopes_found: u8,
    pub is_decoy: bool,
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
    pub feature_score: Vec<f64>,
    /// Name of the run that provided the seed feature for each consensus row.
    pub feature_seed_run: Vec<String>,
    /// Number of runs that contributed at least one detection to each consensus group.
    pub feature_n_contributing_runs: Vec<usize>,
    // Flat [feature * n_runs + run] layout
    pub intensities: Vec<f64>,
    pub q_values: Vec<f64>,
}

impl IntensityMatrix {
    pub fn intensity(&self, feat: usize, run: usize) -> f64 {
        self.intensities[feat * self.n_runs + run]
    }

    pub fn q_value(&self, feat: usize, run: usize) -> f64 {
        self.q_values[feat * self.n_runs + run]
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
    let consensus: Vec<ConsensusFeature> =
        build_consensus(runs, alignment, config);

    let n_features = consensus.len();
    let n_runs = runs.len();

    log::info!(
        "LFQ: {} consensus features × {} runs",
        n_features, n_runs
    );

    // Pre-sort hill indices by mz for each run (binary search in build_grid)
    let hill_indices: Vec<Vec<usize>> = runs
        .iter()
        .map(|run| {
            let mut idx: Vec<usize> = (0..run.hills.len()).collect();
            idx.sort_by(|&a, &b| run.hills[a].mz.partial_cmp(&run.hills[b].mz).unwrap());
            idx
        })
        .collect();

    // Process every (feature, run) pair in parallel.
    // Scratch buffers are allocated once per rayon task (once per feature) and
    // reused across all n_runs iterations and the optional decoy pass.
    let all_entries: Vec<LfqEntry> = (0..n_features)
        .into_par_iter()
        .flat_map(|feat_idx| {
            let cf = &consensus[feat_idx];
            let mut entries = Vec::with_capacity(n_runs * 2);

            // Per-task scratch — allocated once, reused for every (run, target/decoy) pair.
            let mut grid = XicGrid::empty(config.n_isotopes, config.grid_cols, 0.0, 1.0);
            let mut scores = ColumnScores::new(config.grid_cols);
            let mut col_totals = vec![0.0f32; config.grid_cols];
            let mut obs = vec![0.0f64; config.n_isotopes];

            for run_idx in 0..n_runs {
                let run = &runs[run_idx];
                let run_rt_range = run.rt_range();
                let sorted_idx = &hill_indices[run_idx];

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
                build_grid(
                    &mut grid, run, sorted_idx, corr_mz, cf.charge,
                    corr_rt, corr_im, half_window, config,
                );
                let tgt_slots = grid.n_slots_filled;
                score_grid(&grid, &cf.theoretical_pattern, config, &mut scores, &mut col_totals, &mut obs);
                let tgt_peak = integrate(&grid, &scores, &col_totals, config);

                entries.push(LfqEntry {
                    feature_idx: feat_idx,
                    run_idx,
                    intensity: tgt_peak.intensity,
                    hybrid_score: tgt_peak.hybrid_score,
                    spectral_angle: tgt_peak.spectral_angle_at_apex,
                    n_isotopes_found: tgt_slots,
                    is_decoy: false,
                });

                // Decoy grid (shifted +11 Da in mz space, RT back 1% of gradient).
                // Reuses the same scratch buffers sequentially.
                if config.run_tdc {
                    let dec_mz = corr_mz + 11.0 / cf.charge as f64;
                    let dec_rt = corr_rt - 0.01 * (run_rt_range.1 - run_rt_range.0);
                    build_grid(
                        &mut grid, run, sorted_idx, dec_mz, cf.charge,
                        dec_rt, corr_im, half_window, config,
                    );
                    let dec_slots = grid.n_slots_filled;
                    score_grid(&grid, &cf.theoretical_pattern, config, &mut scores, &mut col_totals, &mut obs);
                    let dec_peak = integrate(&grid, &scores, &col_totals, config);

                    entries.push(LfqEntry {
                        feature_idx: feat_idx,
                        run_idx,
                        intensity: dec_peak.intensity,
                        hybrid_score: dec_peak.hybrid_score,
                        spectral_angle: dec_peak.spectral_angle_at_apex,
                        n_isotopes_found: dec_slots,
                        is_decoy: true,
                    });
                }
            }

            entries
        })
        .collect();

    // Target-decoy q-values
    let q_map: HashMap<(usize, usize), f64> = if config.run_tdc {
        tdc::compute_qvalues(&all_entries)
    } else {
        all_entries
            .iter()
            .filter(|e| !e.is_decoy)
            .map(|e| ((e.feature_idx, e.run_idx), 1.0))
            .collect()
    };

    // Per-run detection summary
    for run_idx in 0..n_runs {
        let n_detected = all_entries
            .iter()
            .filter(|e| !e.is_decoy && e.run_idx == run_idx && e.intensity > 0.0)
            .count();
        log::info!(
            "LFQ: run '{}': {}/{} features detected ({:.1}%)",
            runs[run_idx].name,
            n_detected,
            n_features,
            100.0 * n_detected as f64 / n_features as f64,
        );
    }

    // Assemble output matrix
    let mut intensities = vec![0.0f64; n_features * n_runs];
    let mut q_values_out = vec![1.0f64; n_features * n_runs];

    for entry in all_entries.iter().filter(|e| !e.is_decoy) {
        let flat = entry.feature_idx * n_runs + entry.run_idx;
        intensities[flat] = entry.intensity;
        if let Some(&qv) = q_map.get(&(entry.feature_idx, entry.run_idx)) {
            q_values_out[flat] = qv;
        }
    }

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
        feature_score: consensus.iter().map(|f| f.seed_score).collect(),
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
    }
}
