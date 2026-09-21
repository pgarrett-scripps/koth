pub mod consensus;
pub mod grid;
pub mod integrate;
mod ownership;
pub mod rescore;
pub mod score;
pub mod targets;
pub mod tdc;

use crate::scoring::model::IsotopeModelSpec;
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::alignment::{AlignmentResult, RunInput};

use consensus::{build_consensus, ConsensusConfig, ConsensusFeature};
use grid::{build_grid_excluding, SampleId, SortedHills, XicGrid};
use integrate::{integrate, PeakResult};
use ownership::CellCandidate;
use score::{score_grid, ColumnScores};
use std::collections::HashSet;

use crate::models::Hill;
use crate::scoring::elements::K_PATTERN;

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
#[serde(default, deny_unknown_fields)]
pub struct LfqConfig {
    /// m/z tolerance for hill lookup (ppm)
    pub mz_ppm: f64,
    /// Half-window size as a fraction of the run's total RT span.
    /// 0.005 = ±0.5% of the observed RT span → 1% total window, centred on
    /// the alignment-predicted native RT. Default 0.005.
    pub rt_window_pct: f64,
    /// Ion mobility half-window for hill lookup (absolute 1/K0 units).
    /// Default 0.015; inert when the input has no ion mobility.
    pub im_tolerance: f64,
    /// Number of isotopologue rows: 1 = M only, 2 = M+M1, 3 = M+M1+M2
    pub n_isotopes: usize,
    /// Number of RT bins per grid (default 100)
    pub grid_cols: usize,
    /// Minimum spectral **Bhattacharyya** score for a grid column to keep
    /// extending the integration peak (peak-expansion gate). Named for the
    /// metric it actually uses — it is NOT a cosine threshold. Optional
    /// (defaults to 0.1 when omitted); the config was previously named
    /// `spectral_cosine_min`, which now errors (see `deny_unknown_fields`).
    #[serde(default = "default_min_spectral_bhattacharyya")]
    pub min_spectral_bhattacharyya: f64,
    /// Scoring mode used to find the best integration window
    pub score_mode: ScoreMode,
    // NOTE on isotope scoring (this used to be a source of confusion): there
    // are two DISTINCT isotope signals, both ALWAYS computed, no config flags:
    //   * spectral Bhattacharyya — observed vs theoretical isotope PATTERN.
    //   * co-elution cosine       — cosine between the matched isotopes' XIC
    //                               traces across RT (do they co-elute?).
    // hybrid = (rt · intensity · bhattacharyya · coelution)^¼. The old
    // `spectral_bhattacharyya` and `spectral_coelution` flags are gone;
    // leftover values in an existing TOML now error (see `deny_unknown_fields`).
    /// Whether to run target-decoy competition and compute q-values
    pub run_tdc: bool,
    /// Decoy m/z shift in Da, added to the target m/z and divided by charge
    /// (so the decoy m/z = target_mz + decoy_mz_shift_da / charge). Must be
    /// large enough to fall outside any plausible isotopologue or adduct
    /// envelope of the target.
    #[serde(default = "default_decoy_mz_shift_da")]
    pub decoy_mz_shift_da: f64,
    /// Decoy RT shift as a fraction of the run's RT span, subtracted from the
    /// target RT. The default 0.01 is twice the default extraction half-window,
    /// so the target and decoy RT windows meet at one boundary. A shift larger
    /// than twice `rt_window_pct` fully separates their RT intervals.
    #[serde(default = "default_decoy_rt_shift_pct")]
    pub decoy_rt_shift_pct: f64,
    /// Tolerances for building the multi-run consensus feature list.
    #[serde(default)]
    pub consensus: ConsensusConfig,
    /// Cross-run intensity normalization applied to the assembled matrix.
    /// "none"           = raw intensities (original; downstream must normalize).
    /// "median_ratios"  = DESeq/edgeR-style size factors: each run is scaled by
    ///                    the median, over features detected in ALL runs, of its
    ///                    log2 deviation from the per-feature mean. Robust to a
    ///                    fraction of genuinely-changing features, unlike a plain
    ///                    column-median (which is biased when a sizeable subset of
    ///                    peptides changes systematically). Default "none".
    #[serde(default = "default_normalize")]
    pub normalize: String,
    /// Per-cell intensity estimator for the consensus matrix.
    /// "sum"  = integrated peak area (detected: feature `intensitySum`; MBR:
    ///          summed XIC grid over `start..=end`). Original behaviour.
    /// "apex" = peak height (detected: feature `intensityApex`; MBR: total
    ///          isotopologue intensity at the apex grid column). Removes the
    ///          integration-window variance that inflates replicate CV on
    ///          MBR-filled cells, matching the per-run finder's apex estimator.
    #[serde(default = "default_quant_estimator")]
    pub quant_estimator: String,
    /// Quantify EVERY consensus cell (detected and MBR) by the same grid
    /// re-integration, instead of using the per-run feature's own intensity for
    /// detected cells. Keeps detected and MBR cells on one commensurate scale —
    /// essential on timsTOF, where the feature integrates the ion-mobility
    /// dimension but the 2-D XIC grid does not, so mixing the two inflates CV,
    /// and validated as a win on Orbitrap too (PXD003881: gated matrix CV
    /// 14.27 → 12.31 %, ECOLI bias −0.067 → −0.020, HUMAN IQR 0.219 → 0.194).
    /// Default TRUE since 2026-08-27: one estimator, one scale, everywhere.
    /// Legacy false values now warn and use grid intensities: exclusive native
    /// signal ownership cannot validate a feature-intensity override.
    #[serde(default = "default_detected_use_grid")]
    pub detected_use_grid: bool,
    /// Replace the raw RT closeness term in the hybrid score with a
    /// σ-normalised Gaussian RT likelihood `exp(−½ z²)`, where
    /// `z = (column RT − predicted RT) / σ_local(RT)` and `σ_local` is the
    /// post-warp RT-residual spread estimated per run from the RANSAC inlier
    /// anchors (see `alignment::RtSigmaModel`). Makes the RT term of the
    /// target-decoy discriminant region-aware: strict where alignment is
    /// confident, tolerant where it is not. Applied identically to target and
    /// decoy cells so the TDC null stays calibrated. Non-reference runs only;
    /// falls back to the raw term when no σ model is available. Default false.
    #[serde(default)]
    pub rt_spread_scoring: bool,
    /// Experimental search-guided scorer; legacy preserves the release path.
    #[serde(default)]
    pub search_scoring: rescore::search::SearchScoring,
    /// Rescue unambiguous same-run IDs when global RT observations conflict.
    /// Ambiguous targets remain ineligible for cross-run transfer.
    #[serde(default)]
    pub search_rt_rescue: bool,
    /// Decide transfer eligibility per recipient run rather than per target.
    /// Without this, pruning one run's ambiguous IDs withdraws the peptide from
    /// transfer in every run, including runs whose chromatography was clean.
    /// Requires `search_rt_rescue`; runs whose own IDs were withheld as
    /// unsupported stay ineligible, and a target whose retained IDs still
    /// disagree across runs remains globally ineligible because its consensus
    /// coordinate is what a transfer would be centred on.
    #[serde(default)]
    pub search_rt_rescue_per_run_transfers: bool,
    /// Experimental inferred 2+/3+/4+ targets, with explicit charge provenance.
    #[serde(default)]
    pub search_expand_charges: bool,
    /// Apply the extraction q-value gate only to transferred cells. A cell with
    /// an accepted same-run MS2 identification has already passed peptide-level
    /// FDR; gating it again on an exploratory extraction score re-litigates an
    /// identification that was accepted, and discards a quantity for a peptide
    /// known to be present in that run. Transfers, which have no such evidence,
    /// stay gated. FlashLFQ applies FDR to its match-between-runs peaks only.
    /// Identification-free mode is unaffected: every cell there is inferred, so
    /// every cell stays gated.
    #[serde(default)]
    pub search_gate_transfers_only: bool,
    /// Report the *averagine-projected* intensity per cell instead of the raw
    /// box-sum. For each grid column the observed isotopologue vector is passed
    /// through a matched filter for the theoretical averagine pattern (the same
    /// mass-derived fingerprint the Bhattacharyya score already uses):
    /// `proj = <observed, pattern_hat>`, where `pattern_hat` is the L2-normalised
    /// pattern. As a matched filter this maximises SNR for on-pattern signal
    /// (recovering faint peptides the raw sum buries in noise) and rejects the
    /// component of each column orthogonal to the fingerprint (intensity
    /// distributed unlike averagine — interference, chemical noise). It is NOT a
    /// pattern-match gate: gross mismatches (lone monoisotopes, wrong charge)
    /// are still rejected by the Bhattacharyya window gate and the TDC q-value,
    /// unchanged. The per-peptide constant `‖pattern‖` factor cancels in
    /// cross-run ratios and CV, so LFQ ratios/CV are unaffected in scale. Apex
    /// selection and window expansion stay on the raw signal / hybrid score;
    /// only the reported value changes. Default false (byte-identical raw
    /// box-sum). Mass-based, ID-free. Effect on the benchmark is unmeasured —
    /// this is an opt-in knob to be validated, not a new default.
    #[serde(default)]
    pub averagine_projection: bool,
    /// (Experiment B / audit A2) Score each **decoy** cell against an averagine
    /// isotope template and theoretical pattern computed from the DECOY's own
    /// shifted mass (`neutral_mass + decoy_mz_shift_da`) rather than the target's
    /// mass. When false the decoy reuses the target's template. When true
    /// (the default) the +11 Da decoy is judged against the isotope envelope it
    /// actually sits on, removing the target-template freebie a mis-massed
    /// decoy otherwise inherits on the Bhattacharyya (QDA feature 2) and the
    /// co-elution theory weighting (feature 3). Applies only to the TDC decoy
    /// grid; target scoring and all reported intensities are unchanged.
    #[serde(default = "default_decoy_own_template")]
    pub decoy_own_template: bool,
    /// (Experiment B / audit A4) Co-elution value assigned to a grid cell with
    /// fewer than two isotope rows carrying signal (a lone monoisotope — nothing
    /// to co-elute). Applied identically to target and decoy cells. Default 0.5
    /// (neutral). 1.0 gives an un-judgeable cell the maximal, target-like
    /// score, which hands noise-grabbing lone-hill decoys a free target-like
    /// coordinate on QDA feature index 3.
    #[serde(default = "default_lone_coelution")]
    pub lone_coelution: f64,
    /// Which analyte class's average composition the LFQ isotope templates are
    /// built from. Mirrors `[features] isotope_model` and must match the model
    /// the per-run features were detected with, or the consensus cells are
    /// scored against a pattern the detector never used. Default `"peptide"`.
    #[serde(default)]
    pub isotope_model: IsotopeModelSpec,
}

fn default_min_spectral_bhattacharyya() -> f64 {
    0.1
}

fn default_lone_coelution() -> f64 {
    0.5
}

fn default_decoy_own_template() -> bool {
    true
}

fn default_normalize() -> String {
    "none".to_string()
}

fn default_quant_estimator() -> String {
    "sum".to_string()
}

fn default_detected_use_grid() -> bool {
    true
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
            rt_window_pct: 0.005,
            im_tolerance: 0.015,
            n_isotopes: 3,
            grid_cols: 100,
            min_spectral_bhattacharyya: 0.1,
            score_mode: ScoreMode::Hybrid,
            run_tdc: true,
            decoy_mz_shift_da: default_decoy_mz_shift_da(),
            decoy_rt_shift_pct: default_decoy_rt_shift_pct(),
            consensus: ConsensusConfig::default(),
            normalize: default_normalize(),
            quant_estimator: default_quant_estimator(),
            detected_use_grid: default_detected_use_grid(),
            rt_spread_scoring: false,
            search_scoring: rescore::search::SearchScoring::Legacy,
            search_rt_rescue: false,
            search_rt_rescue_per_run_transfers: false,
            search_expand_charges: false,
            search_gate_transfers_only: false,
            averagine_projection: false,
            decoy_own_template: default_decoy_own_template(),
            lone_coelution: default_lone_coelution(),
            isotope_model: IsotopeModelSpec::default(),
        }
    }
}

#[cfg(test)]
mod norm_tests {
    use super::apply_median_ratio_normalization;

    /// A uniform per-run scale factor must be fully removed (columns equalize
    /// for every complete-case feature).
    #[test]
    fn median_ratio_removes_uniform_scale() {
        let (n_features, n_runs) = (3usize, 3usize);
        // run 1 is uniformly 2x brighter than runs 0 and 2.
        let mut t = vec![
            100.0, 200.0, 100.0, //
            300.0, 600.0, 300.0, //
            50.0, 100.0, 50.0,
        ];
        let mut d = vec![0.0f64; n_features * n_runs];
        let sf = apply_median_ratio_normalization(&mut t, &mut d, n_features, n_runs);
        assert!(
            (sf[1] - 1.0).abs() < 1e-9,
            "run 1 size factor should be +1 log2, got {}",
            sf[1]
        );
        for feat in 0..n_features {
            let b = feat * n_runs;
            assert!(
                (t[b] - t[b + 1]).abs() < 1e-6 && (t[b] - t[b + 2]).abs() < 1e-6,
                "columns should equalize after removing the 2x scale"
            );
        }
    }

    /// Features missing in some run (0 intensity) are excluded from the size
    /// factor (complete-case only) and left untouched if zero.
    #[test]
    fn median_ratio_ignores_incomplete_features() {
        let (n_features, n_runs) = (2usize, 2usize);
        let mut t = vec![100.0, 100.0, 500.0, 0.0]; // 2nd feature missing in run 1
        let mut d = vec![0.0f64; 4];
        let sf = apply_median_ratio_normalization(&mut t, &mut d, n_features, n_runs);
        // only the complete feature (equal) drives the factor -> no scaling
        assert!(sf[0].abs() < 1e-9 && sf[1].abs() < 1e-9);
        assert_eq!(t[3], 0.0, "missing cell stays zero");
    }
}

/// Compute DESeq/edgeR-style size factors from the assembled target matrix and
/// apply them in place to both the target and decoy intensities.
///
/// A run's size factor is the median, over features present (> 0) in EVERY run,
/// of `log2(intensity) - mean_over_runs(log2(intensity))`. Factors are centred
/// so the median run is unscaled (overall intensity scale preserved). Returns
/// the centred per-run log2 size factors for logging.
///
/// Robust to a fraction of genuinely-changing features because the median
/// ignores the tails where real biology lives — unlike a plain column-median,
/// which shifts with the changing subset.
fn apply_median_ratio_normalization(
    intensities: &mut [f64],
    decoy_intensities: &mut [f64],
    n_features: usize,
    n_runs: usize,
) -> Vec<f64> {
    if n_runs < 2 {
        return vec![0.0; n_runs];
    }
    let mut ratios_per_run: Vec<Vec<f64>> = vec![Vec::new(); n_runs];
    for feat in 0..n_features {
        let base = feat * n_runs;
        let row = &intensities[base..base + n_runs];
        if row.iter().any(|&v| v <= 0.0) {
            continue; // complete cases only
        }
        let mean_log = row.iter().map(|&v| v.log2()).sum::<f64>() / n_runs as f64;
        for (run, &v) in row.iter().enumerate() {
            ratios_per_run[run].push(v.log2() - mean_log);
        }
    }
    let mut size_factors: Vec<f64> = ratios_per_run
        .iter_mut()
        .map(|v| crate::stats::median(v))
        .collect();
    // Centre so the median run is unscaled.
    let mut sf_copy = size_factors.clone();
    let center = crate::stats::median(&mut sf_copy);
    for s in &mut size_factors {
        *s -= center;
    }
    // Apply to both target and decoy columns: intensity /= 2^size_factor[run].
    for arr in [intensities, decoy_intensities] {
        for feat in 0..n_features {
            let base = feat * n_runs;
            for run in 0..n_runs {
                let v = arr[base + run];
                if v > 0.0 {
                    arr[base + run] = v / 2f64.powf(size_factors[run]);
                }
            }
        }
    }
    size_factors
}

/// Per-feature, per-run quantification result.
///
/// Holds enough information to write the long-format `lfq_details` file:
/// integration outputs (intensity / scores), RT diagnostics, and the
/// observed monoisotopic m/z + ion mobility of the winning M hill.
///
/// Two isotope-quality signals, kept strictly distinct (this is the pair that
/// was historically conflated):
///
/// - `spectral_bhattacharyya` = observed vs theoretical isotope abundance
///   pattern (penalises missing peaks).
/// - `coelution` = cosine between the matched isotopes' XIC traces across RT
///   ("do they co-elute?").
///
/// hybrid_score = (rt_score × int_score × spectral_bhattacharyya × coelution)^¼
/// at the apex column. Neither is `Feature.cosine_score` (the feature-finder's
/// chromatographic cosine between isotope hills).
#[derive(Debug, Clone)]
pub struct LfqEntry {
    /// Exclusive raw-sample accounting, separately for targets and controls.
    pub owned_samples: usize,
    /// Co-eluting preceding-isotope signal unexplained by this mono hypothesis.
    pub preceding_signal_fraction: f32,
    pub excluded_samples: usize,
    pub competing_feature: Option<usize>,
    pub ownership_status: &'static str,
    pub feature_idx: usize,
    pub run_idx: usize,
    pub intensity: f64,
    pub hybrid_score: f32,
    /// Spectral **Bhattacharyya** at the apex column — observed-vs-theoretical
    /// isotope pattern match; penalises expected-but-missing peaks. [0, 1].
    pub spectral_bhattacharyya: f32,
    pub n_isotopes_found: u8,
    /// Individual hybrid components at the apex column, exposed so a downstream
    /// rescorer can weight each independently instead of using only the
    /// composite `hybrid_score`. `coelution` is the inter-isotope cosine
    /// (do the matched isotopes co-elute). All bounded [0, 1].
    pub rt_score: f32,
    pub int_score: f32,
    pub coelution: f32,
    pub is_decoy: bool,
    /// True when this run did NOT contribute a feature to the consensus group,
    /// so the intensity was re-integrated from raw hills at the predicted RT/mz
    /// (match-between-runs). False = an original detection contributed to the
    /// group; both kinds are quantified by the same exclusive grid extraction.
    /// In search-guided mode this instead means no accepted same-run MS2 ID.
    pub is_mbr: bool,
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
    /// Peptide identities and donor evidence, present only in search-guided mode.
    pub search_guidance: Option<targets::SearchGuidance>,
    /// Original consensus membership retained for long-format evidence export.
    pub consensus: Vec<consensus::ConsensusFeature>,
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
    /// Hill count per run (recorded as each run's hills were streamed in).
    pub hills_per_run: Vec<usize>,
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
///
/// 1. Applies the alignment corrections to the seed's reference-space coordinates.
/// 2. Builds a 100-bin × n_isotopes XIC grid from the run's hills.
/// 3. Scores each column and integrates the best peak.
/// 4. Optionally repeats with a decoy feature for target-decoy q-values.
///
/// `runs` carry only features (their `hills` may be empty); the full hills for
/// each run are loaded on demand.
///
/// CPU-nanosecond accumulators for the three hot-loop phases, summed across all
/// rayon tasks. Kept as a small struct so `quantify_cell` can be handed one
/// reference instead of three loose atomics.
struct PhaseTimers {
    build: AtomicU64,
    score: AtomicU64,
    integrate: AtomicU64,
}

impl PhaseTimers {
    fn new() -> Self {
        Self {
            build: AtomicU64::new(0),
            score: AtomicU64::new(0),
            integrate: AtomicU64::new(0),
        }
    }

    fn add_build(&self, d: std::time::Duration) {
        self.build.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }
    fn add_score(&self, d: std::time::Duration) {
        self.score.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }
    fn add_integrate(&self, d: std::time::Duration) {
        self.integrate
            .fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }
}

/// Quantify a single (feature, run) cell — target OR decoy — by the shared
/// build_grid → score_grid → integrate → peak_rt_stats → `LfqEntry` sequence.
/// Grid intensity (apex or summed, per `quant_estimator`) is always used so
/// every reported quantity has original-sample provenance.
#[allow(clippy::too_many_arguments)]
fn quantify_cell(
    grid: &mut XicGrid,
    scores: &mut ColumnScores,
    col_totals: &mut [f32],
    obs: &mut [f64],
    hills_vec: &[Hill],
    sorted: &SortedHills,
    scan_times: &[f64],
    theoretical_pattern: &[f64],
    bc_template: &[f64; K_PATTERN],
    config: &LfqConfig,
    timers: &PhaseTimers,
    feat_idx: usize,
    run_idx: usize,
    charge: u8,
    mz: f64,
    rt: f64,
    im: f64,
    half_window: f64,
    rt_sigma: Option<f32>,
    is_decoy: bool,
    is_mbr: bool,
    excluded: &HashSet<SampleId>,
) -> CellCandidate {
    let t_b = Instant::now();
    build_grid_excluding(
        grid,
        hills_vec,
        sorted,
        scan_times,
        mz,
        charge,
        rt,
        im,
        half_window,
        config,
        excluded,
    );
    timers.add_build(t_b.elapsed());
    let slots = grid.n_slots_filled;
    let obs_mz = grid.winner_mz[0];
    let obs_im = grid.winner_im[0];

    let t_s = Instant::now();
    score_grid(
        grid,
        theoretical_pattern,
        bc_template,
        config,
        scores,
        col_totals,
        obs,
        rt_sigma,
    );
    timers.add_score(t_s.elapsed());

    let t_i = Instant::now();
    let peak = integrate(grid, scores, col_totals, config);
    timers.add_integrate(t_i.elapsed());
    let (apex_rt, peak_width) = peak_rt_stats(grid, config.grid_cols, &peak);

    let use_apex = config.quant_estimator.eq_ignore_ascii_case("apex");
    let intensity = if use_apex {
        peak.apex_intensity
    } else {
        peak.intensity
    };

    let mut support: Vec<_> = grid
        .samples
        .iter()
        .filter(|s| peak.intensity > 0.0 && s.col >= peak.start_bin && s.col <= peak.end_bin)
        .map(|s| s.id)
        .collect();
    support.sort_unstable();
    support.dedup();
    let mut rows = vec![false; config.n_isotopes];
    for s in &grid.samples {
        if peak.intensity > 0.0 && s.col >= peak.start_bin && s.col <= peak.end_bin {
            rows[s.row] = true;
        }
    }
    // Compare mono hypotheses against the preceding isotope position as well.
    // An M+1 alias may fit its truncated right-hand envelope while leaving a
    // substantial co-eluting M peak unexplained. This is ranking evidence only;
    // the preceding row is never included in this candidate's reported quantity.
    let preceding_signal_fraction = if peak.intensity > 0.0 {
        let mut preceding = XicGrid::empty(1, config.grid_cols, 0.0, 1.0);
        let mut one_row = config.clone();
        one_row.n_isotopes = 1;
        build_grid_excluding(
            &mut preceding,
            hills_vec,
            sorted,
            scan_times,
            mz - grid::C13_NEUTRON / charge as f64,
            charge,
            rt,
            im,
            half_window,
            &one_row,
            excluded,
        );
        let (mut dot, mut norm_a, mut norm_b) = (0.0f64, 0.0f64, 0.0f64);
        for (&a, &b) in preceding.intensities[0].iter().zip(&grid.intensities[0]) {
            let (a, b) = (f64::from(a), f64::from(b));
            dot += a * b;
            norm_a += a * a;
            norm_b += b * b;
        }
        let coelution = if norm_a > 0.0 && norm_b > 0.0 {
            (dot / (norm_a * norm_b).sqrt()).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let previous: f64 = preceding.intensities[0][peak.start_bin..=peak.end_bin]
            .iter()
            .map(|&v| f64::from(v))
            .sum();
        let positive: f64 = col_totals[peak.start_bin..=peak.end_bin]
            .iter()
            .map(|&v| f64::from(v))
            .sum();
        (coelution * previous / (previous + positive).max(f64::MIN_POSITIVE)) as f32
    } else {
        0.0
    };
    let bin_width = (grid.rt_max - grid.rt_min) / config.grid_cols as f64;
    CellCandidate {
        peak_start: grid.rt_min + peak.start_bin as f64 * bin_width,
        peak_end: grid.rt_min + (peak.end_bin + 1) as f64 * bin_width,
        peak_rows: rows.iter().filter(|&&v| v).count(),
        entry: LfqEntry {
            owned_samples: support.len(),
            preceding_signal_fraction,
            excluded_samples: 0,
            competing_feature: None,
            ownership_status: "exclusive",

            feature_idx: feat_idx,
            run_idx,
            intensity,
            hybrid_score: peak.hybrid_score,
            spectral_bhattacharyya: peak.bhattacharyya_at_apex,
            n_isotopes_found: slots,
            rt_score: peak.rt_score_at_apex,
            int_score: peak.int_score_at_apex,
            coelution: peak.coelution,
            is_decoy,
            is_mbr,
            expected_rt: rt,
            apex_rt,
            peak_width_rt: peak_width,
            expected_mz: mz,
            observed_mz: obs_mz,
            expected_im: im,
            observed_im: obs_im,
        },
        context: {
            let mut ids: Vec<_> = grid.samples.iter().map(|s| s.id).collect();
            ids.sort_unstable();
            ids.dedup();
            ids
        },
        support,
    }
}

/// Log the per-run detection summary (target/decoy counts, or plain detected
/// counts when TDC is off). Pure logging — no effect on the output matrix.
fn log_detection_summary(
    all_entries: &[LfqEntry],
    runs: &[RunInput],
    n_features: usize,
    n_runs: usize,
    run_tdc: bool,
) {
    for run_idx in 0..n_runs {
        let n_detected = all_entries
            .iter()
            .filter(|e| !e.is_decoy && e.run_idx == run_idx && e.intensity > 0.0)
            .count();
        let n_decoy = all_entries
            .iter()
            .filter(|e| e.is_decoy && e.run_idx == run_idx && e.intensity > 0.0)
            .count();
        if run_tdc {
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
}

/// Fold the per-(feature, run) target/decoy `entries` into the flat intensity /
/// q-value / decoy-intensity matrices, apply optional median-of-ratios
/// normalization, and build the final `IntensityMatrix`. Takes `all_entries` by
/// value — it is stored in the returned matrix.
#[allow(clippy::too_many_arguments)]
fn assemble_matrix(
    consensus: &[ConsensusFeature],
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
    all_entries: Vec<LfqEntry>,
    q_map: &HashMap<(usize, usize), f64>,
    hills_per_run: Vec<usize>,
    n_features: usize,
    n_runs: usize,
) -> IntensityMatrix {
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

    log::info!(
        "[timing] intensity matrix assemble: {:.2?}",
        t_assemble.elapsed()
    );

    if config.normalize.eq_ignore_ascii_case("median_ratios") {
        let sf = apply_median_ratio_normalization(
            &mut intensities,
            &mut decoy_intensities,
            n_features,
            n_runs,
        );
        let lo = sf.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = sf.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        log::info!(
            "LFQ normalization: median-of-ratios size factors applied, \
             per-run log2 range [{:.3}, {:.3}] (~{:.1}% max scale)",
            lo,
            hi,
            (2f64.powf(hi - lo) - 1.0) * 100.0
        );
    }

    IntensityMatrix {
        search_guidance: None,
        consensus: consensus.to_vec(),
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
        feature_n_contributing_runs: consensus.iter().map(|f| f.n_contributing_runs).collect(),
        hills_per_run,
        intensities,
        q_values: q_values_out,
        decoy_intensities,
        entries: all_entries,
    }
}

pub fn quantify(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
    load_hills: impl Fn(usize) -> Vec<crate::models::Hill>,
) -> IntensityMatrix {
    let t_consensus = Instant::now();
    let consensus = build_consensus(runs, alignment, config);
    log::info!("[timing] build_consensus: {:.2?}", t_consensus.elapsed());
    quantify_consensus(runs, alignment, config, consensus, load_hills)
}

/// Quantify an already selected candidate set using the shared exclusive
/// extraction and cell-rescoring pipeline. Used by controlled grouping comparisons.
pub fn quantify_consensus(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
    consensus: Vec<ConsensusFeature>,
    load_hills: impl Fn(usize) -> Vec<crate::models::Hill>,
) -> IntensityMatrix {
    quantify_candidates(runs, alignment, config, consensus, None, load_hills)
}

/// Quantify imported peptide targets with the same signal ownership and estimator
/// as identification-free LFQ. The guidance and consensus rows must correspond.
pub fn quantify_guided(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
    consensus: Vec<ConsensusFeature>,
    guidance: targets::SearchGuidance,
    load_hills: impl Fn(usize) -> Vec<crate::models::Hill>,
) -> IntensityMatrix {
    assert_eq!(consensus.len(), guidance.targets.len());
    let mut matrix = quantify_candidates(
        runs,
        alignment,
        config,
        consensus,
        Some(&guidance),
        load_hills,
    );
    matrix.search_guidance = Some(guidance);
    matrix
}

fn quantify_candidates(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
    consensus: Vec<ConsensusFeature>,
    guidance: Option<&targets::SearchGuidance>,
    load_hills: impl Fn(usize) -> Vec<crate::models::Hill>,
) -> IntensityMatrix {
    // Guard against a degenerate grid: grid_cols / n_isotopes of 0 would make
    // `n_cols - 1` underflow and index into empty rows (panic). Clamp to 1 and
    // warn rather than crash on a misconfigured TOML.
    let config = &if config.grid_cols == 0 || config.n_isotopes == 0 {
        log::warn!(
            "LFQ config has grid_cols={} n_isotopes={}; clamping each to a minimum of 1",
            config.grid_cols,
            config.n_isotopes
        );
        let mut c = config.clone();
        c.grid_cols = c.grid_cols.max(1);
        c.n_isotopes = c.n_isotopes.max(1);
        c
    } else {
        config.clone()
    };

    let n_features = consensus.len();
    let n_runs = runs.len();

    log::info!("LFQ: {} consensus features × {} runs", n_features, n_runs);

    // Full-length averagine Bhattacharyya templates, one per consensus feature.
    // Each depends only on `cf.neutral_mass`, so precompute here instead of
    // rebuilding inside every (run × feature × target/decoy) score_grid call.
    let isotope_model = config.isotope_model.model();
    let bc_templates: Vec<_> = consensus
        .iter()
        .map(|cf| isotope_model.distribution(cf.neutral_mass))
        .collect();

    // (Audit A2) Decoy-specific averagine template + theoretical pattern, computed
    // from the DECOY's shifted neutral mass (target mass + `decoy_mz_shift_da`,
    // the same +Da the decoy grid is placed at). Only consulted when
    // `config.decoy_own_template` is set; otherwise the decoy reuses the target's
    // template below (byte-identical to the shipped paper q-values). Precomputed
    // once per feature, like `bc_templates`, so the hot loop stays allocation-free.
    let decoy_bc_templates: Vec<[f64; K_PATTERN]> = if config.decoy_own_template {
        consensus
            .iter()
            .map(|cf| isotope_model.distribution(cf.neutral_mass + config.decoy_mz_shift_da))
            .collect()
    } else {
        Vec::new()
    };
    let decoy_theoretical_patterns: Vec<Vec<f64>> = if config.decoy_own_template {
        consensus
            .iter()
            .map(|cf| {
                // Normalised averagine pattern at the decoy mass, truncated to the
                // same length as the target's `theoretical_pattern` so the grid's
                // `n_rows` trimming/weighting behaves identically to the target.
                let k = cf.theoretical_pattern.len().min(K_PATTERN);
                let full = isotope_model.distribution(cf.neutral_mass + config.decoy_mz_shift_da);
                let sum: f64 = full[..k].iter().sum();
                if sum > 0.0 {
                    full[..k].iter().map(|&x| x / sum).collect()
                } else {
                    vec![0.0; k]
                }
            })
            .collect()
    } else {
        Vec::new()
    };

    // Precompute each run's RT range once (from features — no hills needed).
    // `RunInput::rt_range` falls back to a full sweep of `run.features` when
    // `scan_times` is empty (it always is for the align pipeline), so calling it
    // inside the hot loop would cost ~10ms × n_features × n_runs.
    let t_rt = Instant::now();
    let run_rt_ranges: Vec<(f64, f64)> = runs.par_iter().map(|r| r.rt_range()).collect();
    log::info!("[timing] run_rt_ranges precompute: {:.2?}", t_rt.elapsed());

    // Quantify RUN-MAJOR: load one run's hills, quantify every consensus feature
    // against it, then drop the hills before loading the next run. Peak memory is
    // O(one run's hills) instead of O(all runs' hills loaded at once).
    let t_quant = Instant::now();
    // Per-phase nanosecond accumulators across all rayon tasks.
    let timers = PhaseTimers::new();
    let mut all_entries: Vec<LfqEntry> = Vec::new();
    let mut hills_per_run: Vec<usize> = vec![0; n_runs];

    if !config.detected_use_grid {
        log::warn!("Exclusive LFQ extraction requires grid intensities; detected_use_grid=false is ignored");
    }
    let mut priority = [Vec::new(), Vec::new()];
    // First inspect the immutable signal in every run. Then use one cross-run
    // preference per hypothesis, avoiding run-specific switching of aliases.
    // Hills remain streamed; only scalar candidate evidence survives pass one.
    for pass in 0..2 {
        if pass == 1 {
            priority = ownership::global_priority(&all_entries, &consensus);
            all_entries.clear();
        }
        for run_idx in 0..n_runs {
            let run = &runs[run_idx];
            let run_rt_range = run_rt_ranges[run_idx];
            let t_load = Instant::now();
            let hills_vec = load_hills(run_idx);
            hills_per_run[run_idx] = hills_vec.len();
            let sorted = SortedHills::from_hills(&hills_vec);
            log::info!(
                "[lfq] run {}/{} '{}': {} hills loaded [{:.2?}]",
                run_idx + 1,
                n_runs,
                run.name,
                hills_vec.len(),
                t_load.elapsed()
            );

            let extract = |feat_idx: usize, is_decoy: bool, excluded: &HashSet<SampleId>| {
                let cf = &consensus[feat_idx];
                let mut grid = XicGrid::empty(config.n_isotopes, config.grid_cols, 0.0, 1.0);
                let mut scores = ColumnScores::new(config.grid_cols);
                let mut col_totals = vec![0.0f32; config.grid_cols];
                let mut obs = vec![0.0f64; config.n_isotopes];
                let (mut corr_rt, corr_mz, mut corr_im) = if run_idx == alignment.reference_idx {
                    (cf.ref_rt, cf.ref_mz, cf.ref_im)
                } else {
                    let al = &alignment.alignments[&run.name];
                    (
                        al.predict_run_rt(cf.ref_rt),
                        al.predict_run_mz(cf.ref_mz, cf.ref_rt),
                        al.predict_run_im(cf.ref_im, cf.ref_rt),
                    )
                };
                if let Some(g) = guidance {
                    // A same-run MS2 observation supplies native RT/IM. Its paired
                    // decoy uses the same centre before the usual coordinate shifts.
                    if let Some(id) = g.native_anchor(feat_idx, run_idx) {
                        corr_rt = id.rt_minutes;
                        corr_im = id.im;
                    } else if cf.ref_im == 0.0 {
                        corr_im = 0.0;
                    }
                }
                let span = run_rt_range.1 - run_rt_range.0;
                let rt = corr_rt
                    - if is_decoy {
                        config.decoy_rt_shift_pct * span
                    } else {
                        0.0
                    };
                let mz = corr_mz
                    + if is_decoy {
                        config.decoy_mz_shift_da / cf.charge as f64
                    } else {
                        0.0
                    };
                let rt_sigma = if config.rt_spread_scoring && run_idx != alignment.reference_idx {
                    let al = &alignment.alignments[&run.name];
                    let norm = if span > 0.0 {
                        ((rt - run_rt_range.0) / span).clamp(0.0, 1.0)
                    } else {
                        0.5
                    };
                    al.rt_sigma.sigma_norm_at(norm).map(|s| {
                        (s * config.grid_cols as f64 / (2.0 * config.rt_window_pct)) as f32
                    })
                } else {
                    None
                };
                let (pattern, template) = if is_decoy && config.decoy_own_template {
                    (
                        decoy_theoretical_patterns[feat_idx].as_slice(),
                        &decoy_bc_templates[feat_idx],
                    )
                } else {
                    (cf.theoretical_pattern.as_slice(), &bc_templates[feat_idx])
                };
                let mut candidate = quantify_cell(
                    &mut grid,
                    &mut scores,
                    &mut col_totals,
                    &mut obs,
                    &hills_vec,
                    &sorted,
                    &run.scan_times,
                    pattern,
                    template,
                    config,
                    &timers,
                    feat_idx,
                    run_idx,
                    cf.charge,
                    mz,
                    rt,
                    corr_im,
                    config.rt_window_pct * span,
                    rt_sigma,
                    is_decoy,
                    !is_decoy
                        && guidance.map_or(cf.per_run_feature[run_idx].is_none(), |g| {
                            !g.is_direct(feat_idx, run_idx)
                        }),
                    excluded,
                );
                // Apply the same minimum envelope support to inferred targets
                // and their paired decoys before ownership and confidence scoring.
                if guidance.is_some_and(|g| g.is_inferred(feat_idx, run_idx))
                    && (candidate.entry.n_isotopes_found < 2 || candidate.entry.coelution < 0.5)
                {
                    candidate.entry.intensity = 0.0;
                    candidate.support.clear();
                    candidate.context.clear();
                }
                candidate
            };
            let empty = HashSet::new();
            let candidates: Vec<CellCandidate> = (0..n_features)
                .into_par_iter()
                .filter(|&i| guidance.is_none_or(|g| g.attempted(i, run_idx)))
                .flat_map_iter(|i| {
                    let mut cells = vec![extract(i, false, &empty)];
                    if config.run_tdc {
                        cells.push(extract(i, true, &empty));
                    }
                    cells
                })
                .collect();
            if pass == 0 {
                all_entries.extend(candidates.into_iter().map(|c| c.entry));
            } else {
                let resolved = ownership::resolve_run(candidates, &priority, &hills_vec, extract);
                let conflicts = resolved.iter().filter(|e| e.excluded_samples > 0).count();
                log::info!(
                    "[lfq ownership] '{}': {} cells reassessed for shared signal",
                    run.name,
                    conflicts
                );
                all_entries.extend(resolved);
            }
            // `hills_vec` and `sorted` are dropped here before the next run loads.
        }
    }
    log::info!(
        "[timing] quantification loop ({} entries, {} runs streamed): {:.2?}",
        all_entries.len(),
        n_runs,
        t_quant.elapsed()
    );
    let ns_b = timers.build.load(Ordering::Relaxed);
    let ns_s = timers.score.load(Ordering::Relaxed);
    let ns_i = timers.integrate.load(Ordering::Relaxed);
    let total_ns = ns_b + ns_s + ns_i;
    let pct = |ns: u64| {
        if total_ns == 0 {
            0.0
        } else {
            100.0 * ns as f64 / total_ns as f64
        }
    };
    log::info!(
        "[timing] hot-loop breakdown (cpu-ns across all threads): \
         build_grid={:.2}s ({:.0}%), score_grid={:.2}s ({:.0}%), integrate={:.2}s ({:.0}%)",
        ns_b as f64 / 1e9,
        pct(ns_b),
        ns_s as f64 / 1e9,
        pct(ns_s),
        ns_i as f64 / 1e9,
        pct(ns_i),
    );

    // Target-decoy q-values
    let t_tdc = Instant::now();
    let q_map: HashMap<(usize, usize), f64> = if config.run_tdc {
        if let Some(g) = guidance {
            let peptides: Vec<_> = g
                .targets
                .iter()
                .map(|t| t.modified_peptide.clone())
                .collect();
            let charges: Vec<_> = g.targets.iter().map(|t| t.charge).collect();
            // Do not let abundant direct IDs dilute the transferred-cell null.
            // Partition paired decoys using the same donor mask as targets.
            let mut q = HashMap::new();
            // Keep same-run inferred-charge evidence separate from both direct
            // PSMs and cross-run transfers, for targets and paired decoys alike.
            let class = |e: &LfqEntry| {
                if g.is_direct(e.feature_idx, e.run_idx) {
                    1
                } else if g.is_inferred(e.feature_idx, e.run_idx)
                    && g.native_anchor(e.feature_idx, e.run_idx).is_some()
                {
                    2
                } else {
                    0
                }
            };
            for evidence in 0..=2 {
                let entries: Vec<_> = all_entries
                    .iter()
                    .filter(|e| class(e) == evidence)
                    .cloned()
                    .collect();
                q.extend(rescore::search::compute(
                    &entries,
                    &peptides,
                    &charges,
                    config.search_scoring,
                ));
            }
            q
        } else {
            rescore::compute_qvalues_qda(&all_entries)
        }
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
    log_detection_summary(&all_entries, runs, n_features, n_runs, config.run_tdc);
    log::info!(
        "[timing] per-run detection summary: {:.2?}",
        t_summary.elapsed()
    );

    // Assemble output matrix
    assemble_matrix(
        &consensus,
        runs,
        alignment,
        config,
        all_entries,
        &q_map,
        hills_per_run,
        n_features,
        n_runs,
    )
}
