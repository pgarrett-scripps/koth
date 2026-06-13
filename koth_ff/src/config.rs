use serde::{Deserialize, Serialize};

use crate::alignment::AlignmentConfig;
use crate::lfq::LfqConfig;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ToleranceType {
    Ppm,
    Da,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ImToleranceType {
    Relative,
    Absolute,
}

/// Per-hill mass-uncertainty model used during isotope-chain extension.
///
/// - `Off` (default): legacy behaviour. Chain extension uses the flat
///   `mz_tolerance` window regardless of per-hill confidence.
/// - `Kish`: each hill carries a `mz_se` (intensity-weighted-mean standard
///   error via Kish's effective sample size). Chain extension combines the
///   instrument tolerance with both endpoints' SEs in quadrature:
///   `tol² = mz_tolerance² + (σ_mult × se_ref)² + (σ_mult × se_cand)²`.
///   Hills with confident m/z get the legacy-tight window; hills with
///   sparse / skewed intensity get a wider window, recovering low-S/N
///   peptides that the flat tolerance would have excluded.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum MzUncertaintyMode {
    #[default]
    Off,
    Kish,
}

/// File-reading and shared tolerance settings.
///
/// Tolerances defined here are used by both hill detection and feature finding.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileConfig {
    /// m/z tolerance (shared by hills and features stages)
    pub mz_tolerance: f64,
    /// "ppm" or "da"
    pub mz_tolerance_type: ToleranceType,
    /// Ion mobility tolerance (shared by hills and features stages)
    pub im_tolerance: f64,
    /// "relative" (fraction of IM value) or "absolute" (1/K0 units)
    pub im_tolerance_type: ImToleranceType,
    /// Discard peaks below this m/z. Defaults to 0.0 (no lower bound).
    #[serde(default = "default_global_min_mz")]
    pub global_min_mz: f64,
    /// Discard peaks above this m/z. Defaults to +inf (no upper bound).
    #[serde(default = "default_global_max_mz")]
    pub global_max_mz: f64,
    /// Fraction of total hill intensity to retain (0.0–1.0). 1.0 = keep all.
    pub intensity_coverage: f64,
    /// Number of threads (None = use all available CPUs)
    pub n_threads: Option<usize>,
    /// Bruker vertical-IM filter: column half-width in TOF-index units.
    /// The filter scans each TOF column over the IM (scan) axis looking for
    /// long vertical streaks of signal; this width sets how many TOF indices
    /// to the left/right of each column are summed into the column profile.
    #[serde(default = "default_bruker_filter_mz_half_width")]
    pub bruker_filter_mz_half_width: u32,
    /// Bruker vertical-IM filter: maximum consecutive empty scans tolerated
    /// inside a kept run (morphological-close radius). Default 1.
    #[serde(default = "default_bruker_filter_max_internal_gap")]
    pub bruker_filter_max_internal_gap: usize,
    /// Bruker vertical-IM filter: minimum run span (gap-inclusive) in scans
    /// for a column feature to survive. Default 5.
    #[serde(default = "default_bruker_filter_min_feature_length")]
    pub bruker_filter_min_feature_length: usize,
    /// Bruker vertical-IM filter: per-scan summed-intensity floor for a scan
    /// to count as "occupied" in the column profile. 0 = keep every nonzero
    /// scan. Default 0.
    #[serde(default)]
    pub bruker_filter_min_window_intensity: u64,
    /// Bruker vertical-IM filter: total summed intensity over the kept span
    /// (including sub-threshold cells inside gaps) required to keep a run.
    /// 0 = no floor. Default 0.
    #[serde(default)]
    pub bruker_filter_min_feature_intensity: u64,
    /// Bruker vertical-IM filter: how many times to re-apply the filter to
    /// its own survivors. Each pass is strictly more aggressive. Default 1.
    #[serde(default = "default_bruker_filter_num_iterations")]
    pub bruker_filter_num_iterations: usize,
    /// Bruker watershed centroider: nearest-neighbour reach on the scan axis
    /// (scans). Two points farther apart than this on either axis cannot
    /// join the same group. Default 10.
    #[serde(default = "default_bruker_watershed_box_scan")]
    pub bruker_watershed_box_scan: u32,
    /// Bruker watershed centroider: nearest-neighbour reach on the TOF axis
    /// (TOF indices). Default 3.
    #[serde(default = "default_bruker_watershed_box_mz_idx")]
    pub bruker_watershed_box_mz_idx: u32,
    /// Bruker watershed centroider: intensity floor for promoting an orphan
    /// (point with no in-box neighbour) to a new seed. Below this it is
    /// dropped without claiming territory. Default 0.
    #[serde(default)]
    pub bruker_watershed_min_seed_intensity: u64,
    /// Bruker watershed centroider: drop final centroids whose summed group
    /// intensity is below this floor. Default 0.
    #[serde(default)]
    pub bruker_watershed_min_centroid_total: u64,
    /// Bruker watershed centroider: hard cap on how far (in TOF-index units)
    /// any member of a group can sit from that group's seed. Prevents a
    /// long follower chain from creeping past the real peak edge. Default 10.
    #[serde(default = "default_bruker_watershed_max_tof_offset")]
    pub bruker_watershed_max_tof_offset: u32,
    /// Per-frame iterative MAD noise filter applied to the centroided peaks
    /// before the Bruker reader emits a `Spectrum`. Same algorithm as
    /// `noise_filter_sigma` (median + sigma * 1.4826 * MAD of the noise
    /// floor, found by iterative sigma-clipping). `None` (default) disables.
    /// Independent of `noise_filter_sigma`, which runs later on every spectrum
    /// (from any input format) during hill detection. Typical value: 3.0.
    #[serde(default)]
    pub bruker_noise_sigma: Option<f64>,
    /// Per-scan iterative sigma-clipping noise filter. When set, peaks whose
    /// intensity falls below `median + sigma * (1.4826 * MAD)` of the estimated
    /// noise floor are discarded before hill detection.
    /// `None` (default) disables the filter; a typical starting value is 3.0.
    pub noise_filter_sigma: Option<f64>,
    /// Decoy mode: shuffle MS1 spectra into a random order before hill and
    /// feature finding. Destroys the chromatographic structure while preserving
    /// the per-scan peak distributions, producing a null (decoy) feature set.
    #[serde(default)]
    pub decoy_mode: bool,
    /// If true, additionally detect MS2 hills (one set of hills per precursor
    /// isolation window) from mzML inputs. Written to `hills_ms2.{ext}`.
    /// Bruker .d MS2 frames are not yet supported.
    #[serde(default)]
    pub ms2_hills_enabled: bool,
    /// Enable two-pass empirical m/z tolerance calibration. Pass 1 runs hill
    /// detection at a widened tolerance (`mz_tolerance × pass1_multiplier`),
    /// records the absolute ppm-delta of every accepted peak-to-hill match,
    /// then sets the pass-2 tolerance to `median + sigma_mult × σ`, capped
    /// by the pass-1 ceiling. Only honoured when `mz_tolerance_type = Ppm`;
    /// silently skipped for Dalton tolerances.
    #[serde(default)]
    pub adaptive_mz_tolerance: bool,
    /// Multiplier applied to `mz_tolerance` for the pass-1 calibration
    /// sweep. The calibrated pass-2 tolerance is hard-capped at this same
    /// `mz_tolerance × pass1_multiplier` value, so it's also the upper
    /// safety bound. Default 2.0.
    #[serde(default = "default_adaptive_mz_tolerance_pass1_multiplier")]
    pub adaptive_mz_tolerance_pass1_multiplier: f64,
    /// `N` in `median + N × σ` when deriving the calibrated pass-2 ppm.
    /// Default 3.0 (matches AlphaPept).
    #[serde(default = "default_adaptive_mz_tolerance_sigma_mult")]
    pub adaptive_mz_tolerance_sigma_mult: f64,
    /// Per-hill mass-uncertainty model for isotope-chain extension. See
    /// `MzUncertaintyMode` for details. Default `Off` (legacy behaviour).
    #[serde(default)]
    pub mz_uncertainty_mode: MzUncertaintyMode,
    /// `σ_mult` in the combined-tolerance formula
    /// `tol² = mz_tolerance² + (σ_mult × se_ref)² + (σ_mult × se_cand)²`.
    /// Default 3.0 — wraps each hill's standard error in a 3σ envelope.
    /// Ignored when `mz_uncertainty_mode = Off`.
    #[serde(default = "default_mz_uncertainty_sigma_mult")]
    pub mz_uncertainty_sigma_mult: f64,
}

impl Default for FileConfig {
    fn default() -> Self {
        Self {
            mz_tolerance: 8.0,
            mz_tolerance_type: ToleranceType::Ppm,
            im_tolerance: 0.05,
            im_tolerance_type: ImToleranceType::Relative,
            global_min_mz: 0.0,
            global_max_mz: f64::INFINITY,
            intensity_coverage: 1.0,
            n_threads: None,
            bruker_filter_mz_half_width: default_bruker_filter_mz_half_width(),
            bruker_filter_max_internal_gap: default_bruker_filter_max_internal_gap(),
            bruker_filter_min_feature_length: default_bruker_filter_min_feature_length(),
            bruker_filter_min_window_intensity: 0,
            bruker_filter_min_feature_intensity: 0,
            bruker_filter_num_iterations: default_bruker_filter_num_iterations(),
            bruker_watershed_box_scan: default_bruker_watershed_box_scan(),
            bruker_watershed_box_mz_idx: default_bruker_watershed_box_mz_idx(),
            bruker_watershed_min_seed_intensity: 0,
            bruker_watershed_min_centroid_total: 0,
            bruker_watershed_max_tof_offset: default_bruker_watershed_max_tof_offset(),
            bruker_noise_sigma: None,
            noise_filter_sigma: None,
            decoy_mode: false,
            ms2_hills_enabled: false,
            adaptive_mz_tolerance: false,
            adaptive_mz_tolerance_pass1_multiplier: default_adaptive_mz_tolerance_pass1_multiplier(),
            adaptive_mz_tolerance_sigma_mult: default_adaptive_mz_tolerance_sigma_mult(),
            mz_uncertainty_mode: MzUncertaintyMode::Off,
            mz_uncertainty_sigma_mult: default_mz_uncertainty_sigma_mult(),
        }
    }
}

fn default_adaptive_mz_tolerance_pass1_multiplier() -> f64 {
    2.0
}

fn default_adaptive_mz_tolerance_sigma_mult() -> f64 {
    3.0
}

fn default_mz_uncertainty_sigma_mult() -> f64 {
    3.0
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HillsConfig {
    pub min_scans: usize,
    pub max_gap: usize,
    pub split_hills: bool,
    pub min_peak_distance: usize,
    pub min_peak_height: f64,
    /// Minimum prominence for a peak to trigger a split, as a fraction of the
    /// hill's maximum intensity (0.0–1.0).  Prominence = peak height minus the
    /// highest valley between the peak and any taller neighbour.  A noise wiggle
    /// sitting on the flank of a larger peak has near-zero prominence even if it
    /// passes the local-maxima check; a genuine co-eluting compound has high
    /// prominence because the valley between the two peaks is deep.
    #[serde(default = "default_min_prominence")]
    pub min_prominence: f64,
    /// Weight for the intensity LFC term in the hill-candidate distance score.
    /// 0.0 disables it. ~0.5 gives intensity consistency roughly half the
    /// influence of m/z proximity when selecting which peak extends a hill.
    pub lfc_weight: f64,
    /// Apply intensity profile smoothing after hill finalization.
    #[serde(default)]
    pub smoothing_enabled: bool,
    /// Half-width of the running-average window. Total window = 2*smoothing_window+1 scans.
    /// 0 = no averaging (gap-fill only). Ignored when smoothing_enabled is false.
    #[serde(default = "default_smoothing_window")]
    pub smoothing_window: usize,
    /// Drop "large baseline-like" hills before isotope chain assembly.
    /// A hill is dropped when its scan length ≥ `large_hill_min_scans` AND
    /// `max(unsmoothed) / smoothed_endpoint < large_hill_peak_factor` on
    /// either end — i.e., the intensity profile lacks a clear apex.
    /// Catches column-bleed contaminants, plasticizers, and baseline
    /// drift centroids that masquerade as hills. Default `false` for
    /// backward compatibility. Ported from AlphaPept's `filter_hills`.
    #[serde(default)]
    pub filter_large_baseline_hills: bool,
    /// Minimum scan span for a hill to be considered for the
    /// baseline-hill filter. Default 40.
    #[serde(default = "default_large_hill_min_scans")]
    pub large_hill_min_scans: usize,
    /// Required ratio of `max(intensity_profile)` to the smoothed
    /// endpoint intensity. A clear chromatographic peak easily exceeds
    /// this; a flat baseline trace fails. Default 2.0.
    #[serde(default = "default_large_hill_peak_factor")]
    pub large_hill_peak_factor: f64,
}

fn default_smoothing_window() -> usize {
    1
}

fn default_large_hill_min_scans() -> usize {
    40
}

fn default_large_hill_peak_factor() -> f64 {
    2.0
}

fn default_global_min_mz() -> f64 {
    0.0
}

fn default_global_max_mz() -> f64 {
    f64::INFINITY
}

fn default_bruker_filter_mz_half_width() -> u32 {
    2
}

fn default_bruker_filter_max_internal_gap() -> usize {
    1
}

fn default_bruker_filter_min_feature_length() -> usize {
    5
}

fn default_bruker_filter_num_iterations() -> usize {
    1
}

fn default_bruker_watershed_box_scan() -> u32 {
    10
}

fn default_bruker_watershed_box_mz_idx() -> u32 {
    3
}

fn default_bruker_watershed_max_tof_offset() -> u32 {
    10
}

fn default_min_prominence() -> f64 {
    0.2
}

impl Default for HillsConfig {
    fn default() -> Self {
        Self {
            min_scans: 3,
            max_gap: 0,
            split_hills: true,
            min_peak_distance: 10,
            min_peak_height: 0.2,
            min_prominence: 0.2,
            lfc_weight: 0.5,
            smoothing_enabled: false,
            smoothing_window: 1,
            filter_large_baseline_hills: false,
            large_hill_min_scans: default_large_hill_min_scans(),
            large_hill_peak_factor: default_large_hill_peak_factor(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeaturesConfig {
    pub min_charge: u8,
    pub max_charge: u8,
    /// Per-extension chromatographic-cosine threshold used while building an
    /// isotope chain. The candidate hill's cosine vs the seed hill must clear
    /// this to be added; otherwise chain extension stops in that direction.
    pub min_chain_cosine: f64,
    pub left_max_decrease: f64,
    pub right_max_decrease: f64,
    pub max_isotopes: usize,
    /// Neutron (C13) mass in Da
    pub neutron_mass: f64,
    /// Drop features whose **isotope_score** (Bhattacharyya vs averagine) is
    /// below this. 0.0 = keep all.
    pub min_isotope_score: f64,
    /// Drop features whose **cosine_score** (mean chromatographic cosine of
    /// adjacent isotope hills) is below this. 0.0 = keep all.
    pub min_cosine_score: f64,
    /// Drop features whose **combined_score** (= isotope × cosine) is below
    /// this. 0.0 = keep all. All three filters are AND-ed.
    pub min_combined_score: f64,
}

impl Default for FeaturesConfig {
    fn default() -> Self {
        Self {
            min_charge: 1,
            max_charge: 7,
            min_chain_cosine: 0.5,
            left_max_decrease: 0.05,
            right_max_decrease: 0.05,
            max_isotopes: 6,
            neutron_mass: 1.003_354_835,
            min_isotope_score: 0.0,
            min_cosine_score: 0.0,
            min_combined_score: 0.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoringConfig {
    /// Range of neutron offsets to test: [min, max] inclusive
    pub isotope_offset_min: i8,
    pub isotope_offset_max: i8,
    /// Bonus added to the Bhattacharyya score when offset == 0, to prefer
    /// the no-reassignment hypothesis when scores are close. Internal only —
    /// never stored on the feature.
    pub offset_zero_bonus: f64,
    /// Features whose best **isotope_score** (Bhattacharyya) is below this
    /// keep `neutron_offset = 0` (no monoisotopic reassignment). Does not
    /// affect retention — use `FeaturesConfig.min_isotope_score` /
    /// `min_cosine_score` / `min_combined_score` for that.
    pub min_isotope_score_for_offset: f64,
}

impl Default for ScoringConfig {
    fn default() -> Self {
        Self {
            isotope_offset_min: -1,
            isotope_offset_max: 1,
            offset_zero_bonus: 0.15,
            min_isotope_score_for_offset: 0.5,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum OutputFormat {
    #[default]
    Tsv,
    Parquet,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputConfig {
    /// Output format for hills and features files: "tsv" or "parquet"
    pub format: OutputFormat,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            format: OutputFormat::Tsv,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct KothConfig {
    pub file: FileConfig,
    pub hills: HillsConfig,
    pub features: FeaturesConfig,
    pub scoring: ScoringConfig,
    pub output: OutputConfig,
}

/// Output settings for the alignment + LFQ stage.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignOutputConfig {
    /// "tsv" or "parquet"
    pub format: OutputFormat,
    /// Only write matrix entries with q-value ≤ this threshold (1.0 = keep all)
    pub max_qvalue: f64,
    /// Write `decoy_intensity_matrix.tsv` (target-decoy decoy LFQ values).
    /// Has no effect when `lfq.run_tdc = false`.
    #[serde(default = "default_export_decoys")]
    pub export_decoys: bool,
    /// Write `lfq_details.tsv` — long-format file with one row per
    /// (feature, run, is_decoy) including scores, RT diff, observed mz/IM.
    #[serde(default = "default_export_details")]
    pub export_details: bool,
}

fn default_export_decoys() -> bool {
    true
}

fn default_export_details() -> bool {
    true
}

impl Default for AlignOutputConfig {
    fn default() -> Self {
        Self {
            format: OutputFormat::Tsv,
            max_qvalue: 1.0,
            export_decoys: default_export_decoys(),
            export_details: default_export_details(),
        }
    }
}

/// Top-level configuration for the `koth_align` binary.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AlignConfig {
    pub alignment: AlignmentConfig,
    pub lfq: LfqConfig,
    pub output: AlignOutputConfig,
}

impl AlignConfig {
    pub fn from_toml(path: &std::path::Path) -> Result<Self, crate::error::KothError> {
        let content = std::fs::read_to_string(path)?;
        toml::from_str(&content).map_err(|e| crate::error::KothError::ConfigError(e.to_string()))
    }

    pub fn to_toml_string(&self) -> Result<String, crate::error::KothError> {
        toml::to_string(self).map_err(|e| crate::error::KothError::ConfigError(e.to_string()))
    }
}

impl KothConfig {
    pub fn from_toml(path: &std::path::Path) -> Result<Self, crate::error::KothError> {
        let content = std::fs::read_to_string(path)?;
        toml::from_str(&content).map_err(|e| crate::error::KothError::ConfigError(e.to_string()))
    }

    pub fn to_toml_string(&self) -> Result<String, crate::error::KothError> {
        toml::to_string(self).map_err(|e| crate::error::KothError::ConfigError(e.to_string()))
    }
}
