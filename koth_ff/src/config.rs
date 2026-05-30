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
    /// For Bruker .d centroiding: m/z tolerance in ppm
    pub bruker_mz_ppm: f64,
    /// For Bruker .d centroiding: ion mobility tolerance in percent
    pub bruker_im_pct: f64,
    /// For Bruker .d centroiding: minimum raw subpeaks required per centroided peak
    pub bruker_min_subpeaks: usize,
    /// Half-window for Bruker raw-cloud smoothing along the scan (IM) axis,
    /// in scan-index units. Full window = 2*N + 1. 0 disables. Default 2 (window 5).
    #[serde(default = "default_bruker_im_smoothing")]
    pub bruker_im_smoothing_window: usize,
    /// Half-window for Bruker raw-cloud smoothing along the TOF (m/z) axis,
    /// in TOF-index units. Full window = 2*N + 1. 0 disables. Default 1 (window 3).
    #[serde(default = "default_bruker_mz_smoothing")]
    pub bruker_mz_smoothing_window: usize,
    /// Half-width (Da) of the satellite-suppression window applied after each
    /// centroid is emitted. Within this window, raw peaks falling under a
    /// linear ramp from the anchor's raw intensity (at d=0) down to
    /// anchor * `bruker_satellite_end_fraction` (at d=window) are marked used
    /// and prevented from seeding their own centroids. 0.0 disables. Default 0.15.
    #[serde(default = "default_bruker_satellite_window")]
    pub bruker_satellite_window_da: f64,
    /// End fraction for the satellite-suppression linear ramp (see
    /// `bruker_satellite_window_da`). 0.0 = full triangle ramp; 0.3 = ramp
    /// floor sits at 30% of the anchor intensity at the window edge.
    #[serde(default = "default_bruker_satellite_end_fraction")]
    pub bruker_satellite_end_fraction: f64,
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
            bruker_mz_ppm: 5.0,
            bruker_im_pct: 3.0,
            bruker_min_subpeaks: 1,
            bruker_im_smoothing_window: default_bruker_im_smoothing(),
            bruker_mz_smoothing_window: default_bruker_mz_smoothing(),
            bruker_satellite_window_da: default_bruker_satellite_window(),
            bruker_satellite_end_fraction: default_bruker_satellite_end_fraction(),
            bruker_noise_sigma: None,
            noise_filter_sigma: None,
            decoy_mode: false,
            ms2_hills_enabled: false,
        }
    }
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
}

fn default_smoothing_window() -> usize {
    1
}

fn default_global_min_mz() -> f64 {
    0.0
}

fn default_global_max_mz() -> f64 {
    f64::INFINITY
}

fn default_bruker_im_smoothing() -> usize {
    2
}

fn default_bruker_mz_smoothing() -> usize {
    1
}

fn default_bruker_satellite_window() -> f64 {
    0.15
}

fn default_bruker_satellite_end_fraction() -> f64 {
    0.3
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
