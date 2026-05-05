use serde::{Deserialize, Serialize};

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
    /// Discard peaks below this m/z
    pub global_min_mz: f64,
    /// Discard peaks above this m/z
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
    /// Per-scan iterative sigma-clipping noise filter. When set, peaks whose
    /// intensity falls below `median + sigma * (1.4826 * MAD)` of the estimated
    /// noise floor are discarded before hill detection.
    /// `None` (default) disables the filter; a typical starting value is 3.0.
    pub noise_filter_sigma: Option<f64>,
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
            noise_filter_sigma: None,
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
    pub min_valley_ratio: f64,
    /// Weight for the intensity LFC term in the hill-candidate distance score.
    /// 0.0 disables it. ~0.5 gives intensity consistency roughly half the
    /// influence of m/z proximity when selecting which peak extends a hill.
    pub lfc_weight: f64,
}

impl Default for HillsConfig {
    fn default() -> Self {
        Self {
            min_scans: 3,
            max_gap: 0,
            split_hills: true,
            min_peak_distance: 10,
            min_peak_height: 0.2,
            min_valley_ratio: 0.6,
            lfc_weight: 0.5,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeaturesConfig {
    pub min_charge: u8,
    pub max_charge: u8,
    pub min_cosine_similarity: f64,
    pub left_max_decrease: f64,
    pub right_max_decrease: f64,
    pub max_isotopes: usize,
    /// Neutron (C13) mass in Da
    pub neutron_mass: f64,
}

impl Default for FeaturesConfig {
    fn default() -> Self {
        Self {
            min_charge: 1,
            max_charge: 7,
            min_cosine_similarity: 0.5,
            left_max_decrease: 0.05,
            right_max_decrease: 0.05,
            max_isotopes: 6,
            neutron_mass: 1.003_354_835,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoringConfig {
    /// Range of neutron offsets to test: [min, max] inclusive
    pub isotope_offset_min: i8,
    pub isotope_offset_max: i8,
    /// Bonus score for offset == 0 (no reassignment needed)
    pub offset_zero_bonus: f64,
    /// Features scoring below this threshold keep offset=0
    pub min_score_threshold: f64,
}

impl Default for ScoringConfig {
    fn default() -> Self {
        Self {
            isotope_offset_min: -1,
            isotope_offset_max: 1,
            offset_zero_bonus: 0.15,
            min_score_threshold: 0.5,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct KothConfig {
    pub file: FileConfig,
    pub hills: HillsConfig,
    pub features: FeaturesConfig,
    pub scoring: ScoringConfig,
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
