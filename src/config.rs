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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HillsConfig {
    pub mz_tolerance: f64,
    pub mz_tolerance_type: ToleranceType,
    pub min_scans: usize,
    pub max_gap: usize,
    pub split_hills: bool,
    pub min_peak_distance: usize,
    pub min_peak_height: f64,
    pub min_valley_ratio: f64,
    pub im_tolerance: f64,
    pub im_tolerance_type: ImToleranceType,
    pub n_threads: Option<usize>,
    pub global_min_mz: f64,
    pub global_max_mz: f64,
    /// For Bruker .d centroiding: m/z tolerance in ppm
    pub bruker_mz_ppm: f64,
    /// For Bruker .d centroiding: ion mobility tolerance in percent
    pub bruker_im_pct: f64,
    /// Fraction of total hill intensity to retain (0.0–1.0).
    /// Scans are kept by expanding outward from the apex, greedily choosing
    /// whichever neighbour adds more intensity, until the target fraction is
    /// reached. Trailing low-intensity tails are discarded. 1.0 = keep all.
    pub intensity_coverage: f64,
}

impl Default for HillsConfig {
    fn default() -> Self {
        Self {
            mz_tolerance: 8.0,
            mz_tolerance_type: ToleranceType::Ppm,
            min_scans: 3,
            max_gap: 0,
            split_hills: true,
            min_peak_distance: 10,
            min_peak_height: 0.2,
            min_valley_ratio: 0.6,
            im_tolerance: 0.05,
            im_tolerance_type: ImToleranceType::Relative,
            n_threads: None,
            global_min_mz: 0.0,
            global_max_mz: f64::INFINITY,
            bruker_mz_ppm: 5.0,
            bruker_im_pct: 3.0,
            intensity_coverage: 1.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeaturesConfig {
    pub mz_tolerance: f64,
    pub mz_tolerance_type: ToleranceType,
    pub min_charge: u8,
    pub max_charge: u8,
    pub min_cosine_similarity: f64,
    pub left_max_decrease: f64,
    pub right_max_decrease: f64,
    pub im_tolerance: f64,
    pub im_tolerance_type: ImToleranceType,
    pub max_isotopes: usize,
    /// Neutron (C13) mass in Da
    pub neutron_mass: f64,
}

impl Default for FeaturesConfig {
    fn default() -> Self {
        Self {
            mz_tolerance: 5.0,
            mz_tolerance_type: ToleranceType::Ppm,
            min_charge: 1,
            max_charge: 7,
            min_cosine_similarity: 0.5,
            left_max_decrease: 0.9,
            right_max_decrease: 0.9,
            im_tolerance: 0.05,
            im_tolerance_type: ImToleranceType::Absolute,
            max_isotopes: 6,
            neutron_mass: 1.003_354_835, // C13 neutron mass
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
