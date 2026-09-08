use serde::{Deserialize, Serialize};

use super::OutputFormat;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ScoringConfig {
    /// Whether to test the ±1 neutron-offset monoisotope reassignment during
    /// scoring. `false` (default) tests only offset 0 (no reassignment);
    /// `true` searches offsets [-1, +1] inclusive.
    pub isotope_offset_enabled: bool,
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
            isotope_offset_enabled: false,
            offset_zero_bonus: 0.15,
            min_isotope_score_for_offset: 0.5,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
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
