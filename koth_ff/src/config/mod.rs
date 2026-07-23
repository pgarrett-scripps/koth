use serde::{Deserialize, Serialize};

use crate::alignment::AlignmentConfig;
use crate::lfq::LfqConfig;

mod align;
mod features;
mod file;
mod hills;
mod scoring;

// Re-export every config type at the `crate::config::*` path so that existing
// import sites (`crate::config::FileConfig`, `crate::config::CosineAnchor`, …)
// keep resolving unchanged after the split into submodules.
pub use align::AlignOutputConfig;
pub use features::{CosineAnchor, FeaturesConfig};
pub use file::{FileConfig, ImToleranceType, ToleranceType};
pub use hills::HillsConfig;
pub use scoring::{OutputConfig, ScoringConfig};

/// Output format for hills / features / matrix files: "tsv" or "parquet".
/// Shared by [`OutputConfig`] (koth_ff) and [`AlignOutputConfig`] (koth_align).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum OutputFormat {
    #[default]
    Tsv,
    Parquet,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct KothConfig {
    pub file: FileConfig,
    pub hills: HillsConfig,
    pub features: FeaturesConfig,
    pub scoring: ScoringConfig,
    pub output: OutputConfig,
}

/// Top-level configuration for the `koth_align` binary.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
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
        let cfg: KothConfig = toml::from_str(&content)
            .map_err(|e| crate::error::KothError::ConfigError(e.to_string()))?;
        cfg.features.validate()?;
        Ok(cfg)
    }

    pub fn to_toml_string(&self) -> Result<String, crate::error::KothError> {
        toml::to_string(self).map_err(|e| crate::error::KothError::ConfigError(e.to_string()))
    }
}

#[cfg(test)]
mod config_parse_tests {
    use super::{AlignConfig, KothConfig};

    /// Every shipped config must parse under `deny_unknown_fields`. If a field
    /// is renamed or removed in the code, the stale key in one of these files
    /// makes this test fail — which is the point: config drift is now caught at
    /// build time rather than silently ignored at runtime.
    fn parse_koth(name: &str, toml_src: &str) {
        if let Err(e) = toml::from_str::<KothConfig>(toml_src) {
            panic!("shipped koth_ff config `{name}` failed to parse: {e}");
        }
    }

    fn parse_align(name: &str, toml_src: &str) {
        if let Err(e) = toml::from_str::<AlignConfig>(toml_src) {
            panic!("shipped koth_align config `{name}` failed to parse: {e}");
        }
    }

    #[test]
    fn shipped_koth_configs_parse() {
        parse_koth("example_config.toml", include_str!("../../../example_config.toml"));
        parse_koth("koth_ff.toml", include_str!("../../../benchmark/config/koth_ff.toml"));
        parse_koth("koth_ff_bruker.toml", include_str!("../../../benchmark/config/koth_ff_bruker.toml"));
        parse_koth("koth_ff_relaxed.toml", include_str!("../../../benchmark/config/koth_ff_relaxed.toml"));
        parse_koth("koth_ff_sulfur_on.toml", include_str!("../../../benchmark/config/koth_ff_sulfur_on.toml"));
        parse_koth("koth_ff_sulfur_off.toml", include_str!("../../../benchmark/config/koth_ff_sulfur_off.toml"));
        parse_koth("koth_ff_alphapept_like.toml", include_str!("../../../benchmark/config/koth_ff_alphapept_like.toml"));
    }

    #[test]
    fn shipped_align_configs_parse() {
        parse_align("example_config_align.toml", include_str!("../../../example_config_align.toml"));
        parse_align("koth_align.toml", include_str!("../../../benchmark/config/koth_align.toml"));
        parse_align("koth_align_bruker.toml", include_str!("../../../benchmark/config/koth_align_bruker.toml"));
    }
}
