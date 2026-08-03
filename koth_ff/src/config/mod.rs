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
pub use hills::{HillsConfig, HillsMs2Overrides};
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
    /// Optional per-field overrides for MS2 (DIA fragment) hill detection.
    /// Absent (`None`) ⇒ MS2 hills use `[hills]` verbatim (byte-identical to the
    /// pre-existing behavior). Present ⇒ MS2 hill detection inherits `[hills]`
    /// and applies only the fields set here. See [`KothConfig::ms2_hills`].
    #[serde(default)]
    pub hills_ms2: Option<HillsMs2Overrides>,
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

    /// Resolve the [`HillsConfig`] to use for MS2 (DIA fragment) hill detection.
    ///
    /// When `[hills_ms2]` is absent this returns `self.hills` unchanged — the MS2
    /// path is byte-for-byte identical to using `[hills]` directly. When present,
    /// MS2 hill detection inherits every `[hills]` value and overwrites only the
    /// fields set under `[hills_ms2]` (see [`HillsMs2Overrides`]).
    ///
    /// Pipeline entry points (`run_pipeline_with_ms2`, `run_pipeline_streaming`,
    /// …) and the binary route MS2 detection through this method, so a single
    /// `KothConfig` drives MS1 and MS2 hill settings independently.
    pub fn ms2_hills(&self) -> HillsConfig {
        match &self.hills_ms2 {
            Some(ov) => ov.apply_to(&self.hills),
            None => self.hills.clone(),
        }
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

    // These used to also parse eight tuned configs from `benchmark/config/`. The
    // benchmark now lives in its own repository (tacular-omics/koth-paper), and
    // `include_str!` resolves at COMPILE time, so leaving them would not merely
    // skip a test -- it would stop this crate building at all.
    //
    // The guard still covers everything this repository ships. The tuned configs
    // are guarded on the other side by `#[serde(deny_unknown_fields)]`: a renamed
    // or removed field makes koth_ff reject them the moment the benchmark runs,
    // which is the only thing that runs them.

    #[test]
    fn shipped_koth_configs_parse() {
        parse_koth("example_config.toml", include_str!("../../../example_config.toml"));
    }

    #[test]
    fn shipped_align_configs_parse() {
        parse_align("example_config_align.toml", include_str!("../../../example_config_align.toml"));
    }
}

#[cfg(test)]
mod hills_ms2_tests {
    use super::KothConfig;

    /// A valid `KothConfig` whose `[hills]` carries a few *non-default* values,
    /// so tests can tell "inherited from `[hills]`" apart from "fell back to
    /// `HillsConfig::default()`" (defaults are min_scans=3, max_gap=0,
    /// split_hills=true, lfc_weight=0.5, smoothing_enabled=false,
    /// smoothing_window=1). Built from `default()` + a valid `[file]` so we
    /// don't have to hardcode every required `[file]` field.
    fn base_config_toml() -> String {
        let mut cfg = KothConfig::default();
        cfg.hills.min_scans = 5;
        cfg.hills.max_gap = 2;
        cfg.hills.split_hills = false;
        cfg.hills.lfc_weight = 0.9;
        cfg.hills.smoothing_enabled = true;
        cfg.hills.smoothing_window = 3;
        cfg.file.mz_tolerance = 8.0;
        cfg.to_toml_string().expect("serialize base config")
    }

    /// (a) Absent `[hills_ms2]` ⇒ `ms2_hills()` is byte-identical to `hills`.
    /// This is the regression guard for the default path.
    #[test]
    fn absent_hills_ms2_equals_hills() {
        let cfg: KothConfig = toml::from_str(&base_config_toml()).expect("parse");
        assert!(cfg.hills_ms2.is_none());
        let ms2 = cfg.ms2_hills();
        // Compare via TOML serialization for a full field-by-field equality that
        // does not depend on HillsConfig deriving PartialEq.
        let a = toml::to_string(&cfg.hills).unwrap();
        let b = toml::to_string(&ms2).unwrap();
        assert_eq!(a, b, "ms2_hills() must equal [hills] when [hills_ms2] absent");
    }

    /// (b) A partial `[hills_ms2]` overrides only the named fields and inherits
    /// every other field from `[hills]` (NOT from HillsConfig::default()).
    #[test]
    fn partial_hills_ms2_overrides_named_inherits_rest() {
        let toml_src =
            format!("{}\n[hills_ms2]\nmin_scans = 2\nsplit_hills = true\n", base_config_toml());
        let cfg: KothConfig = toml::from_str(&toml_src).expect("parse");
        assert!(cfg.hills_ms2.is_some());
        let ms2 = cfg.ms2_hills();

        // Overridden fields take the [hills_ms2] value.
        assert_eq!(ms2.min_scans, 2);
        assert!(ms2.split_hills);

        // Un-overridden fields inherit from [hills] (max_gap=2, lfc_weight=0.9,
        // smoothing on / window=3), which differ from the HillsConfig defaults
        // (0 / 0.5 / false / 1) — proving inheritance, not default fallback.
        assert_eq!(ms2.max_gap, 2, "inherited from [hills], not default");
        assert_eq!(ms2.lfc_weight, 0.9, "inherited from [hills], not default");
        assert!(ms2.smoothing_enabled, "inherited from [hills]");
        assert_eq!(ms2.smoothing_window, 3, "inherited from [hills]");

        // [hills] itself is untouched by the MS2 override.
        assert_eq!(cfg.hills.min_scans, 5);
        assert!(!cfg.hills.split_hills);
    }

    /// (d) `deny_unknown_fields` still rejects a typo inside `[hills_ms2]`.
    #[test]
    fn hills_ms2_rejects_unknown_field() {
        let toml_src = format!("{}\n[hills_ms2]\nmin_scanz = 2\n", base_config_toml());
        let err = toml::from_str::<KothConfig>(&toml_src);
        assert!(err.is_err(), "typo inside [hills_ms2] must be rejected");
    }
}
