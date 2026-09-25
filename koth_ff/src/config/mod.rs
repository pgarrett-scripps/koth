//! Configuration. The `koth_ff` detection settings (`[file]`, `[hills]`,
//! `[hills_ms2]`, `[features]`, `[scoring]`, `[output]` and [`KothConfig`]) live in
//! `koth-core` and are re-exported here at their original paths; this module adds
//! the `koth_align` configuration.

use serde::{Deserialize, Serialize};

use crate::alignment::AlignmentConfig;
use crate::lfq::LfqConfig;

mod align;

pub use align::AlignOutputConfig;
pub use koth_core::config::*;

/// Top-level configuration for the `koth_align` binary.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct AlignConfig {
    pub alignment: AlignmentConfig,
    pub lfq: LfqConfig,
    pub output: AlignOutputConfig,
}

impl AlignConfig {
    pub fn from_toml(path: &std::path::Path) -> Result<Self, crate::error::KothError> {
        let content = std::fs::read_to_string(path)?;
        let cfg: Self = toml::from_str(&content)
            .map_err(|e| crate::error::KothError::ConfigError(e.to_string()))?;
        cfg.lfq.consensus.validate()?;
        Ok(cfg)
    }

    pub fn to_toml_string(&self) -> Result<String, crate::error::KothError> {
        toml::to_string(self).map_err(|e| crate::error::KothError::ConfigError(e.to_string()))
    }
}

#[cfg(test)]
mod config_parse_tests {
    use super::AlignConfig;

    fn parse_align(name: &str, toml_src: &str) {
        if let Err(e) = toml::from_str::<AlignConfig>(toml_src) {
            panic!("shipped koth_align config `{name}` failed to parse: {e}");
        }
    }

    // These used to also parse eight tuned configs from `benchmark/config/`. The
    // benchmark now lives in its own repository (pgarrett-scripps/koth-paper), and
    // `include_str!` resolves at COMPILE time, so leaving them would not merely
    // skip a test -- it would stop this crate building at all.
    //
    // The guard still covers everything this repository ships. The tuned configs
    // are guarded on the other side by `#[serde(deny_unknown_fields)]`: a renamed
    // or removed field makes koth_ff reject them the moment the benchmark runs,
    // which is the only thing that runs them.

    #[test]
    fn shipped_align_configs_parse() {
        parse_align(
            "example_config_align.toml",
            include_str!("../../../example_config_align.toml"),
        );
    }

    #[test]
    fn removed_consensus_gates_are_rejected_instead_of_ignored() {
        for (key, value) in [
            ("min_member_combined_score", "0.5"),
            ("min_seed_combined_score", "0.75"),
            ("min_group_size", "2"),
            ("allow_replicated_weak_seeds", "true"),
        ] {
            let source = format!("[lfq.consensus]\n{key} = {value}\n");
            let err = toml::from_str::<AlignConfig>(&source).unwrap_err();
            assert!(err.to_string().contains(key));
        }
    }

    #[test]
    fn consensus_gate_rejects_invalid_values() {
        let mut cfg = AlignConfig::default();
        for value in [f64::NAN, f64::INFINITY, -0.01, 1.01] {
            cfg.lfq.consensus.max_group_qvalue = value;
            assert!(cfg.lfq.consensus.validate().is_err());
        }
        cfg.lfq.consensus.max_group_qvalue = 0.05;
        cfg.lfq.consensus.rt_window_pct = -0.1;
        assert!(cfg.lfq.consensus.validate().is_err());
    }

    #[test]
    fn omitted_lfq_windows_use_release_defaults_without_changing_grouping() {
        for text in [
            "",
            "[lfq]\n",
            include_str!("../../../example_config_align.toml"),
        ] {
            let cfg: AlignConfig = toml::from_str(text).expect("valid alignment config");
            assert_eq!(cfg.lfq.rt_window_pct, 0.005);
            assert_eq!(cfg.lfq.im_tolerance, 0.015);
            assert_eq!(cfg.lfq.consensus.rt_window_pct, 0.02);
            assert_eq!(cfg.lfq.consensus.im_tolerance, 0.05);
            assert_eq!(cfg.alignment.im_tolerance, 0.05);
        }
    }

    #[test]
    fn explicit_lfq_windows_override_defaults_without_changing_grouping() {
        let cfg: AlignConfig =
            toml::from_str("[lfq]\nrt_window_pct = 0.01\nim_tolerance = 0.025\n")
                .expect("valid custom extraction windows");
        assert_eq!(cfg.lfq.rt_window_pct, 0.01);
        assert_eq!(cfg.lfq.im_tolerance, 0.025);
        assert_eq!(cfg.lfq.consensus.rt_window_pct, 0.02);
        assert_eq!(cfg.lfq.consensus.im_tolerance, 0.05);
    }
}
