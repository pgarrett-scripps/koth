use serde::{Deserialize, Serialize};

use super::OutputFormat;

/// Output settings for the alignment + LFQ stage.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct AlignOutputConfig {
    /// Export the versioned long-format observation and extraction bundle.
    pub export_long: bool,
    /// Requested wide-output format; the alignment writer currently falls back
    /// to TSV with a warning when Parquet is requested.
    pub format: OutputFormat,
    /// Quality threshold for n_runs_detected in consensus_features.tsv.
    /// Does not filter the exported intensity matrix.
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
            export_long: false,
            format: OutputFormat::Tsv,
            max_qvalue: 1.0,
            export_decoys: default_export_decoys(),
            export_details: default_export_details(),
        }
    }
}
