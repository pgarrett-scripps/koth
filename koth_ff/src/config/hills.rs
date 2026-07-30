use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HillsConfig {
    pub min_scans: usize,
    pub max_gap: usize,
    pub split_hills: bool,
    /// [persistence] A split between two adjacent peaks is kept only if the
    /// valley between them drops to <= this fraction of the SMALLER peak.
    /// Scale-free and immune to spike inflation. Default 0.70.
    #[serde(default = "default_split_valley_ratio")]
    pub split_valley_ratio: f64,
    /// [persistence] Absolute notch-depth floor as a multiple of the estimated
    /// noise sigma (MAD of the trace's first difference). Rejects shallow noise
    /// notches that the relative test alone would admit on a low baseline.
    /// Default 4.0.
    #[serde(default = "default_split_sigma_mult")]
    pub split_sigma_mult: f64,
    /// [persistence] Minimum peak height as a fraction of the robust (95th-pct)
    /// max intensity. Lower catches faint minor co-eluting peaks (e.g. 10:1
    /// duals); higher rejects baseline bumps. Default 0.10.
    #[serde(default = "default_split_height_frac")]
    pub split_height_frac: f64,
    /// Weight for the intensity LFC term in the hill-candidate distance score.
    /// 0.0 disables it. ~0.5 gives intensity consistency roughly half the
    /// influence of m/z proximity when selecting which peak extends a hill.
    pub lfc_weight: f64,
    /// Linear-interpolate intensity through internal zero-gap scans during
    /// hill finalization. Independent of `smoothing_enabled`. Default false.
    /// NOTE: turning this on can create artificial local maxima that the
    /// split-hills logic picks up as new peaks — it materially increases
    /// final hill count on data with many short gaps.
    #[serde(default)]
    pub gap_fill_enabled: bool,
    /// Apply a running-average filter to the intensity profile during hill
    /// finalization. Independent of `gap_fill_enabled`. Default false.
    #[serde(default)]
    pub smoothing_enabled: bool,
    /// Half-width of the running-average window. Total window = 2*smoothing_window+1 scans.
    /// 0 = no averaging. Ignored when smoothing_enabled is false.
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
    /// Pre-detection anomaly normalisation. When > 0, peak intensities
    /// in each MS1 scan are scaled by `clamp(local_median(reference) /
    /// scan_reference, tic_norm_min_scale, tic_norm_max_scale)` where
    /// the local median is taken over a centred window of
    /// `tic_norm_window` scans. The reference quantity is selected by
    /// `tic_norm_mode`. Defaults to 0 = off.
    ///
    /// Use case: brief ESI dropouts / sub-AGC-failure events that the
    /// instrument's automatic gain control didn't fully compensate.
    /// Choosing the window too small (≲ peak FWHM in scans) will
    /// flatten real apices, especially in `tic` mode.
    #[serde(default)]
    pub tic_norm_window: usize,
    #[serde(default = "default_tic_norm_min_scale")]
    pub tic_norm_min_scale: f64,
    #[serde(default = "default_tic_norm_max_scale")]
    pub tic_norm_max_scale: f64,
    /// Which per-scan quantity to use as the normalisation reference.
    /// "tic"    = total ion current (sum of all peak intensities). Heavy-
    ///            tailed; dominated by a few intense peaks. Real elution
    ///            apices spike the TIC and get scaled DOWN by the
    ///            normaliser — apex shapes can be flattened.
    /// "median" = median peak intensity per scan. Robust to a few intense
    ///            peaks: apex elution adds tall peaks but barely moves
    ///            the median. Dropouts dim everything → median drops →
    ///            scaling fires only when the whole scan is suppressed.
    ///            Recommended.
    #[serde(default = "default_tic_norm_mode")]
    pub tic_norm_mode: String,
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

fn default_tic_norm_min_scale() -> f64 {
    0.5
}

fn default_tic_norm_max_scale() -> f64 {
    3.0
}

fn default_tic_norm_mode() -> String {
    "median".to_string()
}

fn default_split_valley_ratio() -> f64 {
    0.70
}

fn default_split_sigma_mult() -> f64 {
    4.0
}

fn default_split_height_frac() -> f64 {
    0.10
}

/// Per-field overrides applied to MS2 (DIA fragment) hill detection.
///
/// **Override semantics.** MS2 hill detection *inherits* the user's entire
/// `[hills]` config and then overwrites only the fields set here. A field left
/// `None` (absent from `[hills_ms2]` in the TOML) keeps its `[hills]` value —
/// it does **not** fall back to [`HillsConfig::default`]. When the whole
/// `[hills_ms2]` table is absent, MS2 uses `[hills]` verbatim (see
/// [`super::KothConfig::ms2_hills`]), which is byte-for-byte the pre-existing
/// behavior.
///
/// Every field mirrors the same-named field on [`HillsConfig`]; see there for
/// what each one does.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct HillsMs2Overrides {
    // keep in sync with HillsConfig — one Option<T> per overridable field, and a
    // matching `if let Some(..)` in `apply_to`. The `every_override_field_is_applied`
    // guard test fails if a field here is not wired into `apply_to`.
    pub min_scans: Option<usize>,
    pub max_gap: Option<usize>,
    pub split_hills: Option<bool>,
    pub split_valley_ratio: Option<f64>,
    pub split_sigma_mult: Option<f64>,
    pub split_height_frac: Option<f64>,
    pub lfc_weight: Option<f64>,
    pub gap_fill_enabled: Option<bool>,
    pub smoothing_enabled: Option<bool>,
    pub smoothing_window: Option<usize>,
    pub filter_large_baseline_hills: Option<bool>,
    pub large_hill_min_scans: Option<usize>,
    pub large_hill_peak_factor: Option<f64>,
    pub tic_norm_window: Option<usize>,
    pub tic_norm_min_scale: Option<f64>,
    pub tic_norm_max_scale: Option<f64>,
    pub tic_norm_mode: Option<String>,
}

impl HillsMs2Overrides {
    /// Return a copy of `base` (`[hills]`) with each `Some(..)` override applied.
    /// Fields left `None` keep their `base` value.
    pub fn apply_to(&self, base: &HillsConfig) -> HillsConfig {
        let mut c = base.clone();
        if let Some(v) = self.min_scans {
            c.min_scans = v;
        }
        if let Some(v) = self.max_gap {
            c.max_gap = v;
        }
        if let Some(v) = self.split_hills {
            c.split_hills = v;
        }
        if let Some(v) = self.split_valley_ratio {
            c.split_valley_ratio = v;
        }
        if let Some(v) = self.split_sigma_mult {
            c.split_sigma_mult = v;
        }
        if let Some(v) = self.split_height_frac {
            c.split_height_frac = v;
        }
        if let Some(v) = self.lfc_weight {
            c.lfc_weight = v;
        }
        if let Some(v) = self.gap_fill_enabled {
            c.gap_fill_enabled = v;
        }
        if let Some(v) = self.smoothing_enabled {
            c.smoothing_enabled = v;
        }
        if let Some(v) = self.smoothing_window {
            c.smoothing_window = v;
        }
        if let Some(v) = self.filter_large_baseline_hills {
            c.filter_large_baseline_hills = v;
        }
        if let Some(v) = self.large_hill_min_scans {
            c.large_hill_min_scans = v;
        }
        if let Some(v) = self.large_hill_peak_factor {
            c.large_hill_peak_factor = v;
        }
        if let Some(v) = self.tic_norm_window {
            c.tic_norm_window = v;
        }
        if let Some(v) = self.tic_norm_min_scale {
            c.tic_norm_min_scale = v;
        }
        if let Some(v) = self.tic_norm_max_scale {
            c.tic_norm_max_scale = v;
        }
        if let Some(ref v) = self.tic_norm_mode {
            c.tic_norm_mode = v.clone();
        }
        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Guard against a `HillsMs2Overrides` field that exists but was never wired
    /// into `apply_to`. Build a `[hills]` where every field is non-default, and a
    /// `[hills_ms2]` overriding every field to a value distinct from that base;
    /// after `apply_to`, every field must have changed. An un-wired field keeps
    /// its base value and fails the matching `assert_ne!`.
    #[test]
    fn every_override_field_is_applied() {
        let base = HillsConfig {
            min_scans: 5,
            max_gap: 2,
            split_hills: false,
            split_valley_ratio: 0.5,
            split_sigma_mult: 3.0,
            split_height_frac: 0.2,
            lfc_weight: 0.7,
            gap_fill_enabled: true,
            smoothing_enabled: true,
            smoothing_window: 3,
            filter_large_baseline_hills: true,
            large_hill_min_scans: 30,
            large_hill_peak_factor: 1.5,
            tic_norm_window: 7,
            tic_norm_min_scale: 0.4,
            tic_norm_max_scale: 2.5,
            tic_norm_mode: "tic".to_string(),
        };
        let ov = HillsMs2Overrides {
            min_scans: Some(9),
            max_gap: Some(4),
            split_hills: Some(true),
            split_valley_ratio: Some(0.8),
            split_sigma_mult: Some(6.0),
            split_height_frac: Some(0.3),
            lfc_weight: Some(0.9),
            gap_fill_enabled: Some(false),
            smoothing_enabled: Some(false),
            smoothing_window: Some(6),
            filter_large_baseline_hills: Some(false),
            large_hill_min_scans: Some(50),
            large_hill_peak_factor: Some(3.0),
            tic_norm_window: Some(11),
            tic_norm_min_scale: Some(0.6),
            tic_norm_max_scale: Some(4.0),
            tic_norm_mode: Some("median".to_string()),
        };
        let got = ov.apply_to(&base);

        assert_ne!(got.min_scans, base.min_scans);
        assert_ne!(got.max_gap, base.max_gap);
        assert_ne!(got.split_hills, base.split_hills);
        assert_ne!(got.split_valley_ratio, base.split_valley_ratio);
        assert_ne!(got.split_sigma_mult, base.split_sigma_mult);
        assert_ne!(got.split_height_frac, base.split_height_frac);
        assert_ne!(got.lfc_weight, base.lfc_weight);
        assert_ne!(got.gap_fill_enabled, base.gap_fill_enabled);
        assert_ne!(got.smoothing_enabled, base.smoothing_enabled);
        assert_ne!(got.smoothing_window, base.smoothing_window);
        assert_ne!(
            got.filter_large_baseline_hills,
            base.filter_large_baseline_hills
        );
        assert_ne!(got.large_hill_min_scans, base.large_hill_min_scans);
        assert_ne!(got.large_hill_peak_factor, base.large_hill_peak_factor);
        assert_ne!(got.tic_norm_window, base.tic_norm_window);
        assert_ne!(got.tic_norm_min_scale, base.tic_norm_min_scale);
        assert_ne!(got.tic_norm_max_scale, base.tic_norm_max_scale);
        assert_ne!(got.tic_norm_mode, base.tic_norm_mode);

        // Sanity: applied config equals the override values exactly.
        assert_eq!(got.min_scans, 9);
        assert_eq!(got.tic_norm_mode, "median");
    }

    /// A `None` override leaves the base value untouched (the documented
    /// inherit-not-reset semantics).
    #[test]
    fn none_override_keeps_base_value() {
        let base = HillsConfig {
            min_scans: 5,
            ..HillsConfig::default()
        };
        let got = HillsMs2Overrides::default().apply_to(&base);
        assert_eq!(got.min_scans, 5);
    }
}

impl Default for HillsConfig {
    fn default() -> Self {
        Self {
            min_scans: 3,
            max_gap: 0,
            split_hills: true,
            split_valley_ratio: default_split_valley_ratio(),
            split_sigma_mult: default_split_sigma_mult(),
            split_height_frac: default_split_height_frac(),
            lfc_weight: 0.5,
            gap_fill_enabled: false,
            smoothing_enabled: false,
            smoothing_window: 1,
            filter_large_baseline_hills: false,
            large_hill_min_scans: default_large_hill_min_scans(),
            large_hill_peak_factor: default_large_hill_peak_factor(),
            tic_norm_window: 0,
            tic_norm_min_scale: default_tic_norm_min_scale(),
            tic_norm_max_scale: default_tic_norm_max_scale(),
            tic_norm_mode: default_tic_norm_mode(),
        }
    }
}
