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
