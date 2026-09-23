use serde::{Deserialize, Serialize};

use crate::models::Polarity;

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
#[serde(default, deny_unknown_fields)]
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
    /// Bruker vertical-IM filter: column half-width in TOF-index units.
    /// The filter scans each TOF column over the IM (scan) axis looking for
    /// long vertical streaks of signal; this width sets how many TOF indices
    /// to the left/right of each column are summed into the column profile.
    #[serde(default = "default_bruker_filter_mz_half_width")]
    pub bruker_filter_mz_half_width: u32,
    /// Bruker vertical-IM filter: maximum consecutive empty scans tolerated
    /// inside a kept run (morphological-close radius). Default 1.
    #[serde(default = "default_bruker_filter_max_internal_gap")]
    pub bruker_filter_max_internal_gap: usize,
    /// Ion polarity of the acquisition, which sets how a neutral mass is
    /// recovered from an observed m/z: `"positive"` (the default, `M + zH`) or
    /// `"negative"` (`M − zH`).
    ///
    /// Peptides are acquired in positive mode. Nucleic acids are acquired in
    /// negative mode, and getting this wrong shifts every reported `massCalib`
    /// by `2·z·1.00728` Da, which is 8 Da on a 4-charged oligonucleotide — large
    /// enough to defeat any downstream identification. koth does not read the
    /// polarity out of the file; set it alongside `[features] isotope_model`.
    #[serde(default)]
    pub polarity: Polarity,
    /// Bruker vertical-IM filter: minimum run span (gap-inclusive) in scans
    /// for a column feature to survive. Default 5.
    #[serde(default = "default_bruker_filter_min_feature_length")]
    pub bruker_filter_min_feature_length: usize,
    /// Bruker vertical-IM filter: per-scan summed-intensity floor for a scan
    /// to count as "occupied" in the column profile. 0 = keep every nonzero
    /// scan. Default 0.
    #[serde(default)]
    pub bruker_filter_min_window_intensity: u64,
    /// Bruker vertical-IM filter: total summed intensity over the kept span
    /// (including sub-threshold cells inside gaps) required to keep a run.
    /// 0 = no floor. Default 0.
    #[serde(default)]
    pub bruker_filter_min_feature_intensity: u64,
    /// Bruker vertical-IM filter: how many times to re-apply the filter to
    /// its own survivors. Each pass is strictly more aggressive. Default 1.
    #[serde(default = "default_bruker_filter_num_iterations")]
    pub bruker_filter_num_iterations: usize,
    /// Bruker watershed centroider: nearest-neighbour reach on the scan axis
    /// (scans). Two points farther apart than this on either axis cannot
    /// join the same group. Default 10.
    #[serde(default = "default_bruker_watershed_box_scan")]
    pub bruker_watershed_box_scan: u32,
    /// Bruker watershed centroider: nearest-neighbour reach on the TOF axis
    /// (TOF indices). Default 3.
    #[serde(default = "default_bruker_watershed_box_mz_idx")]
    pub bruker_watershed_box_mz_idx: u32,
    /// Bruker watershed centroider: intensity floor for promoting an orphan
    /// (point with no in-box neighbour) to a new seed. Below this it is
    /// dropped without claiming territory. Default 0.
    #[serde(default)]
    pub bruker_watershed_min_seed_intensity: u64,
    /// Bruker watershed centroider: drop final centroids whose summed group
    /// intensity is below this floor. Default 0.
    #[serde(default)]
    pub bruker_watershed_min_centroid_total: u64,
    /// Bruker watershed centroider: hard cap on how far (in TOF-index units)
    /// any member of a group can sit from that group's seed. Prevents a
    /// long follower chain from creeping past the real peak edge. Default 10.
    #[serde(default = "default_bruker_watershed_max_tof_offset")]
    pub bruker_watershed_max_tof_offset: u32,
    /// Per-frame iterative MAD noise filter applied to the centroided peaks
    /// before the Bruker reader emits a `Spectrum`. Same algorithm as
    /// `noise_filter_sigma` (median + sigma * 1.4826 * MAD of the noise
    /// floor, found by iterative sigma-clipping). `None` (default) disables.
    /// Independent of `noise_filter_sigma`, which runs later on every spectrum
    /// (from any input format) during hill detection. Typical value: 3.0.
    #[serde(default)]
    pub bruker_noise_sigma: Option<f64>,
    /// 1/K0 scale for Bruker `.d` input. `"calibrated"` (default since 0.10.0)
    /// converts each scan number with the run's acquisition calibration
    /// (`TimsCalibration` in `analysis.tdf`), matching Bruker's timsdata SDK,
    /// DataAnalysis and SDK-based tools. `"linear"` restores the straight line
    /// between the acquisition-range bounds that koth used up to 0.9.0 (it can
    /// be off by ~0.03 1/K0, about 2%). The conversion happens per centroid,
    /// before hills are built. No effect on mzML or Thermo input. The streaming
    /// path's dnoise MS1 polygon gate (`bruker_ms1_polygon`) uses the same scale.
    #[serde(default)]
    pub bruker_mobility_scale: crate::io::tims_calibration::MobilityScale,
    /// Use the in-process `dnoise` streaming API instead of the local two-stage
    /// reader. When `true`, each raw frame is run through dnoise's configured
    /// stages (vertical-IM filter -> optional horizontal halo -> optional MS1
    /// polygon -> watershed) in a single pass with no denoised `.d` written to
    /// disk. When `false` (default), the historical local path runs (vertical
    /// filter + watershed only, no halo/polygon) — this is the paper's validated
    /// pipeline. This flag selects preprocessing, not spectrum buffering:
    /// both modes stream MS1 spectra into hill detection.
    #[serde(default)]
    pub bruker_streaming: bool,
    /// Streaming path only: apply dnoise's horizontal-halo filter after the
    /// vertical filter (removes the weak m/z halo flanking bright ions). Default
    /// `true` (matches the standalone `dnoise` MS1 pipeline). No effect unless
    /// `bruker_streaming` is set.
    #[serde(default = "default_bruker_halo")]
    pub bruker_halo: bool,
    /// Streaming halo: drop a peak below this fraction of the off-column box-max
    /// reference. Default 0.15.
    #[serde(default = "default_bruker_halo_peak_fraction")]
    pub bruker_halo_peak_fraction: f64,
    /// Streaming halo: reference-box half-width along the TOF index. Default 80.
    #[serde(default = "default_bruker_halo_mz_idx_half_width")]
    pub bruker_halo_mz_idx_half_width: u32,
    /// Streaming halo: reference-box half-width along the ion-mobility scan axis.
    /// Default 2.
    #[serde(default = "default_bruker_halo_scan_half_width")]
    pub bruker_halo_scan_half_width: usize,
    /// Streaming path only: apply dnoise's ddaPASEF MS1 selection-polygon gate,
    /// dropping MS1 points outside the run's IMS PolygonFilter — the (m/z, 1/K0)
    /// region the acquisition method restricts precursor selection to. Signal
    /// outside it was never a fragmentation candidate (background, the
    /// singly-charged hump, out-of-range ions).
    ///
    /// Auto-detected: a no-op on diaPASEF and on any run that stores no polygon,
    /// so enabling it cannot fail on data that has none. Default `false`, which
    /// matches dnoise's own default and leaves every existing config's output
    /// byte-identical. No effect unless `bruker_streaming` is set.
    #[serde(default)]
    pub bruker_ms1_polygon: bool,
    /// Polygon gate: m/z leniency added to each side of the polygon interior, in
    /// Daltons. Isotopes run to higher m/z, so a pad keeps an edge precursor's
    /// envelope intact. Default 0.0 (the literal polygon).
    #[serde(default)]
    pub bruker_ms1_polygon_mz_pad: f64,
    /// Polygon gate: ion-mobility leniency added to each side, in 1/K0.
    /// Default 0.0 (the literal polygon).
    #[serde(default)]
    pub bruker_ms1_polygon_im_pad: f64,
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
    /// isolation window) from DIA mzML or diaPASEF inputs. Written to
    /// `hills_ms2.{ext}`. DDA inputs produce no MS2 hills.
    #[serde(default)]
    pub ms2_hills_enabled: bool,
    /// Enable ID-free isotope-consistency m/z recalibration. A pass-1 feature
    /// detection collects the signed ppm deviation of every adjacent
    /// isotope-hill spacing from its theoretical `neutron_mass / z` step and
    /// bins those residuals over (m/z, RT). Pass 2 shifts the *expected*
    /// isotope position during chain extension by the learned per-region
    /// median offset, so isotope hills are searched for at their recalibrated
    /// location. Skipped in decoy mode (shuffled spectra carry no real signal)
    /// and when too few isotope-spacing samples are collected. Default off.
    ///
    /// When enabled, the fixed isotope-match tolerance is additionally replaced
    /// with a region-adaptive one derived from the recalibration surface's
    /// per-region residual spread σ:
    /// `tol_ppm = clamp(tol_sigma_mult × σ(m/z,RT), tol_floor_ppm, mz_tolerance)`.
    /// Tightens the search window where the instrument is precise (rejecting
    /// false isotope matches) and relaxes it — never beyond `mz_tolerance` —
    /// where it is noisy. Only honoured for ppm tolerances. Default true.
    #[serde(default = "default_mz_recalibration")]
    pub mz_recalibration: bool,
    /// Number of m/z bins in the recalibration surface. Bin extents are
    /// derived from the observed sample range. Default 20.
    #[serde(default = "default_mz_recalibration_mz_bins")]
    pub mz_recalibration_mz_bins: usize,
    /// Number of RT bins in the recalibration surface. Default 8.
    #[serde(default = "default_mz_recalibration_rt_bins")]
    pub mz_recalibration_rt_bins: usize,
    /// Minimum residual samples a surface cell must hold before its own median
    /// is trusted; below this the cell falls back to the m/z-marginal median,
    /// then the global median. Default 50.
    #[serde(default = "default_mz_recalibration_min_samples")]
    pub mz_recalibration_min_samples: usize,
    /// `N` in the region-adaptive tolerance `N × σ`. Default 4.0.
    #[serde(default = "default_mz_recalibration_tol_sigma_mult")]
    pub mz_recalibration_tol_sigma_mult: f64,
    /// Lower bound (ppm) on the region-adaptive tolerance, so a spuriously
    /// tiny σ in a sparse region can't collapse the window. Default 1.0.
    #[serde(default = "default_mz_recalibration_tol_floor_ppm")]
    pub mz_recalibration_tol_floor_ppm: f64,
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
            bruker_filter_mz_half_width: default_bruker_filter_mz_half_width(),
            bruker_filter_max_internal_gap: default_bruker_filter_max_internal_gap(),
            polarity: Polarity::default(),
            bruker_filter_min_feature_length: default_bruker_filter_min_feature_length(),
            bruker_filter_min_window_intensity: 0,
            bruker_filter_min_feature_intensity: 0,
            bruker_filter_num_iterations: default_bruker_filter_num_iterations(),
            bruker_watershed_box_scan: default_bruker_watershed_box_scan(),
            bruker_watershed_box_mz_idx: default_bruker_watershed_box_mz_idx(),
            bruker_watershed_min_seed_intensity: 0,
            bruker_watershed_min_centroid_total: 0,
            bruker_watershed_max_tof_offset: default_bruker_watershed_max_tof_offset(),
            bruker_noise_sigma: None,
            bruker_mobility_scale: Default::default(),
            bruker_streaming: false,
            bruker_halo: default_bruker_halo(),
            bruker_halo_peak_fraction: default_bruker_halo_peak_fraction(),
            bruker_halo_mz_idx_half_width: default_bruker_halo_mz_idx_half_width(),
            bruker_halo_scan_half_width: default_bruker_halo_scan_half_width(),
            bruker_ms1_polygon: false,
            bruker_ms1_polygon_mz_pad: 0.0,
            bruker_ms1_polygon_im_pad: 0.0,
            noise_filter_sigma: None,
            decoy_mode: false,
            ms2_hills_enabled: false,
            mz_recalibration: default_mz_recalibration(),
            mz_recalibration_mz_bins: default_mz_recalibration_mz_bins(),
            mz_recalibration_rt_bins: default_mz_recalibration_rt_bins(),
            mz_recalibration_min_samples: default_mz_recalibration_min_samples(),
            mz_recalibration_tol_sigma_mult: default_mz_recalibration_tol_sigma_mult(),
            mz_recalibration_tol_floor_ppm: default_mz_recalibration_tol_floor_ppm(),
        }
    }
}

fn default_mz_recalibration() -> bool {
    true
}

fn default_mz_recalibration_mz_bins() -> usize {
    20
}

fn default_mz_recalibration_rt_bins() -> usize {
    8
}

fn default_mz_recalibration_min_samples() -> usize {
    50
}

fn default_mz_recalibration_tol_sigma_mult() -> f64 {
    // Validated by a joint sweep on both benchmark platforms (Orbitrap
    // PXD003881, timsTOF PXD070049): recall peaks at 5.0 on Orbitrap and at
    // 3.0 on Bruker, with CV flat across the whole range on both. 4.0 is the
    // balanced single default — within ~0.2pp of recall of either
    // platform's own optimum, at no precision cost.
    4.0
}

fn default_mz_recalibration_tol_floor_ppm() -> f64 {
    1.0
}

fn default_global_min_mz() -> f64 {
    0.0
}

fn default_global_max_mz() -> f64 {
    f64::INFINITY
}

fn default_bruker_filter_mz_half_width() -> u32 {
    2
}

fn default_bruker_filter_max_internal_gap() -> usize {
    1
}

fn default_bruker_filter_min_feature_length() -> usize {
    5
}

fn default_bruker_filter_num_iterations() -> usize {
    1
}

fn default_bruker_watershed_box_scan() -> u32 {
    10
}

fn default_bruker_watershed_box_mz_idx() -> u32 {
    3
}

fn default_bruker_watershed_max_tof_offset() -> u32 {
    10
}

fn default_bruker_halo() -> bool {
    true
}

fn default_bruker_halo_peak_fraction() -> f64 {
    0.15
}

fn default_bruker_halo_mz_idx_half_width() -> u32 {
    80
}

fn default_bruker_halo_scan_half_width() -> usize {
    2
}
