use serde::{Deserialize, Serialize};

use crate::alignment::AlignmentConfig;
use crate::lfq::LfqConfig;

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
#[serde(deny_unknown_fields)]
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
    /// isolation window) from mzML inputs. Written to `hills_ms2.{ext}`.
    /// Bruker .d MS2 frames are not yet supported.
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
    /// where it is noisy. Only honoured for ppm tolerances.
    #[serde(default)]
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
            noise_filter_sigma: None,
            decoy_mode: false,
            ms2_hills_enabled: false,
            mz_recalibration: false,
            mz_recalibration_mz_bins: default_mz_recalibration_mz_bins(),
            mz_recalibration_rt_bins: default_mz_recalibration_rt_bins(),
            mz_recalibration_min_samples: default_mz_recalibration_min_samples(),
            mz_recalibration_tol_sigma_mult: default_mz_recalibration_tol_sigma_mult(),
            mz_recalibration_tol_floor_ppm: default_mz_recalibration_tol_floor_ppm(),
        }
    }
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

fn default_max_isotope_log2_ratio() -> f64 {
    1.5
}

fn default_chain_predicted_intensity_gate() -> bool {
    true
}

fn default_min_scan_overlap() -> usize {
    3
}

fn default_exhaustive_min_isotope_score() -> f64 {
    0.0
}

fn default_cosine_anchor() -> String {
    // `seed` is the default: anchoring every isotope's chromatographic cosine
    // to the monoisotope seed (as biosaur2 / AlphaPept / Dinosaur all do) beat
    // the former `adjacent` (predecessor) anchor on the full 20-run PXD003881
    // cohort — recall 0.7933 -> 0.7964 (+0.31 pp, +1565 covered PSMs) with no
    // quant regression (median CV, MV rate, and HUMAN FPR all flat-to-better).
    // Set to `adjacent` to reproduce the pre-2026-07 paper feature output.
    "seed".to_string()
}

fn default_exhaustive_assembly() -> bool {
    true
}

fn default_sulfur_aware_scoring() -> bool {
    true
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

/// Resolved form of `FeaturesConfig::cosine_anchor` — which hill the
/// chromatographic-cosine gate is measured against during isotope-chain
/// extension. Parsed (case-insensitively) from the config string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CosineAnchor {
    /// Anchor to the immediate predecessor in the chain (legacy koth default).
    Adjacent,
    /// Anchor every isotope to the monoisotope seed hill.
    Seed,
}

impl CosineAnchor {
    /// Parse a `cosine_anchor` string strictly: an unrecognised value is a
    /// config error rather than a silent fallback.
    pub fn parse(s: &str) -> Result<Self, crate::error::KothError> {
        match s.to_ascii_lowercase().as_str() {
            "adjacent" => Ok(CosineAnchor::Adjacent),
            "seed" => Ok(CosineAnchor::Seed),
            other => Err(crate::error::KothError::ConfigError(format!(
                "invalid features.cosine_anchor `{other}`: expected \"adjacent\" or \"seed\""
            ))),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeaturesConfig {
    pub min_charge: u8,
    pub max_charge: u8,
    /// Per-extension chromatographic-cosine threshold used while building an
    /// isotope chain. The candidate hill's cosine vs the seed hill must clear
    /// this to be added; otherwise chain extension stops in that direction.
    pub min_chain_cosine: f64,
    pub left_max_decrease: f64,
    pub right_max_decrease: f64,
    pub max_isotopes: usize,
    /// Per-extension intensity-ratio gate. After a candidate hill clears
    /// `min_chain_cosine`, also check that its apex intensity vs the
    /// predecessor's matches the averagine ratio within ± this many log2
    /// units. Catches column-bleed / contaminant hills that co-elute (so
    /// pass cosine) but have wrong intensity for an isotope. Set to a
    /// very large value to effectively disable. Typical: 1.5 (within
    /// ~2.83× of expected).
    #[serde(default = "default_max_isotope_log2_ratio")]
    pub max_isotope_log2_ratio: f64,
    /// When true (default, legacy behaviour), isotope-chain extension stops as
    /// soon as the averagine-*predicted* intensity of the next isotope falls
    /// below the per-run noise floor (5th-pct hill intensity × 0.8), computed
    /// from the seed intensity *before the hill is searched for*. This cheaply
    /// trims decayed chain tails, but it also blocks legitimate low-abundance
    /// monoisotopic seeds from ever pairing their M+1: a seed already near the
    /// floor has a predicted M+1 below the floor, so the chain never forms and
    /// the feature collapses to charge 0. When false, this predicted-intensity
    /// break is skipped and chain extension is terminated purely by *evidence*:
    /// a missing hill (`find_neighbors` empty), `min_chain_cosine`, the
    /// `max_isotope_log2_ratio` intensity-ratio gate, `right_max_decrease`
    /// (found hill must be ≥ this fraction of the predecessor), and the
    /// averagine-template-length / `max_isotopes` caps. Setting this false
    /// recovers dim 2+/3+ features that have both isotope hills present but
    /// were never paired — recommended when feeding an MS1 search that needs
    /// feature depth. Re-benchmark FDR/quant when flipping.
    #[serde(default = "default_chain_predicted_intensity_gate")]
    pub chain_predicted_intensity_gate: bool,
    /// Minimum number of mutually-overlapping scans two hills must share
    /// before a chromatographic cosine is computed between them; below it the
    /// cosine is 0 and the isotope chain will not extend across that pair.
    /// Default 3 (the long-standing hard gate — keep for PXD003881).
    ///
    /// On fast gradients hills are only 3–5 scans wide, so a dim isotope hill
    /// frequently overlaps its monoisotope by just 1–2 scans and is rejected
    /// outright regardless of how well the shapes agree — the dominant cause of
    /// koth emitting a low-abundance monoisotope as an unpaired charge-0 feature
    /// on such data (biosaur2 requires only ≥1 shared scan). Lowering this to 2
    /// recovers those pairs; pair it with a shorter `hills.min_scans` (a 2-scan
    /// hill can never reach a 3-scan overlap) and re-benchmark quant, since
    /// short-overlap pairs are noisier.
    #[serde(default = "default_min_scan_overlap")]
    pub min_scan_overlap: usize,
    /// Experimental non-destructive isotope assembler (biosaur2 / AlphaPept
    /// style). When false (default) koth uses the legacy greedy resolver: one
    /// best candidate per seed, claimed all-or-nothing (a seed whose isotope
    /// hills are stolen collapses to charge 0). When true, every (seed, charge)
    /// hypothesis enters an over-complete pool, contested hills are claimed
    /// longest-envelope-first, and a candidate whose isotope hills are partly
    /// claimed is TRUNCATED to its free prefix and re-queued rather than dropped
    /// — recovering charge-2/3 features greedy loses to shorter, higher-cosine
    /// competitors. Now the DEFAULT: validated paper-safe on PXD003881
    /// (recall +0.7pp, cohort CV neutral) and a win for MS1-search depth
    /// downstream (+38 proteins in uno). Set to `false` for the legacy greedy
    /// assembler (byte-identical to the pre-2026-07 paper output).
    #[serde(default = "default_exhaustive_assembly")]
    pub exhaustive_assembly: bool,
    /// (exhaustive_assembly only) Minimum isotope-pattern (Bhattacharyya) score a
    /// candidate — original or truncated — must reach before it may *claim* its
    /// hills; below-bar candidates are dropped so their hills stay free for a
    /// better-fitting feature. Default 0.0 = no gate.
    #[serde(default = "default_exhaustive_min_isotope_score")]
    pub exhaustive_min_isotope_score: f64,
    /// (exhaustive_assembly only) When true, contested-hill claim priority is
    /// ordered by envelope length, then isotope-pattern score, then composite —
    /// so the best averagine fit wins a shared hill within a length class. Default
    /// false.
    #[serde(default)]
    pub exhaustive_isotope_priority: bool,
    /// Chromatographic-cosine **anchor** for isotope-chain extension: which hill
    /// each candidate isotope's cosine gate is measured against.
    ///
    /// `"adjacent"` (default, legacy): anchor to the immediate predecessor in
    /// the chain — for M+1 that is the seed, for M+k≥2 the previously-claimed
    /// isotope hill. This is koth's deliberate choice: chains drift in S/N as
    /// they extend from the seed, so anchoring far isotopes to the seed
    /// over-rejects them. Default = byte-identical to legacy koth.
    ///
    /// `"seed"`: anchor *every* isotope's cosine to the monoisotope seed hill
    /// (the convention biosaur2 / AlphaPept / Dinosaur all use). The m/z step
    /// target and the intensity-ratio predecessor still step from the immediate
    /// predecessor — only the cosine reference changes. Rejects a far isotope
    /// that co-elutes with its neighbour but not with the mono.
    ///
    /// The per-extension cosine that feeds the composite `mean_cosine` is
    /// computed against whichever reference this selects (so the composite is
    /// consistent with the gate). The reported `cosine_score` field and the
    /// exhaustive-resolver rescoring remain adjacent-style regardless.
    #[serde(default = "default_cosine_anchor")]
    pub cosine_anchor: String,
    /// Score isotope chains against multiple averagine templates that vary
    /// the sulfur atom count `{0, avg, avg+2, avg+4}` and keep the best fit.
    /// Corrects the systematic Bhattacharyya penalty on Cys/Met-rich
    /// peptides whose M+2 is elevated by ³⁴S (4.25%, +2 Da).
    /// Default true.
    #[serde(default = "default_sulfur_aware_scoring")]
    pub sulfur_aware_scoring: bool,
    /// Neutron (C13) mass in Da
    pub neutron_mass: f64,
    /// Drop features whose **isotope_score** (Bhattacharyya vs averagine) is
    /// below this. 0.0 = keep all.
    pub min_isotope_score: f64,
    /// Drop features whose **cosine_score** (mean chromatographic cosine of
    /// adjacent isotope hills) is below this. 0.0 = keep all.
    pub min_cosine_score: f64,
    /// Drop features whose **combined_score** (= isotope × cosine) is below
    /// this. 0.0 = keep all. All three filters are AND-ed.
    pub min_combined_score: f64,
}

impl Default for FeaturesConfig {
    fn default() -> Self {
        Self {
            min_charge: 1,
            max_charge: 7,
            min_chain_cosine: 0.5,
            left_max_decrease: 0.05,
            right_max_decrease: 0.05,
            max_isotopes: 6,
            max_isotope_log2_ratio: 1.5,
            chain_predicted_intensity_gate: true,
            min_scan_overlap: 3,
            exhaustive_assembly: default_exhaustive_assembly(),
            exhaustive_min_isotope_score: 0.0,
            exhaustive_isotope_priority: false,
            cosine_anchor: default_cosine_anchor(),
            sulfur_aware_scoring: true,
            neutron_mass: 1.003_354_835,
            min_isotope_score: 0.0,
            min_cosine_score: 0.0,
            min_combined_score: 0.0,
        }
    }
}

impl FeaturesConfig {
    /// Validate string-valued knobs that serde alone cannot check. Called after
    /// TOML deserialization so an unrecognised value is rejected at load time
    /// rather than silently falling back at runtime.
    pub fn validate(&self) -> Result<(), crate::error::KothError> {
        CosineAnchor::parse(&self.cosine_anchor)?;
        Ok(())
    }

    /// Resolve the chromatographic-cosine anchor. Infallible at runtime because
    /// `validate()` has already rejected bad values at config load; the
    /// defensive fallback preserves legacy behaviour if reached.
    pub fn cosine_anchor_mode(&self) -> CosineAnchor {
        CosineAnchor::parse(&self.cosine_anchor).unwrap_or(CosineAnchor::Adjacent)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum OutputFormat {
    #[default]
    Tsv,
    Parquet,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct KothConfig {
    pub file: FileConfig,
    pub hills: HillsConfig,
    pub features: FeaturesConfig,
    pub scoring: ScoringConfig,
    pub output: OutputConfig,
}

/// Output settings for the alignment + LFQ stage.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlignOutputConfig {
    /// "tsv" or "parquet"
    pub format: OutputFormat,
    /// Only write matrix entries with q-value ≤ this threshold (1.0 = keep all)
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
            format: OutputFormat::Tsv,
            max_qvalue: 1.0,
            export_decoys: default_export_decoys(),
            export_details: default_export_details(),
        }
    }
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
        parse_koth("example_config.toml", include_str!("../../example_config.toml"));
        parse_koth("koth_ff.toml", include_str!("../../benchmark/config/koth_ff.toml"));
        parse_koth("koth_ff_bruker.toml", include_str!("../../benchmark/config/koth_ff_bruker.toml"));
        parse_koth("koth_ff_relaxed.toml", include_str!("../../benchmark/config/koth_ff_relaxed.toml"));
        parse_koth("koth_ff_sulfur_on.toml", include_str!("../../benchmark/config/koth_ff_sulfur_on.toml"));
        parse_koth("koth_ff_sulfur_off.toml", include_str!("../../benchmark/config/koth_ff_sulfur_off.toml"));
        parse_koth("koth_ff_alphapept_like.toml", include_str!("../../benchmark/config/koth_ff_alphapept_like.toml"));
    }

    #[test]
    fn shipped_align_configs_parse() {
        parse_align("example_config_align.toml", include_str!("../../example_config_align.toml"));
        parse_align("koth_align.toml", include_str!("../../benchmark/config/koth_align.toml"));
        parse_align("koth_align_bruker.toml", include_str!("../../benchmark/config/koth_align_bruker.toml"));
    }
}
