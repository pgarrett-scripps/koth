use rustc_hash::FxHashMap;
use std::sync::Arc;

use crate::config::{FileConfig, HillsConfig, ImToleranceType, ToleranceType};
use crate::models::{Hill, Spectrum};

use super::active::ActiveHill;
use super::calibration::MzDeltaHistogram;
use super::smooth;

/// Upper edge of the ppm-delta histogram used during adaptive-tolerance
/// calibration. 50 ppm is loose enough to capture even mis-calibrated
/// instruments; samples beyond it land in the overflow bucket and are
/// excluded from the median/σ calculation.
const CALIBRATION_HIST_MAX_PPM: f64 = 50.0;
/// Histogram bin width. 0.01 ppm gives 5000 bins (~20 KB) — plenty of
/// resolution for sub-ppm calibration.
const CALIBRATION_HIST_BIN_PPM: f64 = 0.01;

pub struct HillDetector {
    min_mz: f64,
    max_mz: f64,
    use_ppm: bool,
    tol_mult: f64,
    check_freq: usize,
    max_gap: usize,
    min_scans: usize,
    intensity_coverage: f64,
    im_tolerance: f64,
    im_tolerance_type: ImToleranceType,
    lfc_weight: f64,
    gap_fill_enabled: bool,
    smoothing_enabled: bool,
    smoothing_window: usize,
    /// Isolation window tag applied to every hill emitted by this detector.
    /// `None` for MS1; `Some` when running MS2 hill detection for a DIA channel.
    isolation_window: Option<crate::models::IsolationWindow>,
    active_hills: FxHashMap<usize, ActiveHill>,
    next_id: usize,
    finalized: Vec<Hill>,
    scan_idx: usize,
    /// Scratch buffer: (last_intensity, hill_id) sorted descending by intensity.
    /// Rebuilt each scan — highest-intensity hills get first pick of peaks.
    intensity_sorted: Vec<(f64, usize)>,
    /// Per-scan claimed-peak flags; resized each scan to reuse allocation.
    claimed_peaks: Vec<bool>,
    /// Histogram of accepted peak-to-hill |ppm| deltas, populated only
    /// when adaptive-tolerance calibration is enabled (pass 1).
    mz_delta_histogram: Option<MzDeltaHistogram>,
    /// Reused scratch for stale hill-ids each scan (avoids per-scan alloc).
    stale_buf: Vec<usize>,
    /// When true, accumulate the per-substage `process_scan` timings below.
    /// Gated on `Debug` logging being enabled so the four `Instant::now()`
    /// calls per scan aren't paid on every production run.
    profiling: bool,
    /// Per-substage timing accumulators (process_scan profiling).
    dbg_clean: std::time::Duration,
    dbg_sort: std::time::Duration,
    dbg_match: std::time::Duration,
    dbg_new: std::time::Duration,
}

impl HillDetector {
    pub fn new(config: &HillsConfig, file: &FileConfig) -> Self {
        let use_ppm = matches!(file.mz_tolerance_type, ToleranceType::Ppm);
        let tol_mult = if use_ppm {
            file.mz_tolerance / 1e6
        } else {
            file.mz_tolerance
        };
        // Stale-hill cleanup runs every scan. Skipping scans (the old
        // `check_freq = max_gap + 1` shortcut) lets a stale hill pick up a
        // peak between cleanups and end up with a gap longer than `max_gap`
        // in its profile.
        let check_freq = 1;
        Self {
            min_mz: file.global_min_mz,
            max_mz: file.global_max_mz,
            use_ppm,
            tol_mult,
            check_freq,
            max_gap: config.max_gap,
            min_scans: config.min_scans,
            intensity_coverage: file.intensity_coverage,
            im_tolerance: file.im_tolerance,
            im_tolerance_type: file.im_tolerance_type.clone(),
            lfc_weight: config.lfc_weight,
            gap_fill_enabled: config.gap_fill_enabled,
            smoothing_enabled: config.smoothing_enabled,
            smoothing_window: config.smoothing_window,
            isolation_window: None,
            active_hills: FxHashMap::default(),
            next_id: 0,
            finalized: Vec::new(),
            scan_idx: 0,
            intensity_sorted: Vec::new(),
            claimed_peaks: Vec::new(),
            mz_delta_histogram: None,
            stale_buf: Vec::new(),
            profiling: log::log_enabled!(log::Level::Debug),
            dbg_clean: std::time::Duration::ZERO,
            dbg_sort: std::time::Duration::ZERO,
            dbg_match: std::time::Duration::ZERO,
            dbg_new: std::time::Duration::ZERO,
        }
    }

    /// Tag every finalized hill with this isolation window (used for MS2/DIA).
    pub fn with_isolation_window(mut self, iw: crate::models::IsolationWindow) -> Self {
        self.isolation_window = Some(iw);
        self
    }

    /// Enable ppm-delta recording on every accepted peak-to-hill match.
    /// Use during a pass-1 sweep to gather the calibration distribution,
    /// then call `calibrated_tolerance_ppm` after `finish` (or before) to
    /// derive the pass-2 tolerance.
    pub fn with_calibration_recording(mut self) -> Self {
        self.mz_delta_histogram = Some(MzDeltaHistogram::new(
            CALIBRATION_HIST_MAX_PPM,
            CALIBRATION_HIST_BIN_PPM,
        ));
        self
    }

    /// Median + `sigma_mult × σ` of the recorded ppm-deltas, in ppm.
    /// Returns `None` if calibration recording was not enabled or fewer
    /// than two samples landed in range.
    pub fn calibrated_tolerance_ppm(&self, sigma_mult: f64) -> Option<f64> {
        self.mz_delta_histogram
            .as_ref()
            .and_then(|h| h.calibrated_ppm(sigma_mult))
    }

    /// Total in-range ppm-delta samples recorded during pass 1. Useful for
    /// logging and for deciding whether the calibration is statistically
    /// meaningful before applying it.
    pub fn calibration_sample_count(&self) -> u64 {
        self.mz_delta_histogram
            .as_ref()
            .map(|h| h.count_in_range())
            .unwrap_or(0)
    }

    pub fn process_scan(&mut self, spectrum: &Spectrum) {
        let rt = spectrum.retention_time;
        let use_im = spectrum.has_ion_mobility();
        let scan_idx = self.scan_idx;
        let n_peaks = spectrum.peaks.len();
        let _t_clean = self.profiling.then(std::time::Instant::now);

        // Single pass over the active set: collect stale hills AND build the
        // intensity-sorted index together, so we don't iterate every active
        // hill twice per scan. Stale-hill cleanup runs before matching so a hill
        // that has already exceeded max_gap cannot pick up a new peak and
        // accumulate extra gap zeros. A hill last seen at scan L tolerates
        // `max_gap` missed scans, so it is stale once `L < scan_idx - max_gap - 1`.
        // Buffers are moved out via take() so the loop can borrow active_hills
        // immutably while pushing to the (now-local) scratch vectors.
        let cutoff = scan_idx.saturating_sub(self.max_gap + 1);
        let do_clean = scan_idx % self.check_freq == 0;
        let mut sorted = std::mem::take(&mut self.intensity_sorted);
        let mut stale = std::mem::take(&mut self.stale_buf);
        sorted.clear();
        stale.clear();
        for (&id, h) in &self.active_hills {
            if do_clean && h.last_scan_seen < cutoff {
                stale.push(id);
            } else {
                sorted.push((h.last_real_intensity().unwrap_or(0.0), id));
            }
        }
        for &id in &stale {
            let hill = self.active_hills.remove(&id).unwrap();
            if let Some(h) = Self::finalize_hill(
                hill,
                self.min_scans,
                self.intensity_coverage,
                self.gap_fill_enabled,
                self.smoothing_enabled,
                self.smoothing_window,
                self.isolation_window,
            ) {
                self.finalized.push(h);
            }
        }
        if let Some(t) = _t_clean { self.dbg_clean += t.elapsed(); }
        let _t_sort = self.profiling.then(std::time::Instant::now);
        // Higher-intensity hills get first pick of candidate peaks, preventing
        // noise spikes from stealing matches away from real signal.
        sorted.sort_unstable_by(|a, b| {
            b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal)
        });
        self.intensity_sorted = sorted;
        self.stale_buf = stale;

        if let Some(t) = _t_sort { self.dbg_sort += t.elapsed(); }
        let _t_match = self.profiling.then(std::time::Instant::now);
        // Per-scan claimed-peak flags (reuse allocation across scans).
        self.claimed_peaks.clear();
        self.claimed_peaks.resize(n_peaks, false);

        // Process active hills in descending-intensity order. For each hill,
        // binary-search the mz-sorted peak array for unclaimed candidates and
        // pick the one with the lowest combined distance score.
        for &(_, hill_id) in &self.intensity_sorted {
            // Extract hill properties into locals to drop the immutable borrow
            // before the mutable add_point call below.
            let (hill_mz, hill_im, last_int) = {
                let h = &self.active_hills[&hill_id];
                (h.running_mz_mean, h.running_im_mean, h.last_real_intensity())
            };

            let tol = if self.use_ppm { hill_mz * self.tol_mult } else { self.tol_mult };
            let im_tol = self.calc_im_tol(hill_im);
            let lo = hill_mz - tol;
            let hi = hill_mz + tol;

            // Binary search into the mz-sorted peak array.
            let start = spectrum.peaks.partition_point(|p| (p.mz as f64) < lo);

            let mut best_peak_idx: Option<usize> = None;
            let mut best_dist = f64::INFINITY;

            for (rel, peak) in spectrum.peaks[start..].iter().enumerate() {
                let peak_mz = peak.mz as f64;
                if peak_mz > hi {
                    break;
                }
                let abs_idx = start + rel;
                if self.claimed_peaks[abs_idx] {
                    continue;
                }

                let mz_dist = (peak_mz - hill_mz).abs();
                let peak_im = peak.ion_mobility as f64;
                let dist = if use_im && im_tol.is_finite() && peak_im != 0.0 {
                    let im_dist = (peak_im - hill_im).abs();
                    if im_dist > im_tol {
                        continue;
                    }
                    mz_dist / tol + im_dist / im_tol
                } else {
                    mz_dist / tol
                };
                let lfc_term = if self.lfc_weight > 0.0 {
                    if let Some(prev) = last_int {
                        let lfc = (peak.intensity as f64 / prev).log2().abs();
                        self.lfc_weight * lfc
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                let dist = dist + lfc_term;
                if dist < best_dist {
                    best_dist = dist;
                    best_peak_idx = Some(abs_idx);
                }
            }

            if let Some(pk_idx) = best_peak_idx {
                let peak = &spectrum.peaks[pk_idx];
                if let Some(hist) = self.mz_delta_histogram.as_mut() {
                    if hill_mz > 0.0 {
                        let ppm = (peak.mz as f64 - hill_mz).abs() / hill_mz * 1e6;
                        hist.record(ppm);
                    }
                }
                self.active_hills
                    .get_mut(&hill_id)
                    .unwrap()
                    .add_point(
                        peak.mz as f64,
                        peak.intensity as f64,
                        rt,
                        peak.ion_mobility as f64,
                        scan_idx,
                    );
                self.claimed_peaks[pk_idx] = true;
            }
        }

        if let Some(t) = _t_match { self.dbg_match += t.elapsed(); }
        let _t_new = self.profiling.then(std::time::Instant::now);
        // Start new hills for any unclaimed peaks within the mz bounds.
        for (pk_idx, peak) in spectrum.peaks.iter().enumerate() {
            if self.claimed_peaks[pk_idx] {
                continue;
            }
            let mz = peak.mz as f64;
            if mz < self.min_mz || mz > self.max_mz {
                continue;
            }
            let id = self.next_id;
            self.next_id += 1;
            self.active_hills.insert(
                id,
                ActiveHill::new(mz, peak.intensity as f64, rt, peak.ion_mobility as f64, scan_idx, use_im),
            );
        }

        if let Some(t) = _t_new { self.dbg_new += t.elapsed(); }
        self.scan_idx += 1;
    }

    pub fn finish(mut self) -> Vec<Hill> {
        if self.profiling {
            log::info!(
                "[process_scan profile] clean {:.2?}, sort {:.2?}, match {:.2?}, new {:.2?}",
                self.dbg_clean, self.dbg_sort, self.dbg_match, self.dbg_new
            );
        }
        for (_, hill) in self.active_hills {
            if let Some(h) = Self::finalize_hill(
                hill,
                self.min_scans,
                self.intensity_coverage,
                self.gap_fill_enabled,
                self.smoothing_enabled,
                self.smoothing_window,
                self.isolation_window,
            ) {
                self.finalized.push(h);
            }
        }
        let hill_count = self.finalized.len();
        log::info!(
            "Finalized {} hills (~{:.1} MB in Hill structs)",
            hill_count,
            hill_count as f64 * 136.0 / 1_048_576.0
        );
        self.finalized.shrink_to_fit();
        self.finalized
    }

    fn calc_im_tol(&self, im: f64) -> f64 {
        if im == 0.0 {
            return f64::INFINITY;
        }
        match self.im_tolerance_type {
            ImToleranceType::Absolute => self.im_tolerance,
            ImToleranceType::Relative => im * self.im_tolerance,
        }
    }

    fn finalize_hill(
        mut hill: ActiveHill,
        min_scans: usize,
        intensity_coverage: f64,
        gap_fill_enabled: bool,
        smoothing_enabled: bool,
        smoothing_window: usize,
        isolation_window: Option<crate::models::IsolationWindow>,
    ) -> Option<Hill> {
        hill.trim();
        hill.trim_to_coverage(intensity_coverage);

        // Count raw gaps before any gap-fill so skipped_scans reflects true missing scans.
        let valid_scans = hill.intensity_profile.iter().filter(|&&x| x > 0.0).count();
        if valid_scans < min_scans {
            return None;
        }
        let skipped = hill.intensity_profile.iter().filter(|&&x| x == 0.0).count();

        smooth::apply_intensity_filters(
            &mut hill.intensity_profile,
            gap_fill_enabled,
            smoothing_enabled,
            smoothing_window,
        );

        let apex_idx = hill.apex_index();
        let apex_scan = hill.scan_start + apex_idx;
        let mz_mean = hill.mz_weighted_mean();
        let mz_std = hill.mz_weighted_std();
        let mz_se = hill.mz_kish_se();
        let im_mean = hill.im_weighted_mean();
        let im_std = hill.im_weighted_std();
        let rt_start = hill.min_rt();
        let rt_end = hill.max_rt();
        let rt_apex = hill.apex_rt();
        let intensity_sum: f64 = hill.intensity_profile.iter().map(|&x| x as f64).sum();
        let intensity_max = hill
            .intensity_profile
            .iter()
            .map(|&x| x as f64)
            .fold(0.0f64, f64::max);
        let hill_score = smooth::compute_hill_score(&hill.intensity_profile);
        let n_scans = hill.intensity_profile.len();

        Some(Hill {
            hill_id: 0, // filled in by `assign_hill_ids` after splitting
            mz: mz_mean,
            mz_std,
            mz_se,
            rt: rt_apex,
            rt_start,
            rt_end,
            rt_width: rt_end - rt_start,
            im: im_mean,
            im_std,
            scan_start: hill.scan_start,
            scan_apex: apex_scan,
            scan_end: hill.scan_start + n_scans - 1,
            n_scans,
            skipped_scans: skipped,
            intensity_sum,
            intensity_max,
            hill_score,
            intensity_profile: Arc::from(hill.intensity_profile.as_slice()),
            isolation_window,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{FileConfig, HillsConfig};
    use crate::models::{Peak, Spectrum};

    fn spec(scan_index: usize, rt: f64, mzs: &[f32]) -> Spectrum {
        Spectrum {
            scan_index,
            retention_time: rt,
            peaks: mzs
                .iter()
                .map(|&mz| Peak {
                    mz,
                    intensity: 1000.0,
                    ion_mobility: 0.0,
                })
                .collect(),
            ms_level: 1,
            isolation_window: None,
        }
    }

    fn default_cfgs(max_gap: usize) -> (HillsConfig, FileConfig) {
        let mut h = HillsConfig::default();
        h.max_gap = max_gap;
        h.split_hills = false;
        h.min_scans = 3;
        (h, FileConfig::default())
    }

    /// Regression: with `max_gap = 0`, a hill that matches the same m/z in N
    /// consecutive scans must be retained (previously a stale-cleanup off-by-one
    /// killed every hill before it could extend past 1 scan).
    #[test]
    fn max_gap_zero_retains_consecutive_hill() {
        let (h, f) = default_cfgs(0);
        let mut det = HillDetector::new(&h, &f);
        for i in 0..5 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        let hills = det.finish();
        assert_eq!(hills.len(), 1, "expected 1 hill, got {}", hills.len());
        assert_eq!(hills[0].n_scans, 5);
        assert_eq!(hills[0].skipped_scans, 0);
    }

    /// With `max_gap = 0`, a one-scan gap must terminate the hill.
    #[test]
    fn max_gap_zero_disallows_gaps() {
        let (h, f) = default_cfgs(0);
        let mut det = HillDetector::new(&h, &f);
        for i in 0..3 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        det.process_scan(&spec(3, 0.3, &[]));
        for i in 4..7 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        let hills = det.finish();
        assert_eq!(hills.len(), 2, "expected 2 hills split across the gap, got {}", hills.len());
        for hl in &hills {
            assert_eq!(hl.skipped_scans, 0);
        }
    }

    /// With `max_gap = 1`, a single missed scan must be bridged.
    #[test]
    fn max_gap_one_bridges_single_miss() {
        let (h, f) = default_cfgs(1);
        let mut det = HillDetector::new(&h, &f);
        for i in 0..2 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        det.process_scan(&spec(2, 0.2, &[])); // miss
        for i in 3..5 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        let hills = det.finish();
        assert_eq!(hills.len(), 1, "expected 1 bridged hill, got {}", hills.len());
        assert_eq!(hills[0].skipped_scans, 1);
    }

    /// Calibration recording: when `with_calibration_recording` is set,
    /// every accepted peak-to-hill match contributes one |ppm| delta to
    /// the histogram. A clean signal (same m/z every scan) lands at 0 ppm
    /// so the calibrated tolerance approaches 0 — the calibration loop
    /// has nothing to widen against.
    #[test]
    fn calibration_records_zero_ppm_for_constant_signal() {
        let (h, f) = default_cfgs(0);
        let mut det = HillDetector::new(&h, &f).with_calibration_recording();
        for i in 0..10 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        // 9 accepted matches (peak in scans 1..9 each match the hill
        // started in scan 0). All are at exactly 500.0 m/z so |ppm|=0.
        assert_eq!(det.calibration_sample_count(), 9);
        let cal = det.calibrated_tolerance_ppm(3.0).unwrap();
        assert!(cal < 0.1, "expected calibrated ~0, got {cal}");
    }

    /// Calibration with a drift: m/z walks up by ~1 ppm per scan. The
    /// histogram median + 3σ should land somewhere comfortably above 0
    /// (the drift magnitude × sigma_mult).
    #[test]
    fn calibration_widens_for_drifting_signal() {
        let (h, f) = default_cfgs(0);
        let mut det = HillDetector::new(&h, &f).with_calibration_recording();
        // 1 ppm/scan walk on a 500 Da peak = +0.0005 Da/scan.
        for i in 0..20 {
            let mz = 500.0 + (i as f32) * 5e-4;
            det.process_scan(&spec(i, i as f64 * 0.1, &[mz]));
        }
        assert!(det.calibration_sample_count() > 0);
        let cal = det.calibrated_tolerance_ppm(3.0).unwrap();
        assert!(
            cal > 0.5,
            "expected calibrated > 0.5 ppm for drifting signal, got {cal}"
        );
    }

    /// With `max_gap = 1`, two consecutive missed scans must split the hill.
    #[test]
    fn max_gap_one_splits_at_two_misses() {
        let (h, f) = default_cfgs(1);
        let mut det = HillDetector::new(&h, &f);
        for i in 0..3 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        det.process_scan(&spec(3, 0.3, &[])); // miss
        det.process_scan(&spec(4, 0.4, &[])); // miss
        for i in 5..8 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        let hills = det.finish();
        assert_eq!(hills.len(), 2, "expected 2 hills, got {}", hills.len());
    }
}
