use std::collections::HashMap;
use std::sync::Arc;

use crate::config::{FileConfig, HillsConfig, ImToleranceType, ToleranceType};
use crate::models::{Hill, Spectrum};

use super::active::ActiveHill;

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
    active_hills: HashMap<usize, ActiveHill>,
    next_id: usize,
    finalized: Vec<Hill>,
    scan_idx: usize,
    /// Scratch buffer: (last_intensity, hill_id) sorted descending by intensity.
    /// Rebuilt each scan — highest-intensity hills get first pick of peaks.
    intensity_sorted: Vec<(f64, usize)>,
    /// Per-scan claimed-peak flags; resized each scan to reuse allocation.
    claimed_peaks: Vec<bool>,
}

impl HillDetector {
    pub fn new(config: &HillsConfig, file: &FileConfig) -> Self {
        let use_ppm = matches!(file.mz_tolerance_type, ToleranceType::Ppm);
        let tol_mult = if use_ppm {
            file.mz_tolerance / 1e6
        } else {
            file.mz_tolerance
        };
        let check_freq = if config.max_gap == 0 {
            1
        } else {
            (config.max_gap + 1).max(1)
        };
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
            active_hills: HashMap::new(),
            next_id: 0,
            finalized: Vec::new(),
            scan_idx: 0,
            intensity_sorted: Vec::new(),
            claimed_peaks: Vec::new(),
        }
    }

    pub fn process_scan(&mut self, spectrum: &Spectrum) {
        let rt = spectrum.retention_time;
        let use_im = spectrum.has_ion_mobility();
        let scan_idx = self.scan_idx;
        let n_peaks = spectrum.peaks.len();

        // Build intensity-sorted (descending) index of active hills.
        // Higher-intensity hills get first pick of candidate peaks, preventing
        // noise spikes from stealing matches away from real signal.
        self.intensity_sorted.clear();
        self.intensity_sorted.extend(
            self.active_hills
                .iter()
                .map(|(&id, h)| (h.last_real_intensity().unwrap_or(0.0), id)),
        );
        self.intensity_sorted.sort_unstable_by(|a, b| {
            b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal)
        });

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
                ActiveHill::new(mz, peak.intensity as f64, rt, peak.ion_mobility as f64, scan_idx),
            );
        }

        // Periodic stale-hill cleanup.
        if scan_idx % self.check_freq == 0 {
            let cutoff = scan_idx.saturating_sub(self.max_gap);
            let stale_ids: Vec<usize> = self
                .active_hills
                .iter()
                .filter(|(_, h)| h.last_scan_seen < cutoff)
                .map(|(&id, _)| id)
                .collect();

            for id in stale_ids {
                let hill = self.active_hills.remove(&id).unwrap();
                if let Some(h) = Self::finalize_hill(hill, self.min_scans, self.intensity_coverage) {
                    self.finalized.push(h);
                }
            }
        }

        self.scan_idx += 1;
    }

    pub fn finish(mut self) -> Vec<Hill> {
        for (_, hill) in self.active_hills {
            if let Some(h) = Self::finalize_hill(hill, self.min_scans, self.intensity_coverage) {
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

    fn finalize_hill(mut hill: ActiveHill, min_scans: usize, intensity_coverage: f64) -> Option<Hill> {
        hill.trim();
        hill.trim_to_coverage(intensity_coverage);
        if hill.len() < min_scans {
            return None;
        }

        let apex_idx = hill.apex_index();
        let apex_scan = hill.scan_start + apex_idx;
        let mz_mean = hill.mz_weighted_mean();
        let mz_std = hill.mz_weighted_std();
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
        let skipped = hill.intensity_profile.iter().filter(|&&x| x == 0.0).count();
        let n_scans = hill.intensity_profile.len();

        Some(Hill {
            mz: mz_mean,
            mz_std,
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
            intensity_profile: Arc::from(hill.intensity_profile.as_slice()),
        })
    }
}
