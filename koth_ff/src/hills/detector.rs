use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use crate::config::{HillsConfig, ImToleranceType, ToleranceType};
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
    active_hills: HashMap<usize, ActiveHill>,
    next_id: usize,
    finalized: Vec<Hill>,
    matched_ids: HashSet<usize>,
    scan_idx: usize,
    /// Scratch buffer rebuilt each scan: sorted (mz_mean, hill_id).
    /// Reusing the allocation avoids a malloc per scan.
    mz_sorted: Vec<(f64, usize)>,
}

impl HillDetector {
    pub fn new(config: &HillsConfig) -> Self {
        let use_ppm = matches!(config.mz_tolerance_type, ToleranceType::Ppm);
        let tol_mult = if use_ppm {
            config.mz_tolerance / 1e6
        } else {
            config.mz_tolerance
        };
        let check_freq = if config.max_gap == 0 {
            1
        } else {
            (config.max_gap + 1).max(1)
        };
        Self {
            min_mz: config.global_min_mz,
            max_mz: config.global_max_mz,
            use_ppm,
            tol_mult,
            check_freq,
            max_gap: config.max_gap,
            min_scans: config.min_scans,
            intensity_coverage: config.intensity_coverage,
            im_tolerance: config.im_tolerance,
            im_tolerance_type: config.im_tolerance_type.clone(),
            active_hills: HashMap::new(),
            next_id: 0,
            finalized: Vec::new(),
            matched_ids: HashSet::new(),
            scan_idx: 0,
            mz_sorted: Vec::new(),
        }
    }

    pub fn process_scan(&mut self, spectrum: &Spectrum) {
        let rt = spectrum.retention_time;
        let use_im = spectrum.has_ion_mobility();
        let scan_idx = self.scan_idx;

        // Build sorted (mz_mean, id) index for this scan.
        // Peaks are already sorted by mz, so binary-searching this vec
        // replaces the entire bin HashMap with a single sort + bisect.
        self.mz_sorted.clear();
        self.mz_sorted
            .extend(self.active_hills.iter().map(|(&id, h)| (h.running_mz_mean, id)));
        self.mz_sorted
            .sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        self.matched_ids.clear();

        for peak in &spectrum.peaks {
            let mz = peak.mz as f64;
            let im = peak.ion_mobility as f64;

            if mz < self.min_mz {
                continue;
            }
            if mz > self.max_mz {
                break;
            }

            let tol = if self.use_ppm { mz * self.tol_mult } else { self.tol_mult };
            let im_tol = self.calc_im_tol(im);
            let lo = mz - tol;
            let hi = mz + tol;

            // Binary search for the first active hill with mz >= lo,
            // then scan forward while mz <= hi.
            let start = self.mz_sorted.partition_point(|&(m, _)| m < lo);

            let mut best_id: Option<usize> = None;
            let mut best_dist = f64::INFINITY;

            for &(hill_mz, id) in &self.mz_sorted[start..] {
                if hill_mz > hi {
                    break;
                }
                if self.matched_ids.contains(&id) {
                    continue;
                }
                let mz_dist = (mz - hill_mz).abs();
                let dist = if use_im && im_tol.is_finite() {
                    let hill_im = self.active_hills[&id].running_im_mean;
                    let im_dist = (im - hill_im).abs();
                    if im_dist > im_tol {
                        continue;
                    }
                    mz_dist / tol + im_dist / im_tol
                } else {
                    mz_dist / tol
                };
                if dist < best_dist {
                    best_dist = dist;
                    best_id = Some(id);
                }
            }

            if let Some(hill_id) = best_id {
                self.active_hills
                    .get_mut(&hill_id)
                    .unwrap()
                    .add_point(mz, peak.intensity as f64, rt, im, scan_idx);
                self.matched_ids.insert(hill_id);
            } else {
                let id = self.next_id;
                self.next_id += 1;
                self.active_hills.insert(
                    id,
                    ActiveHill::new(mz, peak.intensity as f64, rt, im, scan_idx),
                );
            }
        }

        // Periodic stale-hill cleanup
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
