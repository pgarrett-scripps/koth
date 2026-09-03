use std::sync::Arc;

use serde::{Deserialize, Serialize};

/// A single centroided MS1 peak.
#[derive(Debug, Clone, Copy)]
pub struct Peak {
    pub mz: f32,
    pub intensity: f32,
    /// Ion mobility value (1/K0 for timsTOF). 0.0 means not available.
    pub ion_mobility: f32,
}

/// Precursor isolation window for an MS2 spectrum (mzML coordinates).
///
/// All values in m/z. `lower` and `upper` are absolute boundaries (already
/// expanded from any target offset in the source file).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct IsolationWindow {
    pub target: f64,
    pub lower: f64,
    pub upper: f64,
}

impl IsolationWindow {
    /// Stable integer key for grouping MS2 scans that share an isolation window.
    /// Encodes (target, lower, upper) to 4 decimal places.
    pub fn key(&self) -> (i64, i64, i64) {
        (
            (self.target * 10_000.0).round() as i64,
            (self.lower * 10_000.0).round() as i64,
            (self.upper * 10_000.0).round() as i64,
        )
    }
}

impl PartialEq for IsolationWindow {
    fn eq(&self, other: &Self) -> bool {
        self.key() == other.key()
    }
}

impl Eq for IsolationWindow {}

/// A single mass spectrum with all its peaks.
#[derive(Debug, Clone)]
pub struct Spectrum {
    /// 0-based sequential index (scan number in the run)
    pub scan_index: usize,
    /// Retention time in minutes
    pub retention_time: f64,
    /// Peaks sorted by mz ascending
    pub peaks: Vec<Peak>,
    /// 1 for MS1, 2 for MS2, etc.
    pub ms_level: u8,
    /// Precursor isolation window (only set for MS2 spectra)
    pub isolation_window: Option<IsolationWindow>,
    /// FAIMS compensation voltage for this scan, in volts. `None` when the
    /// acquisition does not use FAIMS or the source format does not expose it.
    pub faims_cv: Option<f32>,
}

impl Spectrum {
    pub fn has_ion_mobility(&self) -> bool {
        self.peaks.iter().any(|p| p.ion_mobility != 0.0)
    }
}

/// A finalized chromatographic hill (elution peak for a single m/z trace).
///
/// Column names match the Python zenith_feature_finder hills.tsv output exactly.
#[derive(Debug, Clone)]
pub struct Hill {
    /// Stable identifier assigned at the end of hill detection. Matches the
    /// `hill_id` column in `hills.{tsv,parquet}` and the entries in the
    /// `hill_ids` column of `features.{tsv,parquet}`. IDs are unique within a
    /// single hills file (MS1 and MS2 each have their own ID space).
    pub hill_id: u64,
    pub mz: f64,
    pub mz_std: f64,
    /// Standard error of the intensity-weighted m/z estimate (Kish-effective
    /// sample size). Zero for hills with <2 real peaks. Written to the hills
    /// output as a per-hill m/z-confidence diagnostic.
    pub mz_se: f64,
    pub rt: f64,
    pub rt_start: f64,
    pub rt_end: f64,
    pub rt_width: f64,
    pub im: f64,
    pub im_std: f64,
    pub scan_start: usize,
    pub scan_apex: usize,
    pub scan_end: usize,
    pub n_scans: usize,
    pub skipped_scans: usize,
    pub intensity_sum: f64,
    pub intensity_max: f64,
    /// Monotonicity-fraction shape score (0–1). 1.0 = perfectly hill-shaped.
    pub hill_score: f64,
    /// Per-scan intensity profile (length == n_scans, zeros where gaps)
    pub intensity_profile: Arc<[f32]>,
    /// Isolation window the hill was built from (None for MS1 hills,
    /// Some for MS2/DIA hills).
    pub isolation_window: Option<IsolationWindow>,
    /// FAIMS compensation voltage shared by every scan in this hill.
    /// `None` for non-FAIMS data and for formats without FAIMS metadata.
    pub faims_cv: Option<f32>,
}

impl Hill {
    /// Intensity at a specific absolute scan index, or 0.0 if out of range.
    pub fn intensity_at_scan(&self, scan_idx: usize) -> f64 {
        if scan_idx < self.scan_start || scan_idx > self.scan_end {
            return 0.0;
        }
        let rel = scan_idx - self.scan_start;
        if rel < self.intensity_profile.len() {
            self.intensity_profile[rel] as f64
        } else {
            0.0
        }
    }

    /// Returns true if this hill's scan range overlaps with another.
    pub fn overlaps(&self, other: &Hill) -> bool {
        !(self.scan_end < other.scan_start || other.scan_end < self.scan_start)
    }
}

/// A detected isotope feature: a group of hills forming an isotope envelope.
///
/// `hills` are sorted by mz (monoisotopic first).
/// `charge == 0` means unknown charge (single peak, no isotope partners found).
#[derive(Debug, Clone)]
pub struct Feature {
    pub hills: Vec<Hill>,
    pub charge: u8,
    /// Mean **chromatographic** cosine similarity (in scan/RT space) between
    /// adjacent isotope hills in the chain. Bounded [0, 1]. Measures whether
    /// the isotope hills co-elute in time — *not* whether their intensities
    /// match an averagine pattern.
    pub cosine_score: f64,
    pub ppm_error: f64,
}

impl Feature {
    pub fn monoisotopic_hill(&self) -> &Hill {
        &self.hills[0]
    }

    pub fn monoisotopic_mz(&self) -> f64 {
        self.hills[0].mz
    }

    /// FAIMS compensation voltage of this feature, inherited from its
    /// monoisotopic hill. Feature assembly guarantees all member hills share
    /// this value.
    pub fn faims_cv(&self) -> Option<f32> {
        self.hills.first().and_then(|h| h.faims_cv)
    }

    /// Neutral monoisotopic mass (mass = mz * z - z * proton_mass)
    pub fn monoisotopic_neutral_mass(&self) -> Option<f64> {
        if self.charge == 0 {
            return None;
        }
        const PROTON_MASS: f64 = 1.007_276_466_621;
        Some(self.monoisotopic_mz() * self.charge as f64 - self.charge as f64 * PROTON_MASS)
    }

    pub fn rt_apex(&self) -> f64 {
        // RT at the hill with the highest intensity_max
        self.hills
            .iter()
            .max_by(|a, b| a.intensity_max.partial_cmp(&b.intensity_max).unwrap())
            .map(|h| h.rt)
            .unwrap_or(0.0)
    }

    pub fn rt_start(&self) -> f64 {
        self.hills
            .iter()
            .map(|h| h.rt_start)
            .fold(f64::INFINITY, f64::min)
    }

    pub fn rt_end(&self) -> f64 {
        self.hills
            .iter()
            .map(|h| h.rt_end)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    pub fn im_apex(&self) -> f64 {
        self.hills
            .iter()
            .max_by(|a, b| a.intensity_max.partial_cmp(&b.intensity_max).unwrap())
            .map(|h| h.im)
            .unwrap_or(0.0)
    }

    pub fn total_intensity(&self) -> f64 {
        self.hills.iter().map(|h| h.intensity_sum).sum()
    }

    pub fn n_scans_total(&self) -> usize {
        self.hills.iter().map(|h| h.n_scans).sum()
    }

    /// Combined elution profile: sum of all isotope hill intensities per scan.
    /// Returns (min_scan, max_scan, profile_vec).
    pub fn elution_profile(&self) -> (usize, usize, Vec<f64>) {
        let min_scan = self.hills.iter().map(|h| h.scan_start).min().unwrap_or(0);
        let max_scan = self.hills.iter().map(|h| h.scan_end).max().unwrap_or(0);
        let scan_range = max_scan - min_scan + 1;
        let mut profile = vec![0.0f64; scan_range];
        for hill in &self.hills {
            for (i, &intensity) in hill.intensity_profile.iter().enumerate() {
                let idx = hill.scan_start + i - min_scan;
                if idx < scan_range {
                    profile[idx] += intensity as f64;
                }
            }
        }
        (min_scan, max_scan, profile)
    }

    /// Total intensity at the apex scan (summed across all isotope hills).
    pub fn total_intensity_at_apex(&self) -> f64 {
        let (_, _, profile) = self.elution_profile();
        profile.iter().cloned().fold(0.0f64, f64::max)
    }

    /// Apex scan index: scan where the sum of all isotope hill intensities is highest.
    pub fn apex_scan(&self) -> usize {
        let (min_scan, _, profile) = self.elution_profile();
        profile
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| min_scan + i)
            .unwrap_or(0)
    }

    /// Isotope profile at the apex scan: one intensity per hill.
    pub fn isotope_profile_apex(&self) -> Vec<f64> {
        let apex = self.apex_scan();
        self.hills
            .iter()
            .map(|h| h.intensity_at_scan(apex))
            .collect()
    }

    /// Scan list for monoisotopic hill.
    pub fn mono_scan_list(&self) -> Vec<usize> {
        let h = &self.hills[0];
        (h.scan_start..=h.scan_end).collect()
    }

    /// Intensity list for monoisotopic hill.
    pub fn mono_intensity_list(&self) -> &[f32] {
        &self.hills[0].intensity_profile
    }
}

/// Feature with isotope pattern scoring applied.
///
/// Carries three distinct quality scores so downstream stages (alignment,
/// consensus grouping, LFQ filters) can filter on whichever makes sense for
/// their purpose. All three are bounded `[0, 1]`; higher = better.
#[derive(Debug, Clone)]
pub struct ScoredFeature {
    pub feature: Feature,
    /// Neutron offset applied: -1, 0, or 1.
    /// Non-zero means the observed monoisotopic peak is actually M+|offset|.
    pub neutron_offset: i8,
    /// **Isotope-pattern** score: Bhattacharyya coefficient between the
    /// observed apex-scan isotope intensities and the theoretical averagine
    /// template. Measures "does the isotope envelope look like a peptide?".
    pub isotope_score: f64,
    /// **Chromatographic cosine** score: mirror of `feature.cosine_score` —
    /// the mean cosine similarity (in time domain) between adjacent isotope
    /// hills. Measures "do the isotope hills co-elute?".
    pub cosine_score: f64,
    /// **Combined** score: `isotope_score × cosine_score`. The default
    /// quality knob used by alignment, consensus grouping, and feature
    /// retention. High requires both a good pattern match AND tight co-elution.
    pub combined_score: f64,
    /// Theoretical averagine isotope pattern (normalized to sum=1)
    pub theoretical_pattern: Vec<f64>,
}

impl ScoredFeature {
    /// Corrected monoisotopic m/z (adjusted by neutron_offset).
    pub fn monoisotopic_mz(&self) -> f64 {
        if self.feature.charge == 0 {
            return self.feature.monoisotopic_mz();
        }
        const C13_NEUTRON: f64 = 1.003_354_835;
        self.feature.monoisotopic_mz()
            - (self.neutron_offset as f64 * C13_NEUTRON / self.feature.charge as f64)
    }

    /// Corrected neutral monoisotopic mass.
    pub fn monoisotopic_neutral_mass(&self) -> Option<f64> {
        let base = self.feature.monoisotopic_neutral_mass()?;
        const C13_NEUTRON: f64 = 1.003_354_835;
        Some(base - self.neutron_offset as f64 * C13_NEUTRON)
    }
}
