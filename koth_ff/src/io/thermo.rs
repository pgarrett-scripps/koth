//! Native Thermo Fisher `.raw` reader (`thermo` cargo feature, on by default).
//!
//! Built on [`opentfraw`], a pure-Rust parser of the Thermo `.raw` binary
//! format (file versions 8–66). There is no .NET runtime, no vendor assembly
//! and no system library: a `.raw` file is read like any other local file.
//!
//! The MS1 output matches [`crate::io::mzml::read_mzml`] on an mzML produced by
//! msconvert with vendor peak picking: **MS1 spectra only**, **centroided**
//! (koth's hill detector consumes centroids and does no peak picking of its
//! own), peaks filtered to positive intensity and sorted by m/z ascending,
//! spectra sorted by retention time and re-indexed `0..n`. Ion mobility is
//! always 0.0 (Orbitrap has no IM dimension).
//!
//! **Centroids.** FT (Orbitrap) and Astral scans store the instrument's own
//! centroid list next to the profile signal; that list is what msconvert's
//! vendor peak picking and Thermo's `RawFileReader` return, and it is what this
//! reader uses. A scan that carries only a profile signal (older ion-trap or
//! LTQ profile acquisitions) is centroided here with mzdata's quadratic peak
//! picker at S/N >= 1, exactly as the mzML path centroids profile-mode mzML.
//!
//! [`read_thermo_ms2`] is the opt-in **DIA MS2** counterpart: it emits one MS2
//! [`Spectrum`] per MS2 scan, stamped with the scan's precursor isolation window,
//! mirroring the mzML MS2 path ([`crate::io::mzml::stream_mzml_ms2`]). It is
//! **DIA-only** — a DDA `.raw` (or one it cannot confidently classify as DIA)
//! yields an empty set plus a warning and never disturbs MS1. See
//! `is_dia_schedule` for the heuristic and its failure modes.

use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use opentfraw::generic_data::GenericValue;
use opentfraw::{MsPower, RawFileReader};

use crate::error::KothError;
use crate::models::{IsolationWindow, Peak, Spectrum};

#[cfg(test)]
#[path = "thermo_tests.rs"]
mod tests;

/// An opened `.raw`: the parsed index/metadata plus a buffered handle for
/// decoding scan packets on demand.
pub(crate) struct RawSource {
    raw: RawFileReader,
    file: BufReader<File>,
    /// Whether the decoded scan events line up with the scans (see
    /// [`RawSource::check_event_alignment`]).
    events_aligned: bool,
}

impl RawSource {
    pub(crate) fn open(path: &Path) -> Result<Self, KothError> {
        let open_err = |e: &dyn std::fmt::Display| {
            KothError::ThermoError(format!(
                "could not open Thermo .raw file '{}': {e}",
                path.display()
            ))
        };
        let raw = RawFileReader::open_path(path).map_err(|e| open_err(&e))?;
        let file = File::open(path).map_err(|e| open_err(&e))?;
        let mut source = Self {
            raw,
            file: BufReader::new(file),
            events_aligned: true,
        };
        source.events_aligned = source.check_event_alignment();
        if !source.events_aligned {
            log::warn!(
                "Thermo '{}': scan events are inconsistent with the per-scan trailer; \
                 MS levels come from the trailer's Master Scan Number and scan-event \
                 values (precursor m/z, profile calibration) are not used",
                path.display()
            );
        }
        Ok(source)
    }

    /// Check the decoded scan events against the per-scan trailer. A scan the
    /// trailer marks as dependent (`Master Scan Number` > 0) must carry an MSn
    /// event. opentfraw 1.4 can decode the variable-length scan events of an
    /// Orbitrap Fusion file out of step with the scans (PXD003881: 33 127 of
    /// 66 254 dependent scans read as MS1 or undefined), while the trailer
    /// table, read record by record, stays correct. More than 0.1 % such
    /// scans marks the events as misaligned. The converse is not checked:
    /// Astral DIA scans record `Master Scan Number` 0 on MSn scans.
    fn check_event_alignment(&self) -> bool {
        if self.raw.flat_peaks {
            return true;
        }
        let (mut dependent, mut contradicted) = (0usize, 0usize);
        for idx in 0..self.n_scans() {
            if self.master_scan(idx).is_some_and(|m| m > 0) {
                dependent += 1;
                if self.event_level(idx) == 1 {
                    contradicted += 1;
                }
            }
        }
        contradicted * 1000 <= dependent
    }

    pub(crate) fn n_scans(&self) -> u32 {
        self.raw.num_scans
    }

    /// 1-based Thermo scan number of the zero-based scan index `idx`.
    fn scan_number(&self, idx: u32) -> u32 {
        self.raw.run_header.sample_info.first_scan_number + idx
    }

    /// MS level recorded in scan `idx`'s scan event (1 when absent or
    /// undefined).
    fn event_level(&self, idx: u32) -> u32 {
        match self
            .raw
            .scan_events
            .get(idx as usize)
            .and_then(|e| e.preamble.ms_power())
        {
            None | Some(MsPower::Undefined) | Some(MsPower::Ms1) => 1,
            Some(MsPower::Ms2) => 2,
            Some(_) => 3,
        }
    }

    /// MS level of scan `idx` (no signal decode). TSQ/SRM files (flat peak
    /// format) are all MS2 transitions. When the scan events are misaligned the
    /// trailer's `Master Scan Number` decides (0 = survey scan, else MS2).
    pub(crate) fn ms_level(&self, idx: u32) -> u32 {
        if self.raw.flat_peaks {
            return 2;
        }
        if self.events_aligned {
            return self.event_level(idx);
        }
        match self.master_scan(idx) {
            Some(m) if m > 0 => 2,
            _ => 1,
        }
    }

    /// Trailer `Master Scan Number` of scan `idx`, if recorded.
    fn master_scan(&self, idx: u32) -> Option<i32> {
        self.raw
            .scan_params(self.scan_number(idx))
            .and_then(|p| p.master_scan_number())
    }

    /// Scan start time in minutes, from the scan index (no signal decode).
    pub(crate) fn retention_time(&self, idx: u32) -> f64 {
        self.raw
            .scan_index
            .get(idx as usize)
            .map_or(f64::NAN, |e| e.start_time)
    }

    /// FAIMS compensation voltage from the per-scan trailer, if recorded.
    pub(crate) fn faims_cv(&self, idx: u32) -> Option<f32> {
        let record = self.raw.scan_parameters(self.scan_number(idx))?;
        record
            .values
            .iter()
            .find(|(label, _)| label.trim_end_matches(':').trim() == "FAIMS CV")
            .and_then(|(_, value)| faims_value(value))
    }

    /// Decode scan `idx` to centroided koth peaks.
    pub(crate) fn peaks(&mut self, idx: u32) -> Result<Vec<Peak>, KothError> {
        let scan_number = self.scan_number(idx);
        let read_err = |e: opentfraw::Error| {
            KothError::ThermoError(format!("could not read scan {scan_number}: {e}"))
        };
        let raw_peaks = self
            .raw
            .read_peaks_only(&mut self.file, scan_number)
            .map_err(read_err)?;
        if !raw_peaks.is_empty() || self.raw.flat_peaks {
            return Ok(to_koth_peaks(raw_peaks.iter().map(|p| (p.mz, p.abundance))));
        }
        // No stored centroid list: fall back to picking the profile signal.
        // Its m/z calibration comes from the scan event, so with misaligned
        // events the scan is skipped rather than mis-calibrated.
        if !self.events_aligned {
            log::warn!(
                "Thermo scan {scan_number}: no centroid list and no trusted profile \
                 calibration; skipping"
            );
            return Ok(Vec::new());
        }
        let packet = self
            .raw
            .read_scan(&mut self.file, scan_number)
            .map_err(read_err)?;
        let Some(profile) = packet.profile else {
            return Ok(Vec::new());
        };
        let coefficients = self
            .raw
            .scan_events
            .get(idx as usize)
            .map(|e| e.coefficients.as_slice())
            .unwrap_or(&[]);
        let (mz, intensity): (Vec<f64>, Vec<f32>) = profile
            .to_mz_intensity(coefficients)
            .into_iter()
            .filter(|&(m, _)| m > 0.0)
            .map(|(m, i)| (m, i as f32))
            .unzip();
        Ok(pick_profile_peaks(&mz, &intensity))
    }

    /// Isolation window of MS2 scan `idx`, or `None` if the scan carries no
    /// usable precursor. The center is the scan event's reaction precursor m/z
    /// (the value in the scan filter, e.g. `ms2 512.50@hcd25`), which for DIA is
    /// the window center; the width is the trailer's MS2 isolation width.
    pub(crate) fn isolation_window(&self, idx: u32) -> Option<IsolationWindow> {
        let reaction = self
            .raw
            .scan_events
            .get(idx as usize)
            .filter(|_| self.events_aligned)
            .and_then(|e| e.reactions.first());
        let params = self.raw.scan_params(self.scan_number(idx));
        let fallback = params
            .as_ref()
            .and_then(|p| p.monoisotopic_mz())
            .unwrap_or(f64::NAN);
        let center = reaction
            .map(|r| r.precursor_mz)
            .filter(|mz| mz.is_finite() && *mz > 0.0)
            .unwrap_or(fallback);
        let width = params
            .as_ref()
            .and_then(|p| p.isolation_width_mz())
            .unwrap_or(f64::NAN);
        let (lower, upper) = if width.is_finite() && width > 0.0 {
            (center - width / 2.0, center + width / 2.0)
        } else {
            (f64::NAN, f64::NAN)
        };
        derive_isolation_window(center, lower, upper, fallback)
    }
}

/// Parse a Thermo `FAIMS CV` trailer value, stored numerically on most
/// firmware and as text (optionally with a unit) on some.
fn faims_value(value: &GenericValue) -> Option<f32> {
    let cv = match value {
        GenericValue::Float32(v) => *v,
        GenericValue::Float64(v) => *v as f32,
        GenericValue::Int8(v) => f32::from(*v),
        GenericValue::Int16(v) => f32::from(*v),
        GenericValue::Int32(v) => *v as f32,
        GenericValue::String(s) => return parse_faims_cv(s),
        _ => return None,
    };
    cv.is_finite().then_some(if cv == 0.0 { 0.0 } else { cv })
}

/// Parse a textual `FAIMS CV` trailer value. Accepting a trailing unit is
/// harmless and robust.
fn parse_faims_cv(value: &str) -> Option<f32> {
    let token = value.split_whitespace().next()?;
    let cv = token.parse::<f32>().ok()?;
    cv.is_finite().then_some(if cv == 0.0 { 0.0 } else { cv })
}

/// Convert (m/z, intensity) pairs to koth [`Peak`]s: filter to positive
/// intensity, m/z into `f32`, ion mobility 0.0 (Orbitrap has no IM), sorted by
/// m/z ascending. The single scan→peaks conversion shared by the MS1 and MS2
/// paths, matching `mzml::extract_peaks`'s centroided branch.
fn to_koth_peaks(pairs: impl Iterator<Item = (f64, f32)>) -> Vec<Peak> {
    let mut peaks: Vec<Peak> = pairs
        .filter(|&(_, it)| it > 0.0)
        .map(|(m, it)| Peak {
            mz: m as f32,
            intensity: it,
            ion_mobility: 0.0,
        })
        .collect();
    peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));
    peaks
}

/// Centroid a profile signal with mzdata's quadratic peak picker at S/N >= 1,
/// the same picker and threshold the mzML path applies to profile-mode mzML.
fn pick_profile_peaks(mz: &[f64], intensity: &[f32]) -> Vec<Peak> {
    use mzdata::mzsignal::{PeakFitType, PeakPicker};
    if mz.len() < 3 {
        return Vec::new();
    }
    let picker = PeakPicker {
        fit_type: PeakFitType::Quadratic,
        signal_to_noise_threshold: 1.0,
        ..Default::default()
    };
    if mz.windows(2).any(|w| w[1] < w[0]) {
        log::warn!("Thermo profile scan m/z is not sorted; skipping peak picking");
        return Vec::new();
    }
    let mut fitted = Vec::new();
    if let Err(e) = picker.discover_peaks(mz, intensity, &mut fitted) {
        log::warn!("Thermo profile peak picking failed: {e}");
        return Vec::new();
    }
    to_koth_peaks(fitted.iter().map(|p| (p.mz, p.intensity)))
}

/// Read all MS1 (centroided) spectra from a Thermo `.raw` file, returning the same
/// `Vec<Spectrum>` shape as [`crate::io::mzml::read_mzml`].
pub fn read_thermo(path: &Path) -> Result<Vec<Spectrum>, KothError> {
    stream_thermo(path)?.collect()
}

/// Index scan metadata without decoding signal, then decode MS1 spectra
/// lazily in stable retention-time order. Errors stop the stream; empty spectra
/// are skipped and emitted spectra receive contiguous scan indices.
pub fn stream_thermo(path: &Path) -> Result<super::SpectrumStream, KothError> {
    let mut source = RawSource::open(path)?;
    let mut indices = Vec::new();
    for i in 0..source.n_scans() {
        if source.ms_level(i) == 1 {
            let rt = source.retention_time(i);
            if !rt.is_finite() {
                return Err(KothError::ThermoError(format!(
                    "non-finite RT for scan {i}"
                )));
            }
            indices.push((i, rt));
        }
    }
    // Stable sort: acquisition order is kept for equal RTs.
    indices.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
    let mut indices = indices.into_iter();
    let mut count = 0;
    let mut finished = false;
    Ok(Box::new(std::iter::from_fn(move || {
        if finished {
            return None;
        }
        for (i, retention_time) in indices.by_ref() {
            let peaks = match source.peaks(i) {
                Ok(p) => p,
                Err(e) => {
                    finished = true;
                    return Some(Err(e));
                }
            };
            if peaks.is_empty() {
                continue;
            }
            let spectrum = Spectrum {
                scan_index: count,
                retention_time,
                peaks,
                ms_level: 1,
                isolation_window: None,
                faims_cv: source.faims_cv(i),
            };
            count += 1;
            return Some(Ok(spectrum));
        }
        finished = true;
        if count == 0 {
            Some(Err(KothError::NoSpectra))
        } else {
            None
        }
    })))
}

// ---------------------------------------------------------------------------
// DIA MS2 reader
// ---------------------------------------------------------------------------

/// Fewer MS2 scans than this and we refuse to classify DIA-vs-DDA (too little
/// evidence for the recurrence test) and emit no MS2. Any real acquisition has
/// far more; this only guards pathological / truncated files.
const MIN_MS2_FOR_CLASSIFICATION: usize = 8;

/// Upper bound on the number of *distinct* isolation windows a DIA schedule may
/// have. Generous — classic DIA has 8–70 windows, but staggered / gas-phase-
/// fractionation methods can reach several hundred. Above this we assume the
/// per-scan windows are data-dependent (DDA) rather than a fixed schedule.
const MAX_DIA_WINDOWS: usize = 2_000;

/// Minimum mean recurrence (`n_ms2 / distinct_windows`) for a DIA verdict. A
/// fixed DIA schedule revisits each window every cycle, so recurrence is large
/// (tens to thousands). DDA with dynamic exclusion selects each precursor ~once,
/// so recurrence ≈ 1. `3.0` cleanly separates the two while tolerating a short
/// run of only a few cycles.
const MIN_RECURRENCE: f64 = 3.0;

/// Summary statistics of the per-MS2-scan isolation windows, used both to make
/// the DIA verdict and to report it in logs.
#[derive(Debug, Clone, Copy)]
struct ScheduleStats {
    n_ms2: usize,
    distinct: usize,
    recurrence: f64,
}

fn schedule_stats(windows: &[IsolationWindow]) -> ScheduleStats {
    let n_ms2 = windows.len();
    let distinct = windows
        .iter()
        .map(|w| w.key())
        .collect::<BTreeSet<_>>()
        .len();
    let recurrence = if distinct == 0 {
        0.0
    } else {
        n_ms2 as f64 / distinct as f64
    };
    ScheduleStats {
        n_ms2,
        distinct,
        recurrence,
    }
}

/// Heuristic DIA detector over the per-MS2-scan isolation windows of a `.raw`.
///
/// **DIA** acquires a small, fixed set of isolation windows that repeat on a
/// schedule, so a handful of distinct windows each recur many times. **DDA**
/// selects data-dependent precursors, so (with dynamic exclusion) nearly every
/// MS2 scan has a distinct, narrow precursor window that appears roughly once.
/// The discriminating signal is therefore the **mean recurrence** of the windows
/// (`n_ms2 / distinct`): large for DIA, ≈1 for DDA. We require a small-enough
/// distinct-window count and a recurrence at or above [`MIN_RECURRENCE`], with a
/// floor of [`MIN_MS2_FOR_CLASSIFICATION`] scans so the ratio is meaningful. When
/// unsure we return `false` (treat as DDA → emit nothing), which is the safe
/// direction: it never fabricates DIA channels from DDA precursors.
///
/// **Known failure modes** (honest):
/// * **Targeted PRM / tMS2** uses a fixed, repeating inclusion list of narrow
///   windows and *will* be classified as DIA. Feeding its per-window MS2 to the
///   tracer is not DDA precursor reconstruction, but it is not a true DIA tiling
///   either; the caller should not point PRM data at the MS2 surface.
/// * **All-ion / MSᴱ / bbCID** (one very wide window, or a couple alternating)
///   is reported as DIA with 1–2 channels. That is arguably correct (it *is*
///   windowed fragmentation), but the single wide channel is coarse.
/// * **Highly multiplexed / scanning-quadrupole DIA** whose per-scan window
///   centers are (nearly) all distinct would be *missed* (recurrence ≈ 1 →
///   classified DDA → no MS2). Not a Thermo default; noted for completeness.
/// * A **DDA** run with unusually aggressive re-selection (recurrence ≥ 3) would
///   be misread as DIA — not seen in practice with dynamic exclusion on.
fn is_dia_schedule(windows: &[IsolationWindow]) -> bool {
    let s = schedule_stats(windows);
    s.n_ms2 >= MIN_MS2_FOR_CLASSIFICATION
        && s.distinct >= 1
        && s.distinct <= MAX_DIA_WINDOWS
        && s.recurrence >= MIN_RECURRENCE
}

/// Map a Thermo scan's precursor/isolation metadata to koth's [`IsolationWindow`]
/// (absolute m/z bounds). Prefers the scan event's isolation `target` as the
/// center, falling back to the trailer precursor m/z; if `lower`/`upper` do
/// not bracket a positive-width window (unset / degenerate on some files) the
/// window collapses to the center point. Returns `None` when no usable center
/// exists. Pure and unit-testable.
fn derive_isolation_window(
    target: f64,
    lower: f64,
    upper: f64,
    precursor_mz: f64,
) -> Option<IsolationWindow> {
    let center = if target.is_finite() && target > 0.0 {
        target
    } else if precursor_mz.is_finite() && precursor_mz > 0.0 {
        precursor_mz
    } else {
        return None;
    };
    let (lo, hi) = if lower.is_finite() && upper.is_finite() && upper > lower && lower > 0.0 {
        (lower, upper)
    } else {
        (center, center)
    };
    Some(IsolationWindow {
        target: center,
        lower: lo,
        upper: hi,
    })
}

/// Read all **DIA** MS2 scans from a Thermo `.raw`, one [`Spectrum`] per MS2 scan
/// stamped with its precursor isolation window — the `.raw` analog of
/// [`crate::io::mzml::stream_mzml_ms2`]. Orbitrap and Astral scans are 1-D
/// centroids (no ion mobility), so each MS2 scan maps to exactly one Spectrum
/// with no segmentation.
///
/// **DIA-only.** A first, signal-free pass collects the MS2 isolation windows and
/// classifies the acquisition via `is_dia_schedule`. If it is not confidently
/// DIA (DDA, PRM notwithstanding, or too few scans), this logs a warning and
/// returns an empty Vec — it never reconstructs DDA precursors and never touches
/// MS1. Only on a DIA verdict does a second pass decode the MS2 peak arrays.
pub fn read_thermo_ms2(path: &Path) -> Result<Vec<Spectrum>, KothError> {
    let mut source = RawSource::open(path)?;

    // Pass 1 (cheap): classify DIA vs DDA from MS2 isolation windows alone;
    // no peak arrays are decoded for a file we may reject.
    let mut ms2_scans: Vec<(u32, IsolationWindow)> = Vec::new();
    for i in 0..source.n_scans() {
        if source.ms_level(i) != 2 {
            continue;
        }
        if let Some(w) = source.isolation_window(i) {
            ms2_scans.push((i, w));
        }
    }

    if ms2_scans.is_empty() {
        log::warn!(
            "Thermo .raw '{}' has no MS2 scans with an isolation window; emitting no MS2 \
             (MS1 output is unaffected)",
            path.display()
        );
        return Ok(Vec::new());
    }

    let scan_windows: Vec<IsolationWindow> = ms2_scans.iter().map(|&(_, w)| w).collect();
    let stats = schedule_stats(&scan_windows);
    if !is_dia_schedule(&scan_windows) {
        log::warn!(
            "Thermo .raw '{}' does not look like DIA ({} MS2 scans, {} distinct isolation \
             windows, {:.1}x mean recurrence) — treating as DDA and emitting no MS2 \
             (DDA precursor reconstruction is out of scope; MS1 output is unaffected)",
            path.display(),
            stats.n_ms2,
            stats.distinct,
            stats.recurrence,
        );
        return Ok(Vec::new());
    }

    // Pass 2: DIA confirmed — decode centroided MS2 peaks and build one Spectrum
    // per scan, reusing the exact scan→peaks conversion the MS1 path uses.
    let mut spectra: Vec<Spectrum> = Vec::with_capacity(ms2_scans.len());
    for (i, window) in ms2_scans {
        let peaks = source.peaks(i)?;
        if peaks.is_empty() {
            continue;
        }
        spectra.push(Spectrum {
            // Source scan index, for traceability. The MS2 hill detector
            // re-indexes scans per isolation window, so this is not used for
            // gap tracking (mirrors the Bruker diaPASEF reader).
            scan_index: i as usize,
            retention_time: source.retention_time(i), // minutes
            peaks,
            ms_level: 2,
            isolation_window: Some(window),
            faims_cv: None,
        });
    }

    // Ascending RT so each per-window detector sees its cycles in acquisition
    // order (scans are already in acquisition order; this is a safety net).
    spectra.sort_by(|a, b| {
        a.retention_time
            .partial_cmp(&b.retention_time)
            .unwrap_or(Ordering::Equal)
    });

    let n_windows = spectra
        .iter()
        .filter_map(|s| s.isolation_window.map(|w| w.key()))
        .collect::<BTreeSet<_>>()
        .len();
    log::info!(
        "Read {} MS2 spectra across {} isolation windows from Thermo .raw \
         (DIA, {:.1}x mean recurrence)",
        spectra.len(),
        n_windows,
        stats.recurrence,
    );
    Ok(spectra)
}
