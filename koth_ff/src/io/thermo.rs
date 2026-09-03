//! Native Thermo Fisher `.raw` reader (opt-in `thermo` cargo feature).
//!
//! Wraps the [`thermorawfilereader`] crate — a self-hosted .NET runtime over
//! Thermo's `RawFileReader` assemblies — so a `.raw` file can be fed to koth_ff
//! directly without a prior mzML conversion. The .NET 8 runtime must be installed
//! on the host for reading to succeed at run time.
//!
//! The MS1 output matches [`crate::io::mzml::read_mzml`] exactly: **MS1 spectra
//! only**, **centroided** (the reader is asked for centroids, since koth's hill
//! detector consumes centroids and does no peak picking of its own), peaks
//! filtered to positive intensity and sorted by m/z ascending, spectra sorted by
//! retention time and re-indexed `0..n`. Ion mobility is always 0.0 (Orbitrap has
//! no IM dimension). Reading is **local-file only** — the Thermo API takes a
//! filesystem path and manipulates file locks, so it cannot consume streams or
//! remote stores.
//!
//! [`read_thermo_ms2`] is the opt-in **DIA MS2** counterpart: it emits one MS2
//! [`Spectrum`] per MS2 scan, stamped with the scan's precursor isolation window,
//! mirroring the mzML MS2 path ([`crate::io::mzml::stream_mzml_ms2`]). It is
//! **DIA-only** — a DDA `.raw` (or one it cannot confidently classify as DIA)
//! yields an empty set plus a warning and never disturbs MS1. See
//! [`is_dia_schedule`] for the heuristic and its failure modes.

use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::path::Path;

use thermorawfilereader::{RawFileReader, RawSpectrum};

use crate::error::KothError;
use crate::models::{IsolationWindow, Peak, Spectrum};

#[cfg(test)]
#[path = "thermo_tests.rs"]
mod tests;

/// Best-effort discovery of the .NET runtime for the `nethost`/`hostfxr` layer.
///
/// The Thermo reader hosts a .NET 8 runtime, which `nethost` locates via
/// `DOTNET_ROOT` or a small set of standard install paths. A user-local install
/// (the `dotnet-install.sh` default `~/.dotnet`) is *not* on that list, so if
/// `DOTNET_ROOT` is unset we probe the common locations and point it at the first
/// that actually carries a `Microsoft.NETCore.App` runtime. An explicit
/// `DOTNET_ROOT` is always respected; if nothing is found we leave it unset and
/// let [`RawFileReader::open`] surface the "missing .NET runtime" error.
fn ensure_dotnet_root() {
    if std::env::var_os("DOTNET_ROOT").is_some() {
        return;
    }
    let mut candidates: Vec<std::path::PathBuf> = Vec::new();
    if let Some(home) = std::env::var_os("HOME") {
        candidates.push(Path::new(&home).join(".dotnet"));
    }
    for p in ["/usr/share/dotnet", "/usr/lib/dotnet", "/opt/dotnet"] {
        candidates.push(Path::new(p).to_path_buf());
    }
    for root in candidates {
        if root.join("shared/Microsoft.NETCore.App").is_dir() {
            log::info!(
                "Thermo reader: auto-detected .NET runtime, setting DOTNET_ROOT={}",
                root.display()
            );
            std::env::set_var("DOTNET_ROOT", &root);
            return;
        }
    }
}

/// Read all MS1 (centroided) spectra from a Thermo `.raw` file, returning the same
/// `Vec<Spectrum>` shape as [`crate::io::mzml::read_mzml`].
pub fn read_thermo(path: &Path) -> Result<Vec<Spectrum>, KothError> {
    ensure_dotnet_root();
    let mut reader = RawFileReader::open(path).map_err(|e| {
        KothError::ThermoError(format!(
            "could not open Thermo .raw file '{}': {e} \
             (a .NET 8 runtime must be installed for native .raw reading)",
            path.display()
        ))
    })?;
    // koth's hill detector works on centroids and does no peak picking; ask the
    // Thermo reader to centroid so profile-mode MS1 scans arrive as peak lists.
    reader.set_centroid_spectra(true);

    let n = reader.len();
    let mut spectra: Vec<Spectrum> = Vec::with_capacity(n);

    for i in 0..n {
        let Some(spec) = reader.get(i) else { continue };
        if spec.ms_level() != 1 {
            continue; // MS1 only, matching read_mzml
        }

        let peaks = centroid_peaks(&spec);
        if peaks.is_empty() {
            continue;
        }
        let faims_cv = reader.get_raw_trailers_for(i).and_then(|trailers| {
            trailers
                .get_label("FAIMS CV")
                .and_then(|v| parse_faims_cv(v.value))
        });

        spectra.push(Spectrum {
            scan_index: 0,               // assigned after the RT sort below
            retention_time: spec.time(), // Thermo reports scan time in minutes
            peaks,
            ms_level: 1,
            isolation_window: None,
            faims_cv,
        });
    }

    if spectra.is_empty() {
        return Err(KothError::NoSpectra);
    }

    // Sort by retention time and re-index, identical to `read_mzml`.
    spectra.sort_by(|a, b| {
        a.retention_time
            .partial_cmp(&b.retention_time)
            .unwrap_or(Ordering::Equal)
    });
    for (i, s) in spectra.iter_mut().enumerate() {
        s.scan_index = i;
    }

    log::info!("Read {} MS1 spectra from Thermo .raw", spectra.len());
    Ok(spectra)
}

/// Parse the Thermo `FAIMS CV` trailer value. Trailer values normally contain
/// just the number, but accepting a trailing unit is harmless and robust.
fn parse_faims_cv(value: &str) -> Option<f32> {
    let token = value.split_whitespace().next()?;
    let cv = token.parse::<f32>().ok()?;
    cv.is_finite().then_some(if cv == 0.0 { 0.0 } else { cv })
}

/// Convert a spectrum's centroid arrays to koth [`Peak`]s: filter to positive
/// intensity, m/z into `f32`, ion mobility 0.0 (Orbitrap has no IM), sorted by
/// m/z ascending. This is the **single** scan→peaks conversion shared by the MS1
/// ([`read_thermo`]) and MS2 ([`read_thermo_ms2`]) paths, matching
/// `mzml::extract_peaks`'s centroided branch. The reader must have
/// `set_centroid_spectra(true)` so profile scans arrive as peak lists.
fn centroid_peaks(spec: &RawSpectrum) -> Vec<Peak> {
    let mut peaks: Vec<Peak> = match spec.data() {
        Some(d) => {
            let mz = d.mz();
            let intensity = d.intensity();
            mz.iter()
                .zip(intensity.iter())
                .filter(|(_, &it)| it > 0.0)
                .map(|(&m, &it)| Peak {
                    mz: m as f32,
                    intensity: it,
                    ion_mobility: 0.0,
                })
                .collect()
        }
        None => Vec::new(),
    };
    peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));
    peaks
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
/// (absolute m/z bounds). Prefers the recorded isolation-window `target` as the
/// center, falling back to the precursor m/z; if the recorded `lower`/`upper` do
/// not bracket a positive-width window (unset / degenerate on some files) the
/// window collapses to the center point. Returns `None` when no usable center
/// exists. Pure and unit-testable — the FFI bridge lives in [`precursor_window`].
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

/// Bridge a [`RawSpectrum`]'s precursor to an [`IsolationWindow`], or `None` if
/// the scan carries no precursor (e.g. an MS1 scan, or a malformed MS2 scan).
fn precursor_window(spec: &RawSpectrum) -> Option<IsolationWindow> {
    let p = spec.precursor()?;
    let iso = p.isolation_window();
    derive_isolation_window(iso.target(), iso.lower(), iso.upper(), p.mz())
}

/// Read all **DIA** MS2 scans from a Thermo `.raw`, one [`Spectrum`] per MS2 scan
/// stamped with its precursor isolation window — the `.raw` analog of
/// [`crate::io::mzml::stream_mzml_ms2`]. Orbitrap scans are 1-D centroids (no ion
/// mobility), so each MS2 scan maps to exactly one Spectrum with no segmentation.
///
/// **DIA-only.** A first, signal-free pass collects the MS2 isolation windows and
/// classifies the acquisition via [`is_dia_schedule`]. If it is not confidently
/// DIA (DDA, PRM notwithstanding, or too few scans), this logs a warning and
/// returns an empty Vec — it never reconstructs DDA precursors and never touches
/// MS1. Only on a DIA verdict does a second pass decode the MS2 peak arrays.
pub fn read_thermo_ms2(path: &Path) -> Result<Vec<Spectrum>, KothError> {
    ensure_dotnet_root();
    let mut reader = RawFileReader::open(path).map_err(|e| {
        KothError::ThermoError(format!(
            "could not open Thermo .raw file '{}': {e} \
             (a .NET 8 runtime must be installed for native .raw reading)",
            path.display()
        ))
    })?;

    // Pass 1 (cheap): classify DIA vs DDA from MS2 isolation windows alone, with
    // signal loading OFF so no peak arrays are decoded for a file we may reject.
    reader.set_signal_loading(false);
    let n = reader.len();
    let mut scan_windows: Vec<IsolationWindow> = Vec::new();
    for i in 0..n {
        let Some(spec) = reader.get(i) else { continue };
        if spec.ms_level() != 2 {
            continue;
        }
        if let Some(w) = precursor_window(&spec) {
            scan_windows.push(w);
        }
    }

    if scan_windows.is_empty() {
        log::warn!(
            "Thermo .raw '{}' has no MS2 scans with an isolation window; emitting no MS2 \
             (MS1 output is unaffected)",
            path.display()
        );
        return Ok(Vec::new());
    }

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
    reader.set_signal_loading(true);
    reader.set_centroid_spectra(true);
    let mut spectra: Vec<Spectrum> = Vec::with_capacity(scan_windows.len());
    for i in 0..n {
        let Some(spec) = reader.get(i) else { continue };
        if spec.ms_level() != 2 {
            continue;
        }
        let Some(window) = precursor_window(&spec) else {
            continue;
        };
        let peaks = centroid_peaks(&spec);
        if peaks.is_empty() {
            continue;
        }
        spectra.push(Spectrum {
            // Source scan index, for traceability. The MS2 hill detector
            // re-indexes scans per isolation window, so this is not used for
            // gap tracking (mirrors the Bruker diaPASEF reader).
            scan_index: i,
            retention_time: spec.time(), // Thermo reports scan time in minutes
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
