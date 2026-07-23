//! Native Thermo Fisher `.raw` reader (opt-in `thermo` cargo feature).
//!
//! Wraps the [`thermorawfilereader`] crate — a self-hosted .NET runtime over
//! Thermo's `RawFileReader` assemblies — so a `.raw` file can be fed to koth_ff
//! directly without a prior mzML conversion. The .NET 8 runtime must be installed
//! on the host for reading to succeed at run time.
//!
//! The output matches [`crate::io::mzml::read_mzml`] exactly: **MS1 spectra only**,
//! **centroided** (the reader is asked for centroids, since koth's hill detector
//! consumes centroids and does no peak picking of its own), peaks filtered to
//! positive intensity and sorted by m/z ascending, spectra sorted by retention
//! time and re-indexed `0..n`. Ion mobility is always 0.0 (Orbitrap has no IM
//! dimension). Reading is **local-file only** — the Thermo API takes a filesystem
//! path and manipulates file locks, so it cannot consume streams or remote stores.

use std::cmp::Ordering;
use std::path::Path;

use thermorawfilereader::RawFileReader;

use crate::error::KothError;
use crate::models::{Peak, Spectrum};

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
        KothError::TdfError(format!(
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

        // Centroid peaks: filter to positive intensity and sort by m/z ascending,
        // mirroring `mzml::extract_peaks`.
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
        if peaks.is_empty() {
            continue;
        }
        peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));

        spectra.push(Spectrum {
            scan_index: 0, // assigned after the RT sort below
            retention_time: spec.time(), // Thermo reports scan time in minutes
            peaks,
            ms_level: 1,
            isolation_window: None,
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
