pub mod mzml;
pub mod prefetch;

#[cfg(feature = "tdf")]
pub mod bruker;

#[cfg(feature = "thermo")]
pub mod thermo;

use std::path::Path;

use crate::config::FileConfig;
use crate::error::KothError;
use crate::models::Spectrum;

/// Classification of an input path into one of koth_ff's supported raw formats.
///
/// The variant reflects the *shape* of the input, independent of which cargo
/// features are compiled in — `BrukerD`/`ThermoRaw` are returned even in a build
/// without `tdf`/`thermo`, and the caller is responsible for surfacing the
/// "rebuild with --features …" error in that case.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    /// Plain mzML (`.mzml`, case-insensitive).
    Mzml,
    /// Gzip-compressed mzML (`.mzml.gz`, case-insensitive).
    MzmlGz,
    /// Bruker timsTOF `.d` folder (requires the `tdf` feature to read).
    BrukerD,
    /// Thermo Fisher `.raw` file (requires the `thermo` feature to read).
    ThermoRaw,
    /// Nothing koth_ff can read.
    Unsupported,
}

/// Detect a path's input format from its name (and, for Bruker, directory-ness).
///
/// **Canonical predicates** — this is the single source of truth for the
/// `.mzML`/`.d`/`.raw` decision; every dispatch site routes through it so the
/// edge cases stay consistent. All extension tests are **case-insensitive**:
///
/// * `MzmlGz`   — file name ends `.mzml.gz`.
/// * `Mzml`     — file name ends `.mzml` (and not `.mzml.gz`).
/// * `BrukerD`  — the path **is a directory** AND its name ends `.d`. A `.d`
///   name that is not a directory, or any directory whose name does not end
///   `.d`, is **not** Bruker.
/// * `ThermoRaw`— file name ends `.raw`.
/// * `Unsupported` — anything else.
pub fn detect_format(path: &Path) -> InputFormat {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .to_lowercase();

    if name.ends_with(".mzml.gz") {
        InputFormat::MzmlGz
    } else if name.ends_with(".mzml") {
        InputFormat::Mzml
    } else if name.ends_with(".d") && path.is_dir() {
        InputFormat::BrukerD
    } else if name.ends_with(".raw") {
        InputFormat::ThermoRaw
    } else {
        InputFormat::Unsupported
    }
}

/// Detect input format and read all MS1 spectra.
///
/// Supported formats:
/// - `.mzML` / `.mzml` / `.mzML.gz` → mzML via mzdata
/// - directory ending in `.d` → Bruker timsTOF TDF (requires "tdf" feature)
/// - `.raw` → native Thermo Fisher reader (requires "thermo" feature + .NET runtime)
pub fn read_spectra(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
    match detect_format(path) {
        InputFormat::Mzml | InputFormat::MzmlGz => mzml::read_mzml(path),
        InputFormat::BrukerD => read_bruker_inner(path, file),
        InputFormat::ThermoRaw => read_thermo_inner(path),
        InputFormat::Unsupported => Err(KothError::UnsupportedFormat(format!(
            "'{}' is not a supported input format. Use .mzML, a Bruker .d directory, or a Thermo .raw file.",
            path.display()
        ))),
    }
}

#[cfg(feature = "thermo")]
fn read_thermo_inner(path: &Path) -> Result<Vec<Spectrum>, KothError> {
    thermo::read_thermo(path)
}

#[cfg(not(feature = "thermo"))]
fn read_thermo_inner(_path: &Path) -> Result<Vec<Spectrum>, KothError> {
    Err(KothError::UnsupportedFormat(
        "Thermo .raw support requires the 'thermo' feature flag. Rebuild with --features thermo.".into(),
    ))
}

#[cfg(feature = "tdf")]
fn read_bruker_inner(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
    bruker::read_bruker(path, file)
}

#[cfg(not(feature = "tdf"))]
fn read_bruker_inner(_path: &Path, _file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
    Err(KothError::UnsupportedFormat(
        "Bruker .d support requires the 'tdf' feature flag. Rebuild with --features tdf.".into(),
    ))
}
