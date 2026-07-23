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

/// Detect input format and read all MS1 spectra.
///
/// Supported formats:
/// - `.mzML` / `.mzml` / `.mzML.gz` → mzML via mzdata
/// - directory ending in `.d` → Bruker timsTOF TDF (requires "tdf" feature)
/// - `.raw` → native Thermo Fisher reader (requires "thermo" feature + .NET runtime)
pub fn read_spectra(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .to_lowercase();

    if name.ends_with(".mzml") || name.ends_with(".mzml.gz") {
        mzml::read_mzml(path)
    } else if path.is_dir() && name.ends_with(".d") {
        read_bruker_inner(path, file)
    } else if name.ends_with(".raw") {
        read_thermo_inner(path)
    } else {
        Err(KothError::UnsupportedFormat(format!(
            "'{}' is not a supported input format. Use .mzML, a Bruker .d directory, or a Thermo .raw file.",
            path.display()
        )))
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
