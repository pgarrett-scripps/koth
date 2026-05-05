pub mod mzml;

#[cfg(feature = "tdf")]
pub mod bruker;

use std::path::Path;

use crate::config::FileConfig;
use crate::error::KothError;
use crate::models::Spectrum;

/// Detect input format and read all MS1 spectra.
///
/// Supported formats:
/// - `.mzML` / `.mzml` / `.mzML.gz` → mzML via mzdata
/// - directory ending in `.d` → Bruker timsTOF TDF (requires "tdf" feature)
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
    } else {
        Err(KothError::UnsupportedFormat(format!(
            "'{}' is not a supported input format. Use .mzML or a Bruker .d directory.",
            path.display()
        )))
    }
}

#[cfg(feature = "tdf")]
fn read_bruker_inner(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
    bruker::read_bruker(path, file.bruker_mz_ppm, file.bruker_im_pct, file.bruker_min_subpeaks)
}

#[cfg(not(feature = "tdf"))]
fn read_bruker_inner(_path: &Path, _file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
    Err(KothError::UnsupportedFormat(
        "Bruker .d support requires the 'tdf' feature flag. Rebuild with --features tdf.".into(),
    ))
}
