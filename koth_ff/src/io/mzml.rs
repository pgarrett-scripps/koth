use std::path::Path;
use std::{fs, io};

use flate2::read::GzDecoder;
use mzdata::prelude::*;
use mzdata::spectrum::{MultiLayerSpectrum, SignalContinuity};
use mzdata::MZReader;

use crate::error::KothError;
use crate::models::{IsolationWindow, Peak, Spectrum};

/// Read an mzML file (plain or gzip-compressed) and return all MS1 spectra
/// sorted by retention time.
pub fn read_mzml(path: &Path) -> Result<Vec<Spectrum>, KothError> {
    let mut spectra: Vec<Spectrum> = collect_ms1(open_reader(path)?);

    if spectra.is_empty() {
        return Err(KothError::NoSpectra);
    }

    spectra.sort_by(|a, b| {
        a.retention_time
            .partial_cmp(&b.retention_time)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for (i, s) in spectra.iter_mut().enumerate() {
        s.scan_index = i;
    }

    log::info!("Read {} MS1 spectra from mzML", spectra.len());
    Ok(spectra)
}

/// Stream MS1 spectra from an mzML file (plain or gzip-compressed) one at a time.
///
/// Spectra are yielded in file order (assumed to be RT order for standard
/// LC-MS acquisitions). Each spectrum is dropped after the caller processes it.
pub fn stream_mzml(path: &Path) -> Result<Box<dyn Iterator<Item = Spectrum> + Send>, KothError> {
    crate::mem::log_mem("before open_reader (stream_mzml)");
    Ok(Box::new(ms1_stream(open_reader(path)?)))
}

/// Stream MS2 spectra (with isolation window metadata) from an mzML file.
///
/// Only spectra with `ms_level == 2` and a parsable precursor isolation window
/// are yielded. The `scan_index` on each yielded spectrum is its absolute
/// position in the file, so callers can later re-index per-isolation-window.
pub fn stream_mzml_ms2(
    path: &Path,
) -> Result<Box<dyn Iterator<Item = Spectrum> + Send>, KothError> {
    crate::mem::log_mem("before open_reader (stream_mzml_ms2)");
    Ok(Box::new(ms2_stream(open_reader(path)?)))
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Boxed iterator over MultiLayerSpectrum that works for both compressed and
/// uncompressed readers.
type BoxedRawIter = Box<dyn Iterator<Item = MultiLayerSpectrum> + Send>;

fn open_reader(path: &Path) -> Result<BoxedRawIter, KothError> {
    // Route the compressed-vs-plain decision through the single canonical format
    // predicate. `open_reader` is only reached for paths already classified as
    // mzML, so anything not `MzmlGz` here is plain (uncompressed) mzML.
    if super::detect_format(path) == super::InputFormat::MzmlGz {
        // Decompress the full file into memory so mzdata gets a seekable Cursor.
        // open_gzipped_read uses a limited PreBufferedStream that corrupts large
        // base64 arrays (e.g. zlib-compressed Orbitrap data) during streaming.
        let file = fs::File::open(path).map_err(|e| KothError::MzmlError(e.to_string()))?;
        let mut decoder = GzDecoder::new(io::BufReader::new(file));
        let mut buf = Vec::new();
        io::Read::read_to_end(&mut decoder, &mut buf)
            .map_err(|e| KothError::MzmlError(format!("gzip decompress: {e}")))?;
        let cursor = io::Cursor::new(buf);
        let reader =
            MZReader::open_read_seek(cursor).map_err(|e| KothError::MzmlError(e.to_string()))?;
        Ok(Box::new(reader))
    } else {
        let reader = MZReader::open_path(path).map_err(|e| KothError::MzmlError(e.to_string()))?;
        Ok(Box::new(reader))
    }
}

fn collect_ms1(reader: BoxedRawIter) -> Vec<Spectrum> {
    let mut out = Vec::new();
    for (scan_index, mut spectrum) in reader.enumerate() {
        if spectrum.ms_level() != 1 {
            continue;
        }
        let faims_cv = extract_faims_cv(&spectrum);
        let peaks = extract_peaks(&mut spectrum);
        if peaks.is_empty() {
            continue;
        }
        out.push(Spectrum {
            scan_index,
            retention_time: spectrum.start_time(),
            peaks,
            ms_level: 1,
            isolation_window: None,
            faims_cv,
        });
    }
    out
}

fn ms1_stream(reader: BoxedRawIter) -> impl Iterator<Item = Spectrum> {
    let mut scan_index = 0usize;
    let mut ms1_count = 0usize;

    reader.filter_map(move |mut spectrum| {
        let idx = scan_index;
        scan_index += 1;

        if spectrum.ms_level() != 1 {
            return None;
        }

        if ms1_count == 0 {
            let continuity = spectrum.signal_continuity();
            log::info!("First MS1 scan signal continuity: {:?}", continuity);
            if matches!(continuity, SignalContinuity::Profile) {
                log::info!(
                    "Profile-mode data detected — centroiding on the fly with mzdata peak picker."
                );
            }
        }

        let retention_time = spectrum.start_time();
        let faims_cv = extract_faims_cv(&spectrum);
        let peaks = extract_peaks(&mut spectrum);

        if ms1_count == 0 {
            log::info!(
                "First MS1 scan: scan_index={} rt={:.2} min  peaks={}",
                idx,
                retention_time,
                peaks.len()
            );
        }

        ms1_count += 1;
        if peaks.is_empty() {
            return None;
        }
        Some(Spectrum {
            scan_index: idx,
            retention_time,
            peaks,
            ms_level: 1,
            isolation_window: None,
            faims_cv,
        })
    })
}

fn ms2_stream(reader: BoxedRawIter) -> impl Iterator<Item = Spectrum> {
    let mut scan_index = 0usize;
    let mut ms2_count = 0usize;
    let mut warned_missing_precursor = false;

    reader.filter_map(move |mut spectrum| {
        let idx = scan_index;
        scan_index += 1;

        if spectrum.ms_level() != 2 {
            return None;
        }

        let iw = match spectrum.precursor() {
            Some(p) => {
                let raw = p.isolation_window();
                let target = raw.target as f64;
                let lower = match raw.flags {
                    mzdata::spectrum::IsolationWindowState::Offset => {
                        target - raw.lower_bound as f64
                    }
                    _ => raw.lower_bound as f64,
                };
                let upper = match raw.flags {
                    mzdata::spectrum::IsolationWindowState::Offset => {
                        target + raw.upper_bound as f64
                    }
                    _ => raw.upper_bound as f64,
                };
                IsolationWindow {
                    target,
                    lower,
                    upper,
                }
            }
            None => {
                if !warned_missing_precursor {
                    log::warn!(
                        "MS2 spectrum at scan_index {} has no precursor info; skipping",
                        idx
                    );
                    warned_missing_precursor = true;
                }
                return None;
            }
        };

        let retention_time = spectrum.start_time();
        let peaks = extract_peaks(&mut spectrum);

        if ms2_count == 0 {
            log::info!(
                "First MS2 scan: scan_index={} rt={:.2} min isolation={:.4}({:.4}-{:.4}) peaks={}",
                idx,
                retention_time,
                iw.target,
                iw.lower,
                iw.upper,
                peaks.len()
            );
        }
        ms2_count += 1;

        if peaks.is_empty() {
            return None;
        }
        Some(Spectrum {
            scan_index: idx,
            retention_time,
            peaks,
            ms_level: 2,
            isolation_window: Some(iw),
            faims_cv: None,
        })
    })
}

/// Read the PSI-MS FAIMS compensation-voltage parameter.
///
/// Converters emit MS:1001581 either directly on the spectrum (including the
/// ThermoRawFileParser mzML used by the benchmark) or on the first acquisition
/// scan event. Accept both encodings. Treat malformed and non-finite values as
/// absent so they cannot become unstable detector-channel keys.
fn extract_faims_cv(spectrum: &MultiLayerSpectrum) -> Option<f32> {
    spectrum
        .get_param_by_accession("MS:1001581")
        .or_else(|| {
            spectrum
                .acquisition()
                .first_scan()
                .and_then(|scan| scan.get_param_by_accession("MS:1001581"))
        })
        .and_then(|param| param.value().to_f32().ok())
        .filter(|cv| cv.is_finite())
        .map(|cv| if cv == 0.0 { 0.0 } else { cv })
}

/// Extract peaks from a spectrum, centroiding on the fly if profile-mode.
fn extract_peaks(spectrum: &mut MultiLayerSpectrum) -> Vec<Peak> {
    if matches!(spectrum.signal_continuity(), SignalContinuity::Profile) {
        // Pick peaks using mzdata's built-in quadratic peak fitter (SNR >= 1.0).
        if let Err(e) = spectrum.pick_peaks(1.0) {
            log::warn!("Peak picking failed: {e}");
            return Vec::new();
        }
        if let Some(peaks) = spectrum.peaks.as_ref() {
            let mut result: Vec<Peak> = peaks
                .iter()
                .filter(|p| p.intensity() > 0.0)
                .map(|p| Peak {
                    mz: p.mz() as f32,
                    intensity: p.intensity(),
                    ion_mobility: 0.0,
                })
                .collect();
            result.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));
            return result;
        }
        return Vec::new();
    }

    // Already centroided — read directly from the raw arrays.
    if let Some(arrays) = spectrum.raw_arrays() {
        if let (Ok(mzs), Ok(intensities)) = (arrays.mzs(), arrays.intensities()) {
            if !mzs.is_empty() {
                let mut peaks: Vec<Peak> = mzs
                    .iter()
                    .zip(intensities.iter())
                    .filter(|(_, &i)| i > 0.0)
                    .map(|(&mz, &intensity)| Peak {
                        mz: mz as f32,
                        intensity,
                        ion_mobility: 0.0,
                    })
                    .collect();
                peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));
                return peaks;
            }
        }
    }

    Vec::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use mzdata::params::{ControlledVocabulary, Param};
    use mzdata::spectrum::ScanEvent;

    fn faims_param(value: f64) -> Param {
        Param::builder()
            .name("FAIMS compensation voltage")
            .value(value)
            .controlled_vocabulary(ControlledVocabulary::MS)
            .accession(1001581)
            .build()
    }

    #[test]
    fn extracts_spectrum_level_faims_compensation_voltage() {
        let mut spectrum = MultiLayerSpectrum::default();
        spectrum.add_param(faims_param(-50.0));

        assert_eq!(extract_faims_cv(&spectrum), Some(-50.0));
    }

    #[test]
    fn extracts_scan_level_faims_compensation_voltage() {
        let mut spectrum = MultiLayerSpectrum::default();
        spectrum.description.acquisition.scans.push(ScanEvent {
            params: Some(Box::new(vec![faims_param(-65.0)])),
            ..ScanEvent::default()
        });

        assert_eq!(extract_faims_cv(&spectrum), Some(-65.0));
    }

    #[test]
    fn missing_or_non_finite_faims_voltage_is_absent() {
        let spectrum = MultiLayerSpectrum::default();
        assert_eq!(extract_faims_cv(&spectrum), None);

        let mut spectrum = MultiLayerSpectrum::default();
        spectrum.description.acquisition.scans.push(ScanEvent {
            params: Some(Box::new(vec![faims_param(f64::NAN)])),
            ..ScanEvent::default()
        });
        assert_eq!(extract_faims_cv(&spectrum), None);
    }
}
