use std::path::Path;
use std::sync::{Arc, Mutex};
use std::{fs, io};

use flate2::bufread::MultiGzDecoder;
use mzdata::io::mzml::MzMLReaderType;
use mzdata::prelude::*;
use mzdata::spectrum::{MultiLayerSpectrum, SignalContinuity};
use mzdata::MZReader;

use crate::error::KothError;
use crate::models::{IsolationWindow, Peak, Spectrum};

/// Read an mzML file (plain or gzip-compressed) and return all MS1 spectra
/// sorted by retention time.
pub fn read_mzml(path: &Path) -> Result<Vec<Spectrum>, KothError> {
    let (reader, status) = open_reader(path)?;
    let mut spectra: Vec<Spectrum> = collect_ms1(reader);
    status.check()?;

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
///
/// A gzip error part-way through cannot be returned from this iterator; it is
/// logged and the stream ends early. [`crate::io::stream_spectra`] returns it
/// as an `Err` item instead.
pub fn stream_mzml(path: &Path) -> Result<Box<dyn Iterator<Item = Spectrum> + Send>, KothError> {
    crate::mem::log_mem("before open_reader (stream_mzml)");
    let (reader, status) = open_reader(path)?;
    Ok(Box::new(status.log_at_end(ms1_stream(reader))))
}

/// [`stream_mzml`] with a read error as a final `Err` item.
pub(crate) fn stream_mzml_results(path: &Path) -> Result<super::SpectrumStream, KothError> {
    crate::mem::log_mem("before open_reader (stream_mzml)");
    let (reader, status) = open_reader(path)?;
    let tail = std::iter::once_with(move || status.check()).filter_map(Result::err);
    Ok(Box::new(ms1_stream(reader).map(Ok).chain(tail.map(Err))))
}

/// Stream MS2 spectra (with isolation window metadata) from an mzML file.
///
/// Only spectra with `ms_level == 2` and a parsable precursor isolation window
/// are yielded. The `scan_index` on each yielded spectrum is its absolute
/// position in the file, so callers can later re-index per-isolation-window.
/// A gzip error part-way through is logged and ends the stream early, as for
/// [`stream_mzml`].
pub fn stream_mzml_ms2(
    path: &Path,
) -> Result<Box<dyn Iterator<Item = Spectrum> + Send>, KothError> {
    let (iter, status) = stream_mzml_ms2_checked(path)?;
    Ok(Box::new(status.log_at_end(iter)))
}

/// [`stream_mzml_ms2`] plus the [`ReadStatus`] to check once it is drained.
pub(crate) fn stream_mzml_ms2_checked(
    path: &Path,
) -> Result<(Box<dyn Iterator<Item = Spectrum> + Send>, ReadStatus), KothError> {
    crate::mem::log_mem("before open_reader (stream_mzml_ms2)");
    let (reader, status) = open_reader(path)?;
    Ok((Box::new(ms2_stream(reader)), status))
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Boxed iterator over MultiLayerSpectrum that works for both compressed and
/// uncompressed readers.
type BoxedRawIter = Box<dyn Iterator<Item = MultiLayerSpectrum> + Send>;

/// The first I/O error the gzip stream hit, if any.
///
/// mzdata's streaming reader ends iteration on an I/O error without reporting
/// it, so a truncated or corrupt `.mzML.gz` would otherwise read as a shorter
/// run. The decoder is wrapped so the error is kept here, and every consumer
/// checks it once the spectra are drained. Always empty for plain mzML.
#[derive(Clone, Default)]
pub(crate) struct ReadStatus {
    path: std::path::PathBuf,
    error: Arc<Mutex<Option<String>>>,
}

impl ReadStatus {
    pub(crate) fn check(&self) -> Result<(), KothError> {
        match self.error.lock().ok().and_then(|e| e.clone()) {
            Some(e) => Err(KothError::MzmlError(format!(
                "gzip decompress of '{}' failed part-way through: {e}",
                self.path.display()
            ))),
            None => Ok(()),
        }
    }

    fn log_at_end<I: Iterator<Item = Spectrum> + Send>(
        self,
        iter: I,
    ) -> impl Iterator<Item = Spectrum> + Send {
        iter.chain(std::iter::from_fn(move || {
            if let Err(e) = self.check() {
                log::error!("{e}");
            }
            None
        }))
    }
}

/// A reader that records its first error in a [`ReadStatus`].
struct TrackedRead<R> {
    inner: R,
    error: Arc<Mutex<Option<String>>>,
}

impl<R: io::Read> io::Read for TrackedRead<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf).inspect_err(|e| {
            if let Ok(mut slot) = self.error.lock() {
                slot.get_or_insert_with(|| e.to_string());
            }
        })
    }
}

fn open_reader(path: &Path) -> Result<(BoxedRawIter, ReadStatus), KothError> {
    let status = ReadStatus {
        path: path.to_path_buf(),
        ..ReadStatus::default()
    };
    // Route the compressed-vs-plain decision through the single canonical format
    // predicate. `open_reader` is only reached for paths already classified as
    // mzML, so anything not `MzmlGz` here is plain (uncompressed) mzML.
    if super::detect_format(path) == super::InputFormat::MzmlGz {
        // Stream the decompressed bytes straight into the parser: the spectra
        // are only ever iterated forward, and `MzMLReaderType`'s `Iterator`
        // needs `io::Read` alone. Every caller (MS1 read, MS1 stream, the DIA
        // MS2 pass) opens the path afresh, so nothing seeks or reopens this
        // handle. This replaces a `read_to_end` of the whole decompressed file
        // into a `Cursor`, which held the whole decompressed file in memory
        // (2.5 GB for one PXD066701 Astral DIA run).
        //
        // Do not route this through mzdata's `open_gzipped_read`: its limited
        // `PreBufferedStream` corrupts large base64 arrays (e.g. zlib-compressed
        // Orbitrap data). `MultiGzDecoder` is what mzdata's own conversion
        // pipeline wraps (via `RestartableGzDecoder`), and it also reads
        // multi-member (bgzip-style) files to the end.
        let file = fs::File::open(path).map_err(|e| KothError::MzmlError(e.to_string()))?;
        let decoder = TrackedRead {
            inner: MultiGzDecoder::new(io::BufReader::new(file)),
            error: Arc::clone(&status.error),
        };
        let reader: MzMLReaderType<_> = MzMLReaderType::new(decoder);
        Ok((Box::new(reader), status))
    } else {
        let reader = MZReader::open_path(path).map_err(|e| KothError::MzmlError(e.to_string()))?;
        Ok((Box::new(reader), status))
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

    /// A `.mzML.gz` cut off mid-stream must fail, not read as a shorter run:
    /// mzdata's streaming reader ends iteration on an I/O error silently.
    #[test]
    fn truncated_gzip_is_an_error() {
        use flate2::{write::GzEncoder, Compression};
        use std::io::Write;

        let body = "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<mzML>\n".to_string()
            + &"<!-- padding so the member spans several deflate blocks -->\n".repeat(20_000);
        let mut enc = GzEncoder::new(Vec::new(), Compression::fast());
        enc.write_all(body.as_bytes()).unwrap();
        let gz = enc.finish().unwrap();

        let dir = std::env::temp_dir().join(format!("koth_mzml_gz_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("truncated.mzML.gz");
        fs::write(&path, &gz[..gz.len() / 2]).unwrap();

        let err = read_mzml(&path).expect_err("truncated gzip must not read cleanly");
        let _ = fs::remove_dir_all(&dir);
        assert!(
            matches!(&err, KothError::MzmlError(m) if m.contains("gzip")),
            "unexpected error: {err}"
        );
    }
}
