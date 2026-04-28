use std::path::Path;
use std::{fs, io};

use flate2::read::GzDecoder;
use mzdata::prelude::*;
use mzdata::spectrum::MultiLayerSpectrum;
use mzdata::MZReader;

use crate::error::KothError;
use crate::models::{Peak, Spectrum};

/// Read an mzML file (plain or gzip-compressed) and return all MS1 spectra
/// sorted by retention time.
pub fn read_mzml(path: &Path) -> Result<Vec<Spectrum>, KothError> {
    let mut spectra: Vec<Spectrum> = collect_ms1(open_reader(path)?);

    if spectra.is_empty() {
        return Err(KothError::NoSpectra);
    }

    spectra.sort_by(|a, b| a.retention_time.partial_cmp(&b.retention_time).unwrap());
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
pub fn stream_mzml(path: &Path) -> Result<Box<dyn Iterator<Item = Spectrum>>, KothError> {
    crate::mem::log_mem("before open_reader (stream_mzml)");
    Ok(Box::new(ms1_stream(open_reader(path)?)))
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Boxed iterator over MultiLayerSpectrum that works for both compressed and
/// uncompressed readers.
type BoxedRawIter = Box<dyn Iterator<Item = MultiLayerSpectrum>>;

fn open_reader(path: &Path) -> Result<BoxedRawIter, KothError> {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .to_lowercase();

    if name.ends_with(".gz") {
        // Decompress the full file into memory so mzdata gets a seekable Cursor.
        // open_gzipped_read uses a limited PreBufferedStream that corrupts large
        // base64 arrays (e.g. zlib-compressed Orbitrap data) during streaming.
        let file = fs::File::open(path).map_err(|e| KothError::MzmlError(e.to_string()))?;
        let mut decoder = GzDecoder::new(io::BufReader::new(file));
        let mut buf = Vec::new();
        io::Read::read_to_end(&mut decoder, &mut buf)
            .map_err(|e| KothError::MzmlError(format!("gzip decompress: {e}")))?;
        let cursor = io::Cursor::new(buf);
        let reader = MZReader::open_read_seek(cursor)
            .map_err(|e| KothError::MzmlError(e.to_string()))?;
        Ok(Box::new(reader))
    } else {
        let reader = MZReader::open_path(path)
            .map_err(|e| KothError::MzmlError(e.to_string()))?;
        Ok(Box::new(reader))
    }
}

fn collect_ms1(reader: BoxedRawIter) -> Vec<Spectrum> {
    let mut out = Vec::new();
    for (scan_index, spectrum) in reader.enumerate() {
        if spectrum.ms_level() != 1 {
            continue;
        }
        let peaks = extract_peaks(&spectrum);
        if peaks.is_empty() {
            continue;
        }
        out.push(Spectrum {
            scan_index,
            retention_time: spectrum.start_time(),
            peaks,
        });
    }
    out
}

fn ms1_stream(reader: BoxedRawIter) -> impl Iterator<Item = Spectrum> {
    let mut scan_index = 0usize;
    let mut ms1_count = 0usize;

    reader.filter_map(move |spectrum| {
        let idx = scan_index;
        scan_index += 1;

        if spectrum.ms_level() != 1 {
            return None;
        }

        if ms1_count == 0 {
            let continuity = spectrum.signal_continuity();
            log::info!("First MS1 scan signal continuity: {:?}", continuity);
            if format!("{:?}", continuity).to_lowercase().contains("profile") {
                log::error!(
                    "PROFILE-MODE DATA DETECTED. This algorithm requires centroided spectra. \
                     Convert with: msconvert input.raw --filter \"peakPicking true 1-\" --mzML"
                );
            }
        }

        let retention_time = spectrum.start_time();
        let peaks = extract_peaks(&spectrum);

        if ms1_count == 0 {
            log::info!(
                "First MS1 scan: scan_index={} rt={:.2} min  peaks={}",
                idx, retention_time, peaks.len()
            );
            if peaks.len() > 50_000 {
                log::error!(
                    "First MS1 scan has {} peaks — this is characteristic of profile-mode data \
                     and will cause a memory explosion. Centroid first.",
                    peaks.len()
                );
            }
        }

        ms1_count += 1;
        if peaks.is_empty() {
            return None;
        }
        Some(Spectrum { scan_index: idx, retention_time, peaks })
    })
}

fn extract_peaks(spectrum: &impl SpectrumLike) -> Vec<Peak> {
    if let Some(arrays) = spectrum.raw_arrays() {
        if let (Ok(mzs), Ok(intensities)) = (arrays.mzs(), arrays.intensities()) {
            if !mzs.is_empty() {
                let mut peaks: Vec<Peak> = mzs
                    .iter()
                    .zip(intensities.iter())
                    .filter(|(_, &i)| i > 0.0)
                    .map(|(&mz, &intensity)| Peak {
                        mz: mz as f32,
                        intensity: intensity as f32,
                        ion_mobility: 0.0,
                    })
                    .collect();
                peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
                return peaks;
            }
        }
    }

    Vec::new()
}
