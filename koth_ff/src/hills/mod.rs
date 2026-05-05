pub mod active;
pub mod detector;
pub mod noise;
pub mod split;

use crate::config::{FileConfig, HillsConfig};
use crate::models::{Hill, Spectrum};

use detector::HillDetector;

/// Detect hills from an iterator of spectra.
///
/// Spectra are processed one at a time and dropped immediately —
/// the full dataset is never held in memory simultaneously.
/// Assumes spectra arrive in ascending retention-time order.
pub fn detect_hills_from_iter<I>(spectra: I, config: &HillsConfig, file: &FileConfig) -> Vec<Hill>
where
    I: Iterator<Item = Spectrum>,
{
    let mut det = HillDetector::new(config, file);
    let mut count = 0usize;
    let mut total_peaks = 0usize;
    for mut spectrum in spectra {
        if let Some(sigma) = file.noise_filter_sigma {
            noise::filter_spectrum(&mut spectrum, sigma);
        }
        total_peaks += spectrum.peaks.len();
        det.process_scan(&spectrum);
        count += 1;
    }
    if count > 0 {
        let avg_peaks = total_peaks / count;
        log::info!(
            "Hill detection: processed {} spectra, {:.1}M total peaks, avg {}/scan",
            count,
            total_peaks as f64 / 1_000_000.0,
            avg_peaks
        );
        if avg_peaks > 50_000 {
            log::warn!(
                "Average {} peaks/scan is very high — input may be profile-mode mzML. \
                 Consider centroiding with msconvert --filter 'peakPicking true 1-' first.",
                avg_peaks
            );
        }
    }

    let mut hills = det.finish();
    log::info!("{} hills before splitting", hills.len());

    if config.split_hills {
        hills = split::split_coeluting(
            hills,
            config.min_peak_distance,
            config.min_peak_height,
            config.min_scans,
            config.min_valley_ratio,
        );
        log::info!("{} hills after co-elution splitting", hills.len());
    }

    log::info!("Hill detection complete: {} hills", hills.len());
    hills
}

/// Detect chromatographic hills from a pre-loaded slice of spectra.
/// Prefer `detect_hills_from_iter` for large files to avoid holding
/// all spectra in memory.
pub fn detect_hills(spectra: &[Spectrum], config: &HillsConfig, file: &FileConfig) -> Vec<Hill> {
    log::info!("Starting hill detection on {} spectra", spectra.len());
    detect_hills_from_iter(spectra.iter().cloned(), config, file)
}
