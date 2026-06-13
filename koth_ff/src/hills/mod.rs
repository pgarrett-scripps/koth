pub mod active;
pub mod calibration;
pub mod detector;
pub mod filter;
pub mod noise;
pub mod smooth;
pub mod split;

use std::collections::HashMap;

use crate::config::{FileConfig, HillsConfig};
use crate::models::{Hill, IsolationWindow, Spectrum};

use detector::HillDetector;

/// Assign sequential `hill_id`s (0..hills.len()) in the order the slice is
/// currently in. Called once after detection (and splitting) completes so the
/// IDs are stable for the lifetime of the run — features then reference these
/// IDs via the `hill_ids` column.
pub fn assign_hill_ids(hills: &mut [Hill]) {
    for (i, h) in hills.iter_mut().enumerate() {
        h.hill_id = i as u64;
    }
}

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
            config.min_prominence,
        );
        log::info!("{} hills after co-elution splitting", hills.len());
    }

    if config.filter_large_baseline_hills {
        hills = filter::filter_large_baseline_hills(
            hills,
            config.large_hill_min_scans,
            config.large_hill_peak_factor,
        );
    }

    assign_hill_ids(&mut hills);
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

/// Detect chromatographic hills from MS2 spectra, grouping spectra by their
/// precursor isolation window. Each isolation window gets its own detector
/// (with its own per-channel scan indexing), so a hill never spans two windows.
///
/// Spectra without an isolation window or with `ms_level != 2` are ignored.
pub fn detect_ms2_hills_from_iter<I>(
    spectra: I,
    config: &HillsConfig,
    file: &FileConfig,
) -> Vec<Hill>
where
    I: Iterator<Item = Spectrum>,
{
    let mut detectors: HashMap<(i64, i64, i64), (IsolationWindow, HillDetector, usize)> =
        HashMap::new();
    let mut total = 0usize;
    let mut total_peaks = 0usize;

    for mut spectrum in spectra {
        if spectrum.ms_level != 2 {
            continue;
        }
        let Some(iw) = spectrum.isolation_window else {
            continue;
        };
        if let Some(sigma) = file.noise_filter_sigma {
            noise::filter_spectrum(&mut spectrum, sigma);
        }
        total_peaks += spectrum.peaks.len();
        total += 1;

        let key = iw.key();
        let entry = detectors.entry(key).or_insert_with(|| {
            (
                iw,
                HillDetector::new(config, file).with_isolation_window(iw),
                0,
            )
        });
        // Re-index scan_index per isolation window so consecutive cycles
        // are consecutive scans for the gap-tracking logic.
        let per_channel_idx = entry.2;
        entry.2 += 1;
        let local = Spectrum {
            scan_index: per_channel_idx,
            ..spectrum
        };
        entry.1.process_scan(&local);
    }

    let n_channels = detectors.len();
    log::info!(
        "MS2 hill detection: {} spectra across {} isolation windows ({:.1}M peaks)",
        total,
        n_channels,
        total_peaks as f64 / 1_000_000.0
    );

    let mut hills = Vec::new();
    for (_key, (_iw, det, _n)) in detectors {
        let mut window_hills = det.finish();
        if config.split_hills {
            window_hills = split::split_coeluting(
                window_hills,
                config.min_peak_distance,
                config.min_peak_height,
                config.min_scans,
                config.min_prominence,
            );
        }
        if config.filter_large_baseline_hills {
            window_hills = filter::filter_large_baseline_hills(
                window_hills,
                config.large_hill_min_scans,
                config.large_hill_peak_factor,
            );
        }
        hills.extend(window_hills);
    }

    assign_hill_ids(&mut hills);
    log::info!(
        "MS2 hill detection complete: {} hills across {} isolation windows",
        hills.len(),
        n_channels
    );
    hills
}
