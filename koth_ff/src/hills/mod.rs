pub mod active;
pub mod detector;
pub mod filter;
pub mod noise;
pub mod smooth;
pub mod split;
pub mod tic_norm;

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
    // TIC anomaly normalisation needs the full set of MS1 spectra in advance
    // (centred rolling median), so it forces collection into memory. If the
    // option is off, we keep the streaming path unchanged.
    if config.tic_norm_window > 0 {
        let mut buf: Vec<Spectrum> = spectra.collect();
        tic_norm::normalize_tic_in_place(
            &mut buf,
            config.tic_norm_window,
            config.tic_norm_min_scale,
            config.tic_norm_max_scale,
            tic_norm::RefMode::parse(&config.tic_norm_mode),
        );
        return detect_hills_from_iter_inner(buf.into_iter(), config, file);
    }
    detect_hills_from_iter_inner(spectra, config, file)
}

fn detect_hills_from_iter_inner<I>(spectra: I, config: &HillsConfig, file: &FileConfig) -> Vec<Hill>
where
    I: Iterator<Item = Spectrum>,
{
    let mut det = HillDetector::new(config, file);
    let mut count = 0usize;
    let mut total_peaks = 0usize;
    // Split the streaming loop into read (decompress + XML parse, inside the
    // iterator's next()) vs detect (process_scan) to locate the hotspot.
    let mut read_ns = std::time::Duration::ZERO;
    let mut detect_ns = std::time::Duration::ZERO;
    let mut iter = spectra;
    loop {
        let t0 = std::time::Instant::now();
        let next = iter.next();
        read_ns += t0.elapsed();
        let Some(mut spectrum) = next else { break };
        let t1 = std::time::Instant::now();
        if let Some(sigma) = file.noise_filter_sigma {
            noise::filter_spectrum(&mut spectrum, sigma);
        }
        total_peaks += spectrum.peaks.len();
        det.process_scan(&spectrum);
        detect_ns += t1.elapsed();
        count += 1;
    }
    if let Some(avg_peaks) = total_peaks.checked_div(count) {
        log::info!(
            "Hill detection: processed {} spectra, {:.1}M total peaks, avg {}/scan \
             [read(decompress+parse) {:.2?}, detect(process_scan) {:.2?}]",
            count,
            total_peaks as f64 / 1_000_000.0,
            avg_peaks,
            read_ns,
            detect_ns
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
        hills = split::split_coeluting(hills, config);
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
    // `detect_hills_from_iter` performs TIC normalization itself when
    // `tic_norm_window > 0`; delegate directly so the spectra are normalized
    // exactly once (normalizing here too would double-scale every intensity).
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
    // Sort by isolation-window key so hill collection order (and the IDs
    // assigned by position below) is deterministic across runs — HashMap
    // iteration order is randomized per process.
    let mut detectors: Vec<_> = detectors.into_iter().collect();
    detectors.sort_by_key(|(key, _)| *key);
    for (_key, (_iw, det, _n)) in detectors {
        let mut window_hills = det.finish();
        if config.split_hills {
            window_hills = split::split_coeluting(window_hills, config);
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
