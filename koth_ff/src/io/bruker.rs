/// Bruker .d folder reader. Two per-frame paths, selected by `bruker_streaming`:
///
/// * **local** (default, `read_bruker_local`): timsrust frame iteration, then
///   dnoise's `filter_iterated` (vertical-IM) + `watershed_centroid`, converting
///   each centroid to `(m/z, 1/K0, intensity)` via the timsrust converters. This
///   is the paper's validated pipeline (no halo).
/// * **streaming** (opt-in, `read_bruker_streaming`): dnoise's `RunContext` runs
///   the vertical filter + horizontal halo + watershed in a single in-process
///   pass over the raw frames — the exact stage code the standalone `dnoise` tool
///   runs — with no denoised `.d` written to disk. Calibration comes from the
///   context. This is the future default once re-validated on the Bruker cohort.
///
/// Both paths then optionally apply the `bruker_noise_sigma` MAD filter, sort by
/// retention time, and reassign sequential scan indices (`finalize`). All
/// denoising / centroiding logic lives in the `dnoise` crate.
#[cfg(feature = "tdf")]
pub mod inner {
    use std::cmp::Ordering;
    use std::path::Path;
    use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

    use rayon::prelude::*;
    use timsrust::converters::ConvertableDomain;
    use timsrust::readers::FrameReader;

    use dnoise::{
        filter_iterated, watershed_centroid, FilterParams, FlatFrame, HaloParams, RunContext,
        Stages, WatershedParams,
    };

    use crate::config::FileConfig;
    use crate::error::KothError;
    use crate::hills::noise;
    use crate::models::{Peak, Spectrum};

    const MAX_PEAKS: usize = 10_000;

    pub fn read_bruker(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
        if file.bruker_streaming {
            read_bruker_streaming(path, file)
        } else {
            read_bruker_local(path, file)
        }
    }

    /// Sort spectra by retention time and reassign sequential scan indices. Shared
    /// by both reader paths so the downstream hill detector sees a canonical order.
    fn finalize(mut spectra: Vec<Spectrum>) -> Result<Vec<Spectrum>, KothError> {
        if spectra.is_empty() {
            return Err(KothError::NoSpectra);
        }
        spectra.sort_by(|a, b| {
            a.retention_time
                .partial_cmp(&b.retention_time)
                .unwrap_or(Ordering::Equal)
        });
        for (i, s) in spectra.iter_mut().enumerate() {
            s.scan_index = i;
        }
        Ok(spectra)
    }

    /// In-process streaming path (opt-in via `bruker_streaming`): drive dnoise's
    /// [`RunContext`], which runs the vertical-IM filter + horizontal halo +
    /// watershed in a single pass over each raw MS1 frame — the exact stage code
    /// the standalone `dnoise` tool runs, with no denoised `.d` written to disk.
    /// Calibration comes straight from the context, so no timsrust converters here.
    fn read_bruker_streaming(
        path: &Path,
        file: &FileConfig,
    ) -> Result<Vec<Spectrum>, KothError> {
        let filter_params = FilterParams {
            mz_half_width: file.bruker_filter_mz_half_width,
            min_feature_length: file.bruker_filter_min_feature_length,
            max_internal_gap: file.bruker_filter_max_internal_gap,
            min_window_intensity: file.bruker_filter_min_window_intensity,
            min_feature_intensity: file.bruker_filter_min_feature_intensity,
            num_iterations: file.bruker_filter_num_iterations,
        };
        let watershed_params = WatershedParams {
            box_scan: file.bruker_watershed_box_scan,
            box_mz_idx: file.bruker_watershed_box_mz_idx,
            min_seed_intensity: file.bruker_watershed_min_seed_intensity,
            min_centroid_total: file.bruker_watershed_min_centroid_total,
            max_tof_offset: file.bruker_watershed_max_tof_offset,
        };
        let halo_params = HaloParams {
            peak_fraction: file.bruker_halo_peak_fraction,
            mz_idx_half_width: file.bruker_halo_mz_idx_half_width,
            scan_half_width: file.bruker_halo_scan_half_width,
        };
        // Vertical filter + optional halo + watershed, all in one in-process pass.
        let stages = Stages {
            halo: file.bruker_halo.then_some(&halo_params),
            watershed: Some(&watershed_params),
            ..Stages::default()
        };

        let ctx = RunContext::open(path, &filter_params, &stages)
            .map_err(|e| KothError::TdfError(e.to_string()))?;
        let cal = ctx.calibration();
        let noise_sigma = file.bruker_noise_sigma;

        let spectra: Vec<Spectrum> = (0..ctx.len())
            .into_par_iter()
            .filter(|&i| ctx.is_ms1(i))
            .map(|i| -> Result<Spectrum, KothError> {
                let decoded = ctx
                    .process(i)
                    .map_err(|e| KothError::TdfError(e.to_string()))?;
                let mut peaks: Vec<Peak> = decoded
                    .survivors
                    .iter()
                    .map(|&(scan, tof, intensity)| Peak {
                        mz: cal.tof_to_mz(tof) as f32,
                        intensity: intensity as f32,
                        ion_mobility: cal.scan_to_im(scan) as f32,
                    })
                    .collect();
                peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));

                let mut spectrum = Spectrum {
                    scan_index: decoded.frame_id,
                    retention_time: decoded.rt_seconds / 60.0,
                    peaks,
                    ms_level: 1,
                    isolation_window: None,
                };
                if let Some(sigma) = noise_sigma {
                    noise::filter_spectrum(&mut spectrum, sigma);
                }
                Ok(spectrum)
            })
            .collect::<Result<Vec<_>, _>>()?;

        log::info!(
            "Read {} MS1 spectra from Bruker .d (dnoise streaming)",
            spectra.len()
        );
        finalize(spectra)
    }

    /// Historical local path (default): timsrust frame iteration + dnoise's
    /// vertical filter and watershed (no halo), matching the paper's validated
    /// pipeline. Kept as the default until streaming is re-validated on the
    /// Bruker cohort.
    fn read_bruker_local(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
        let path_str = path
            .to_str()
            .ok_or_else(|| KothError::UnsupportedFormat("non-UTF8 path".into()))?;

        let frame_reader =
            FrameReader::new(path_str).map_err(|e| KothError::TdfError(e.to_string()))?;
        let metadata = timsrust::readers::MetadataReader::new(path_str)
            .map_err(|e| KothError::TdfError(e.to_string()))?;

        let mz_converter = metadata.mz_converter;
        let ims_converter = metadata.im_converter;

        let filter_params = FilterParams {
            mz_half_width: file.bruker_filter_mz_half_width,
            min_feature_length: file.bruker_filter_min_feature_length,
            max_internal_gap: file.bruker_filter_max_internal_gap,
            min_window_intensity: file.bruker_filter_min_window_intensity,
            min_feature_intensity: file.bruker_filter_min_feature_intensity,
            num_iterations: file.bruker_filter_num_iterations,
        };
        let watershed_params = WatershedParams {
            box_scan: file.bruker_watershed_box_scan,
            box_mz_idx: file.bruker_watershed_box_mz_idx,
            min_seed_intensity: file.bruker_watershed_min_seed_intensity,
            min_centroid_total: file.bruker_watershed_min_centroid_total,
            max_tof_offset: file.bruker_watershed_max_tof_offset,
        };
        let noise_sigma = file.bruker_noise_sigma;

        let processed = AtomicUsize::new(0);
        const PROGRESS_INTERVAL: usize = 100;

        let spectra: Vec<Spectrum> = frame_reader
            .parallel_filter(|f| f.ms_level == timsrust::MSLevel::MS1)
            .filter_map(|frame_result| match frame_result {
                Ok(frame) => {
                    let flat = FlatFrame::from_frame(&frame);
                    let n_raw = flat.len();
                    let keep = filter_iterated(&flat, &filter_params);
                    let survivors = flat.survivors(&keep);
                    let n_filtered = survivors.len();

                    let centroids = watershed_centroid(
                        &survivors,
                        &watershed_params,
                        MAX_PEAKS,
                    );
                    let n_centroids = centroids.len();

                    let mut out_peaks: Vec<Peak> = centroids
                        .into_iter()
                        .map(|(scan, tof, intensity)| Peak {
                            mz: mz_converter.convert(tof as f64) as f32,
                            intensity: intensity as f32,
                            ion_mobility: ims_converter.convert(scan as f64) as f32,
                        })
                        .collect();
                    out_peaks.sort_by(|a, b| {
                        a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal)
                    });

                    let mut spectrum = Spectrum {
                        scan_index: frame.index,
                        retention_time: frame.rt_in_seconds as f64 / 60.0,
                        peaks: out_peaks,
                        ms_level: 1,
                        isolation_window: None,
                    };

                    if let Some(sigma) = noise_sigma {
                        noise::filter_spectrum(&mut spectrum, sigma);
                    }
                    let n_after_noise = spectrum.peaks.len();

                    let n = processed.fetch_add(1, AtomicOrdering::Relaxed) + 1;
                    if n % PROGRESS_INTERVAL == 0 {
                        if noise_sigma.is_some() {
                            log::info!(
                                "Bruker frame {}: {} raw -> {} filtered -> {} centroids -> {} post-MAD",
                                n,
                                n_raw,
                                n_filtered,
                                n_centroids,
                                n_after_noise,
                            );
                        } else {
                            log::info!(
                                "Bruker frame {}: {} raw -> {} filtered -> {} centroids",
                                n,
                                n_raw,
                                n_filtered,
                                n_centroids,
                            );
                        }
                    }

                    Some(spectrum)
                }
                Err(e) => {
                    log::error!("Error parsing Bruker frame: {:?}", e);
                    None
                }
            })
            .collect();

        log::info!("Read {} MS1 spectra from Bruker .d", spectra.len());
        finalize(spectra)
    }
}

#[cfg(feature = "tdf")]
pub use inner::read_bruker;
