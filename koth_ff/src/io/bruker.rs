/// Bruker .d folder reader using timsrust.
///
/// Per-frame pipeline:
///   1. `dnoise::FlatFrame::from_frame` — flatten the timsrust frame into
///      integer `(scan_idx, tof_idx, intensity)` arrays.
///   2. `dnoise::filter_iterated` — vertical-IM feature filter (Stage 1).
///   3. `dnoise::watershed_centroid` — watershed centroider (Stage 3),
///      emitting integer `(scan, tof, intensity)` triples.
///   4. Convert each centroid to `(mz, 1/K0, intensity)` via the timsrust
///      `Tof2MzConverter` / `Scan2ImConverter`. One conversion call per
///      centroid (vs. per raw member in the previous local implementation).
///   5. Optional `bruker_noise_sigma` MAD filter on the emitted Spectrum.
///
/// All denoising / centroiding logic now lives in the `dnoise` crate
/// — this file is thin I/O glue around it.
#[cfg(feature = "tdf")]
pub mod inner {
    use std::cmp::Ordering;
    use std::path::Path;
    use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

    use rayon::prelude::*;
    use timsrust::converters::ConvertableDomain;
    use timsrust::readers::FrameReader;

    use dnoise::{
        filter_iterated, watershed_centroid, FilterParams, FlatFrame, WatershedParams,
    };

    use crate::config::FileConfig;
    use crate::error::KothError;
    use crate::hills::noise;
    use crate::models::{Peak, Spectrum};

    const MAX_PEAKS: usize = 10_000;

    pub fn read_bruker(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
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

        let mut spectra: Vec<Spectrum> = frame_reader
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

        if spectra.is_empty() {
            return Err(KothError::NoSpectra);
        }

        spectra.sort_by(|a, b| a.retention_time.partial_cmp(&b.retention_time).unwrap());
        for (i, s) in spectra.iter_mut().enumerate() {
            s.scan_index = i;
        }

        log::info!("Read {} MS1 spectra from Bruker .d", spectra.len());
        Ok(spectra)
    }
}

#[cfg(feature = "tdf")]
pub use inner::read_bruker;
