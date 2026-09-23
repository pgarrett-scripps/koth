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
///   context. The flag selects preprocessing; both modes stream into detection.
///
/// Both paths then optionally apply the `bruker_noise_sigma` MAD filter, emit in
/// retention-time order, and assign sequential scan indices. All
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
        filter_iterated, watershed_centroid, FilterParams, FlatFrame, HaloParams, Ms1PolygonParams,
        RunContext, Stages, WatershedParams,
    };

    use crate::config::FileConfig;
    use crate::error::KothError;
    use crate::hills::noise;
    use crate::io::tims_calibration::{self, MobilityScale, ScanToMobility};
    use crate::models::{IsolationWindow, Peak, Spectrum};

    const MAX_PEAKS: usize = 10_000;

    /// The run's acquisition calibration when `bruker_mobility_scale` is
    /// `calibrated` (the default); `None` keeps the legacy linear converter of
    /// the reader in use. A run whose calibration cannot be read is an error,
    /// never a silent fallback to the linear scale.
    fn calibrated_mobility(
        path: &Path,
        file: &FileConfig,
    ) -> Result<Option<ScanToMobility>, KothError> {
        match file.bruker_mobility_scale {
            MobilityScale::Linear => Ok(None),
            MobilityScale::Calibrated => tims_calibration::load(path, MobilityScale::Calibrated)
                .map(Some)
                .map_err(KothError::TdfError),
        }
    }

    /// Collect the streaming reader for callers that explicitly need all spectra.
    pub fn read_bruker(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
        stream_bruker(path, file).collect()
    }

    /// Stream MS1 spectra in retention-time order with bounded parallel decoding.
    /// Only frame metadata is indexed up front. Both preprocessing modes use
    /// batches of at most 32 frames and a queue of at most 32 spectra.
    /// Opening/decoding failures are yielded as errors; dropping the iterator
    /// stops the producer and joins its thread.
    pub fn stream_bruker(
        path: &Path,
        file: &FileConfig,
    ) -> super::super::prefetch::Prefetch<Result<Spectrum, KothError>> {
        let path = path.to_owned();
        let file = file.clone();
        super::super::prefetch::try_prefetch(
            move |emit| {
                if file.bruker_streaming {
                    read_bruker_streaming(&path, &file, emit)
                } else {
                    read_bruker_local(&path, &file, emit)
                }
            },
            MS1_BATCH_SIZE,
        )
    }

    const MS1_BATCH_SIZE: usize = 32;

    /// Sorting lightweight metadata preserves the old batch reader's stable RT
    /// order even for files whose acquisition indices are not chronological.
    fn ms1_indices(reader: &FrameReader) -> Result<Vec<usize>, KothError> {
        let mut frames = Vec::new();
        for i in 0..reader.len() {
            let frame = reader
                .get_frame_without_coordinates(i)
                .map_err(|e| KothError::TdfError(e.to_string()))?;
            if frame.ms_level == timsrust::MSLevel::MS1 {
                if !frame.rt_in_seconds.is_finite() {
                    return Err(KothError::TdfError(format!(
                        "non-finite RT for frame {}",
                        frame.index
                    )));
                }
                frames.push((i, frame.rt_in_seconds));
            }
        }
        frames.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        Ok(frames.into_iter().map(|(i, _)| i).collect())
    }

    /// An indexed Rayon map preserves order within each bounded batch. No later
    /// batch is decoded until this one has been delivered to the consumer.
    fn emit_batches<F>(
        indices: &[usize],
        process: F,
        emit: &mut dyn FnMut(Spectrum) -> bool,
    ) -> Result<(), KothError>
    where
        F: Fn(usize) -> Result<Spectrum, KothError> + Sync,
    {
        if indices.is_empty() {
            return Err(KothError::NoSpectra);
        }
        for (batch_index, batch) in indices.chunks(MS1_BATCH_SIZE).enumerate() {
            let spectra: Vec<_> = batch.par_iter().map(|&i| process(i)).collect();
            for (offset, spectrum) in spectra.into_iter().enumerate() {
                let mut spectrum = spectrum?;
                spectrum.scan_index = batch_index * MS1_BATCH_SIZE + offset;
                if !emit(spectrum) {
                    return Ok(());
                }
            }
        }
        Ok(())
    }

    #[cfg(test)]
    mod streaming_tests {
        use super::*;

        fn spectrum(i: usize) -> Spectrum {
            Spectrum {
                scan_index: usize::MAX,
                retention_time: i as f64,
                peaks: Vec::new(),
                ms_level: 1,
                isolation_window: None,
                faims_cv: None,
            }
        }

        #[test]
        fn ordered_batches_stop_decoding_when_consumer_stops() {
            let processed = AtomicUsize::new(0);
            let indices: Vec<_> = (0..1000).collect();
            emit_batches(
                &indices,
                |i| {
                    processed.fetch_add(1, AtomicOrdering::Relaxed);
                    Ok(spectrum(i))
                },
                &mut |s| {
                    assert_eq!(s.scan_index, 0);
                    false
                },
            )
            .unwrap();
            assert_eq!(processed.load(AtomicOrdering::Relaxed), MS1_BATCH_SIZE);
        }

        #[test]
        fn batch_boundaries_preserve_order_and_indices() {
            let indices: Vec<_> = (0..MS1_BATCH_SIZE * 3 + 1).rev().collect();
            let mut seen = Vec::new();
            emit_batches(&indices, |i| Ok(spectrum(i)), &mut |s| {
                assert_eq!(s.scan_index, seen.len());
                seen.push(s.retention_time as usize);
                true
            })
            .unwrap();
            assert_eq!(seen, indices);
        }

        #[test]
        fn decode_error_stops_delivery_and_empty_input_errors() {
            let mut seen = Vec::new();
            let result = emit_batches(
                &[0, 1, 2],
                |i| {
                    if i == 1 {
                        Err(KothError::TdfError("corrupt frame".into()))
                    } else {
                        Ok(spectrum(i))
                    }
                },
                &mut |s| {
                    seen.push(s.scan_index);
                    true
                },
            );
            assert!(matches!(result, Err(KothError::TdfError(_))));
            assert_eq!(seen, vec![0]);
            assert!(matches!(
                emit_batches(&[], |i| Ok(spectrum(i)), &mut |_| true),
                Err(KothError::NoSpectra)
            ));
        }
    }

    /// dnoise vertical-IM filter parameters from the run's `[file]` config. Shared
    /// by every Bruker path (MS1 local, MS1 streaming, diaPASEF MS2) so the values
    /// are defined once.
    fn filter_params(file: &FileConfig) -> FilterParams {
        FilterParams {
            mz_half_width: file.bruker_filter_mz_half_width,
            min_feature_length: file.bruker_filter_min_feature_length,
            max_internal_gap: file.bruker_filter_max_internal_gap,
            min_window_intensity: file.bruker_filter_min_window_intensity,
            min_feature_intensity: file.bruker_filter_min_feature_intensity,
            num_iterations: file.bruker_filter_num_iterations,
        }
    }

    /// dnoise watershed-centroiding parameters from the run's `[file]` config.
    /// Shared by every Bruker path (see [`filter_params`]).
    fn watershed_params(file: &FileConfig) -> WatershedParams {
        WatershedParams {
            box_scan: file.bruker_watershed_box_scan,
            box_mz_idx: file.bruker_watershed_box_mz_idx,
            min_seed_intensity: file.bruker_watershed_min_seed_intensity,
            min_centroid_total: file.bruker_watershed_min_centroid_total,
            max_tof_offset: file.bruker_watershed_max_tof_offset,
        }
    }

    /// In-process streaming path (opt-in via `bruker_streaming`): drive dnoise's
    /// [`RunContext`], which runs the vertical-IM filter + horizontal halo +
    /// watershed in a single pass over each raw MS1 frame — the exact stage code
    /// the standalone `dnoise` tool runs, with no denoised `.d` written to disk.
    /// Calibration comes straight from the context, so no timsrust converters here.
    fn read_bruker_streaming(
        path: &Path,
        file: &FileConfig,
        emit: &mut dyn FnMut(Spectrum) -> bool,
    ) -> Result<(), KothError> {
        let filter_params = filter_params(file);
        let watershed_params = watershed_params(file);
        let halo_params = HaloParams {
            peak_fraction: file.bruker_halo_peak_fraction,
            mz_idx_half_width: file.bruker_halo_mz_idx_half_width,
            scan_half_width: file.bruker_halo_scan_half_width,
        };
        let polygon_params = Ms1PolygonParams {
            mz_pad: file.bruker_ms1_polygon_mz_pad,
            im_pad: file.bruker_ms1_polygon_im_pad,
        };
        // Vertical filter + optional halo + optional MS1 selection-polygon gate +
        // watershed, all in one in-process pass. `Stages` carries more stages than
        // koth wires up (MS/MS denoising, smoothing, the box centroider, and the
        // two diaPASEF window gates); those stay at their dnoise defaults, which
        // are off. Any stage added here needs a `bruker_*` key beside it, or it is
        // unreachable from a config file.
        let stages = Stages {
            halo: file.bruker_halo.then_some(&halo_params),
            ms1_polygon: file.bruker_ms1_polygon.then_some(&polygon_params),
            watershed: Some(&watershed_params),
            ..Stages::default()
        };

        let ctx = RunContext::open(path, &filter_params, &stages)
            .map_err(|e| KothError::TdfError(e.to_string()))?;
        let cal = ctx.calibration();
        let mobility = calibrated_mobility(path, file)?;
        let noise_sigma = file.bruker_noise_sigma;

        let indices =
            ms1_indices(&FrameReader::new(path).map_err(|e| KothError::TdfError(e.to_string()))?)?;
        emit_batches(
            &indices,
            |i| {
                let decoded = ctx
                    .process(i)
                    .map_err(|e| KothError::TdfError(e.to_string()))?;
                let mut peaks: Vec<Peak> = decoded
                    .survivors
                    .iter()
                    .map(|&(scan, tof, intensity)| Peak {
                        mz: cal.tof_to_mz(tof) as f32,
                        intensity: intensity as f32,
                        ion_mobility: match &mobility {
                            Some(m) => m.convert(decoded.frame_id, scan),
                            None => cal.scan_to_im(scan),
                        } as f32,
                    })
                    .collect();
                peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));

                let mut spectrum = Spectrum {
                    scan_index: decoded.frame_id,
                    retention_time: decoded.rt_seconds / 60.0,
                    peaks,
                    ms_level: 1,
                    isolation_window: None,
                    faims_cv: None,
                };
                if let Some(sigma) = noise_sigma {
                    noise::filter_spectrum(&mut spectrum, sigma);
                }
                Ok(spectrum)
            },
            emit,
        )
    }

    /// Historical local path (default): timsrust frame iteration + dnoise's
    /// vertical filter and watershed (no halo), matching the paper's validated
    /// pipeline. Spectrum delivery is bounded in both preprocessing modes.
    fn read_bruker_local(
        path: &Path,
        file: &FileConfig,
        emit: &mut dyn FnMut(Spectrum) -> bool,
    ) -> Result<(), KothError> {
        let path_str = path
            .to_str()
            .ok_or_else(|| KothError::UnsupportedFormat("non-UTF8 path".into()))?;

        let frame_reader =
            FrameReader::new(path_str).map_err(|e| KothError::TdfError(e.to_string()))?;
        let metadata = timsrust::readers::MetadataReader::new(path_str)
            .map_err(|e| KothError::TdfError(e.to_string()))?;

        let mz_converter = metadata.mz_converter;
        let ims_converter = metadata.im_converter;
        let mobility = calibrated_mobility(path, file)?;

        let filter_params = filter_params(file);
        let watershed_params = watershed_params(file);
        let noise_sigma = file.bruker_noise_sigma;

        let processed = AtomicUsize::new(0);
        const PROGRESS_INTERVAL: usize = 100;

        let indices = ms1_indices(&frame_reader)?;
        emit_batches(
            &indices,
            |i| {
                let frame = frame_reader
                    .get(i)
                    .map_err(|e| KothError::TdfError(e.to_string()))?;
                let flat = FlatFrame::from_frame(&frame);
                let n_raw = flat.len();
                let keep = filter_iterated(&flat, &filter_params);
                let survivors = flat.survivors(&keep);
                let n_filtered = survivors.len();

                let centroids = watershed_centroid(&survivors, &watershed_params, MAX_PEAKS);
                let n_centroids = centroids.len();

                let mut out_peaks: Vec<Peak> = centroids
                    .into_iter()
                    .map(|(scan, tof, intensity)| Peak {
                        mz: mz_converter.convert(tof as f64) as f32,
                        intensity: intensity as f32,
                        ion_mobility: match &mobility {
                            Some(m) => m.convert(frame.index, scan),
                            None => ims_converter.convert(scan as f64),
                        } as f32,
                    })
                    .collect();
                out_peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));

                let mut spectrum = Spectrum {
                    scan_index: frame.index,
                    retention_time: frame.rt_in_seconds / 60.0,
                    peaks: out_peaks,
                    ms_level: 1,
                    isolation_window: None,
                    faims_cv: None,
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

                Ok(spectrum)
            },
            emit,
        )
    }

    // -----------------------------------------------------------------------
    // MS2 (diaPASEF) reader
    // -----------------------------------------------------------------------

    /// Map a diaPASEF MS2 frame's [`timsrust::QuadrupoleSettings`] to its fixed
    /// isolation-window segments. The struct stores PARALLEL vectors: segment `i`
    /// covers scan range `[scan_starts[i], scan_ends[i])` (half-open — adjacent
    /// diaPASEF windows share their boundary scan, verified against the
    /// `DiaFrameMsMsWindows` table), isolated at `isolation_mz[i]` with total
    /// width `isolation_width[i]`.
    ///
    /// Returns `(scan_start, scan_end, window)` triples; degenerate rows
    /// (empty scan range, non-finite / non-positive center) are dropped. On a
    /// ddaPASEF / MS1 / unknown frame the settings are empty, so this yields an
    /// empty Vec — the reader below never even calls it in that case because it
    /// gates on `AcquisitionType` first.
    pub(crate) fn window_segments(
        qs: &timsrust::QuadrupoleSettings,
    ) -> Vec<(usize, usize, IsolationWindow)> {
        (0..qs.len())
            .filter_map(|i| {
                let scan_start = *qs.scan_starts.get(i)?;
                let scan_end = *qs.scan_ends.get(i)?;
                let center = *qs.isolation_mz.get(i)?;
                let half = qs.isolation_width.get(i).copied().unwrap_or(0.0) / 2.0;
                if scan_end <= scan_start
                    || !center.is_finite()
                    || center <= 0.0
                    || !half.is_finite()
                    || half < 0.0
                {
                    return None;
                }
                Some((
                    scan_start,
                    scan_end,
                    IsolationWindow {
                        target: center,
                        lower: center - half,
                        upper: center + half,
                    },
                ))
            })
            .collect()
    }

    /// Collapse one isolation-window segment of an MS2 frame down the ion-mobility
    /// axis into centroided `(m/z, 1/K0, intensity)` peaks, using the **same**
    /// dnoise vertical-IM filter + watershed centroiding as the MS1 path.
    ///
    /// The segment is materialized as a sub-[`FlatFrame`] holding only the points
    /// whose scan index is in `[scan_start, scan_end)`. Absolute scan indices and
    /// the parent frame's `num_scans` are preserved, so the filter's per-scan IM
    /// profile and the watershed grid operate in the frame's native scan space
    /// (the extra empty scans outside the segment are inert).
    #[allow(clippy::too_many_arguments)]
    fn segment_peaks(
        flat: &FlatFrame,
        scan_start: usize,
        scan_end: usize,
        filter_params: &FilterParams,
        watershed_params: &WatershedParams,
        mz_converter: &timsrust::converters::Tof2MzConverter,
        ims_converter: &timsrust::converters::Scan2ImConverter,
        mobility: Option<&ScanToMobility>,
    ) -> Vec<Peak> {
        let mut scan = Vec::new();
        let mut tof = Vec::new();
        let mut intensity = Vec::new();
        for i in 0..flat.len() {
            let s = flat.scan[i] as usize;
            if s >= scan_start && s < scan_end {
                scan.push(flat.scan[i]);
                tof.push(flat.tof[i]);
                intensity.push(flat.intensity[i]);
            }
        }
        let seg = FlatFrame {
            frame_id: flat.frame_id,
            num_scans: flat.num_scans,
            scan,
            tof,
            intensity,
        };
        if seg.is_empty() {
            return Vec::new();
        }

        let keep = filter_iterated(&seg, filter_params);
        let survivors = seg.survivors(&keep);
        let centroids = watershed_centroid(&survivors, watershed_params, MAX_PEAKS);

        let mut peaks: Vec<Peak> = centroids
            .into_iter()
            .map(|(scan, tof, intensity)| Peak {
                mz: mz_converter.convert(tof as f64) as f32,
                intensity: intensity as f32,
                ion_mobility: match mobility {
                    Some(m) => m.convert(flat.frame_id, scan),
                    None => ims_converter.convert(scan as f64),
                } as f32,
            })
            .collect();
        peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));
        peaks
    }

    /// diaPASEF MS2 reader. Iterates the `.d`'s MS2 frames; for each frame, for
    /// each fixed isolation-window segment, runs the same IM-collapse as MS1 on
    /// just that segment's peaks and emits one MS2 [`Spectrum`] carrying the
    /// segment's [`IsolationWindow`]. One MS2 frame -> one Spectrum per window
    /// segment (mirroring how mzML yields one MS2 spectrum per isolation window).
    ///
    /// **diaPASEF only.** If the acquisition is ddaPASEF (per-precursor selection)
    /// or unknown, this emits an empty Vec plus a warning and does not touch MS1.
    /// The MS1 Bruker path and every non-`.d` path are completely unaffected; this
    /// only runs when the caller opts into MS2 on a `.d`.
    pub fn read_bruker_ms2(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
        let path_str = path
            .to_str()
            .ok_or_else(|| KothError::UnsupportedFormat("non-UTF8 path".into()))?;

        let frame_reader =
            FrameReader::new(path_str).map_err(|e| KothError::TdfError(e.to_string()))?;

        // diaPASEF-only guard. ddaPASEF MS2 frames describe per-precursor
        // selection (driven by a precursor table), not fixed windows, so we do
        // not try to reconstruct them here — emit nothing and leave MS1 alone.
        let acquisition = frame_reader.get_acquisition();
        if acquisition != timsrust::AcquisitionType::DIAPASEF {
            log::warn!(
                "Bruker MS2 hill detection supports diaPASEF only; '{}' is {:?} — \
                 emitting no MS2 spectra (MS1 output is unaffected)",
                path.display(),
                acquisition,
            );
            return Ok(Vec::new());
        }

        let metadata = timsrust::readers::MetadataReader::new(path_str)
            .map_err(|e| KothError::TdfError(e.to_string()))?;
        let mz_converter = metadata.mz_converter;
        let ims_converter = metadata.im_converter;
        let mobility = calibrated_mobility(path, file)?;

        let filter_params = filter_params(file);
        let watershed_params = watershed_params(file);
        let noise_sigma = file.bruker_noise_sigma;

        // One Vec<Spectrum> per MS2 frame, in ascending frame (== RT) order so
        // the flattened result is deterministic (rayon collect preserves the
        // input index order).
        let per_frame: Vec<Vec<Spectrum>> = frame_reader
            .parallel_filter(|f| f.ms_level == timsrust::MSLevel::MS2)
            .map(|frame_result| -> Vec<Spectrum> {
                let frame = match frame_result {
                    Ok(f) => f,
                    Err(e) => {
                        log::error!("Error parsing Bruker MS2 frame: {:?}", e);
                        return Vec::new();
                    }
                };
                let segments = window_segments(&frame.quadrupole_settings);
                if segments.is_empty() {
                    return Vec::new();
                }
                let flat = FlatFrame::from_frame(&frame);
                let rt_minutes = frame.rt_in_seconds / 60.0;

                segments
                    .into_iter()
                    .filter_map(|(scan_start, scan_end, window)| {
                        let peaks = segment_peaks(
                            &flat,
                            scan_start,
                            scan_end,
                            &filter_params,
                            &watershed_params,
                            &mz_converter,
                            &ims_converter,
                            mobility.as_ref(),
                        );
                        let mut spectrum = Spectrum {
                            // Source frame id, for traceability. The MS2 hill
                            // detector re-indexes scans per isolation window, so
                            // this value is not used for gap tracking.
                            scan_index: frame.index,
                            retention_time: rt_minutes,
                            peaks,
                            ms_level: 2,
                            isolation_window: Some(window),
                            faims_cv: None,
                        };
                        if let Some(sigma) = noise_sigma {
                            noise::filter_spectrum(&mut spectrum, sigma);
                        }
                        if spectrum.peaks.is_empty() {
                            None
                        } else {
                            Some(spectrum)
                        }
                    })
                    .collect()
            })
            .collect();

        let mut spectra: Vec<Spectrum> = per_frame.into_iter().flatten().collect();
        // Ascending RT so each per-window detector sees its cycles in acquisition
        // order. Already in frame order from the parallel collect; this is a
        // deterministic safety net (stable sort keeps a frame's segments together).
        spectra.sort_by(|a, b| {
            a.retention_time
                .partial_cmp(&b.retention_time)
                .unwrap_or(Ordering::Equal)
        });

        let n_windows = spectra
            .iter()
            .filter_map(|s| s.isolation_window.map(|w| w.key()))
            .collect::<std::collections::BTreeSet<_>>()
            .len();
        log::info!(
            "Read {} MS2 spectra across {} isolation windows from Bruker .d (diaPASEF)",
            spectra.len(),
            n_windows,
        );
        Ok(spectra)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn qs(
            scan_starts: Vec<usize>,
            scan_ends: Vec<usize>,
            isolation_mz: Vec<f64>,
            isolation_width: Vec<f64>,
        ) -> timsrust::QuadrupoleSettings {
            let n = isolation_mz.len();
            timsrust::QuadrupoleSettings {
                index: 1,
                scan_starts,
                scan_ends,
                isolation_mz,
                isolation_width,
                collision_energy: vec![0.0; n],
            }
        }

        #[test]
        fn window_segments_maps_diapasef_group_to_windows() {
            // A real 5-min diaPASEF window group (3 windows, 25 Th wide, half-open
            // scan ranges that share boundaries: 530, 728).
            let settings = qs(
                vec![125, 530, 728],
                vec![530, 728, 935],
                vec![812.5, 612.5, 412.5],
                vec![25.0, 25.0, 25.0],
            );
            let segs = window_segments(&settings);
            assert_eq!(segs.len(), 3);

            assert_eq!(segs[0].0, 125);
            assert_eq!(segs[0].1, 530);
            assert!((segs[0].2.target - 812.5).abs() < 1e-9);
            assert!((segs[0].2.lower - 800.0).abs() < 1e-9);
            assert!((segs[0].2.upper - 825.0).abs() < 1e-9);

            // Boundary scan 530 is the exclusive end of window 0 and the inclusive
            // start of window 1 — half-open, so it is not double-owned.
            assert_eq!(segs[1].0, 530);
            assert_eq!(segs[1].1, 728);
            assert!((segs[1].2.target - 612.5).abs() < 1e-9);

            assert_eq!(segs[2].0, 728);
            assert_eq!(segs[2].1, 935);
            assert!((segs[2].2.lower - 400.0).abs() < 1e-9);
            assert!((segs[2].2.upper - 425.0).abs() < 1e-9);
        }

        #[test]
        fn window_segments_empty_for_dda_or_unset_settings() {
            // ddaPASEF / MS1 / unknown frames carry Default (empty) settings.
            let settings = timsrust::QuadrupoleSettings::default();
            assert!(window_segments(&settings).is_empty());
        }

        #[test]
        fn window_segments_drops_degenerate_rows() {
            // Empty scan range and a non-positive center are both rejected; the
            // one valid row survives.
            let settings = qs(
                vec![100, 300, 500],
                vec![100, 400, 500], // row0 start==end, row2 start==end
                vec![500.0, 600.0, -1.0],
                vec![25.0, 25.0, 25.0],
            );
            let segs = window_segments(&settings);
            assert_eq!(segs.len(), 1);
            assert_eq!(segs[0].0, 300);
            assert_eq!(segs[0].1, 400);
            assert!((segs[0].2.target - 600.0).abs() < 1e-9);
        }
    }
}

#[cfg(feature = "tdf")]
pub use inner::{read_bruker, read_bruker_ms2, stream_bruker};
