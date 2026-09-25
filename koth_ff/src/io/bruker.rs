/// Bruker .d folder reader. Two per-frame paths, selected by `bruker_streaming`:
///
/// * **local** (default, `read_bruker_local`): timsrust frame iteration, then
///   dnoise's `filter_iterated` (vertical-IM) + `watershed_centroid`, converting
///   each centroid to `(m/z, 1/K0, intensity)` via the `dnoise::tsr` converters
///   (timsrust 0.4.2's lines). This is the paper's validated pipeline (no halo).
/// * **streaming** (opt-in, `read_bruker_streaming`): dnoise's `RunContext` runs
///   the vertical filter + horizontal halo + watershed in a single in-process
///   pass over the raw frames — the exact stage code the standalone `dnoise` tool
///   runs — with no denoised `.d` written to disk. Calibration comes from the
///   context. The flag selects preprocessing; both modes stream into detection.
///
/// Both paths then optionally apply the `bruker_noise_sigma` MAD filter, emit in
/// retention-time order, and assign sequential scan indices. All
/// denoising / centroiding logic lives in the `dnoise` crate. Frames are decoded
/// with `dnoise::tsr::FrameReader` (timsrust-tdf 0.6 behind timsrust 0.4.2's
/// frame order); frame RT, MS level and diaPASEF windows come from
/// `analysis.tdf` directly (`inner::frame_table`), read the way timsrust 0.4.2 did.
#[cfg(feature = "tdf")]
pub mod inner {
    use std::cmp::Ordering;
    use std::path::Path;
    use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

    use rayon::prelude::*;

    use dnoise::tsr::{ConvertableDomain, FrameReader, MetadataReader};
    use dnoise::{
        filter_iterated, watershed_centroid, FilterParams, FlatFrame, HaloParams, Ms1PolygonParams,
        RunContext, Stages, WatershedParams,
    };

    use crate::config::FileConfig;
    use crate::error::KothError;
    use crate::hills::noise;
    use crate::io::tims_calibration::{self, MobilityScale, ScanToMobility};
    use crate::models::{IsolationWindow, Peak, Spectrum};

    use frame_table::{Acquisition, FrameRow, FrameTable, MsLevel, QuadrupoleSettings};

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
    /// Indices are 0-based positions in `Frames.Id` order, as the readers take.
    fn ms1_indices(table: &FrameTable) -> Result<Vec<usize>, KothError> {
        let mut frames = Vec::new();
        for (i, frame) in table.frames.iter().enumerate() {
            if frame.ms_level == MsLevel::Ms1 {
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

        fn polygon_tuple(p: &Ms1PolygonParams) -> (f64, f64, bool) {
            (p.mz_pad, p.im_pad, p.overlap)
        }

        #[test]
        fn polygon_defaults_match_point_by_point_gating() {
            // The parameters koth passed to dnoise before the overlap keys
            // existed. Defaults, and a config that omits the keys, must give
            // exactly these, so existing outputs are unchanged; in particular
            // dnoise 0.5's padded, feature-level defaults must not leak in.
            let before = (0.0, 0.0, false);
            assert_eq!(
                polygon_tuple(&ms1_polygon_params(&FileConfig::default())),
                before
            );
            let parsed: FileConfig = toml::from_str("bruker_ms1_polygon = true").unwrap();
            assert!(!parsed.bruker_ms1_polygon_overlap);
            assert_eq!(polygon_tuple(&ms1_polygon_params(&parsed)), before);
            // The deprecated reach still parses and changes nothing.
            let reach_only: FileConfig =
                toml::from_str("bruker_ms1_polygon_overlap_reach = 0.3").unwrap();
            assert_eq!(polygon_tuple(&ms1_polygon_params(&reach_only)), before);
        }

        #[test]
        fn polygon_overlap_opt_in_maps_onto_dnoise() {
            let on: FileConfig = toml::from_str(
                "bruker_ms1_polygon_overlap = true\nbruker_ms1_polygon_mz_pad = 5.0",
            )
            .unwrap();
            assert_eq!(polygon_tuple(&ms1_polygon_params(&on)), (5.0, 0.0, true));
        }

        #[test]
        fn halo_default_pins_dnoise_040_peak_fraction() {
            // dnoise 0.5 changed its own default to 0.10; koth passes its value
            // explicitly, and its default stays 0.4.0's 0.15.
            assert_eq!(FileConfig::default().bruker_halo_peak_fraction, 0.15);
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

    /// dnoise MS1 selection-polygon parameters from the run's `[file]` config.
    /// Every field is set from koth's config, so dnoise's own defaults (0.5:
    /// padded 3.0 Th / 0.015 1/K0, feature-level gating) never apply: koth
    /// gates point by point on the literal polygon unless configured otherwise.
    /// `bruker_ms1_polygon_overlap` opts in to feature-level gating.
    fn ms1_polygon_params(file: &FileConfig) -> Ms1PolygonParams {
        Ms1PolygonParams {
            mz_pad: file.bruker_ms1_polygon_mz_pad,
            im_pad: file.bruker_ms1_polygon_im_pad,
            overlap: file.bruker_ms1_polygon_overlap,
        }
    }

    /// dnoise 0.5 removed the overlap reach: a feature kept by feature-level
    /// gating is kept over its whole extent. Warn when a config still asks for
    /// a finite reach, since its output differs from dnoise 0.4's.
    fn warn_ignored_overlap_reach(file: &FileConfig) {
        if file.bruker_ms1_polygon
            && file.bruker_ms1_polygon_overlap
            && file.bruker_ms1_polygon_overlap_reach != 0.0
        {
            log::warn!(
                "bruker_ms1_polygon_overlap_reach = {} is ignored: dnoise 0.5 keeps each \
                 feature kept by the polygon gate over its whole extent",
                file.bruker_ms1_polygon_overlap_reach
            );
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
        let polygon_params = ms1_polygon_params(file);
        warn_ignored_overlap_reach(file);
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
            // The polygon gate converts scans on the same 1/K0 scale as koth.
            mobility_scale: crate::io::tims_calibration::dnoise_scale(file.bruker_mobility_scale),
            ..Stages::default()
        };

        let ctx = RunContext::open(path, &filter_params, &stages)
            .map_err(|e| KothError::TdfError(e.to_string()))?;
        let cal = ctx.calibration();
        let mobility = calibrated_mobility(path, file)?;
        let noise_sigma = file.bruker_noise_sigma;

        let indices = ms1_indices(&FrameTable::read(path)?)?;
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
        let metadata =
            MetadataReader::new(path_str).map_err(|e| KothError::TdfError(e.to_string()))?;
        let table = FrameTable::read(path)?;

        let mz_converter = metadata.mz_converter;
        let ims_converter = metadata.im_converter;
        let mobility = calibrated_mobility(path, file)?;

        let filter_params = filter_params(file);
        let watershed_params = watershed_params(file);
        let noise_sigma = file.bruker_noise_sigma;

        let processed = AtomicUsize::new(0);
        const PROGRESS_INTERVAL: usize = 100;

        let indices = ms1_indices(&table)?;
        emit_batches(
            &indices,
            |i| {
                let frame = frame_reader
                    .get(i)
                    .map_err(|e| KothError::TdfError(e.to_string()))?;
                let rt_in_seconds = table.frames[i].rt_in_seconds;
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
                    retention_time: rt_in_seconds / 60.0,
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

    /// Map a diaPASEF MS2 frame's [`QuadrupoleSettings`] to its fixed
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
    pub(crate) fn window_segments(qs: &QuadrupoleSettings) -> Vec<(usize, usize, IsolationWindow)> {
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
        mz_converter: &dnoise::tsr::Tof2MzConverter,
        ims_converter: &dnoise::tsr::Scan2ImConverter,
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

        let table = FrameTable::read(path)?;

        // diaPASEF-only guard. ddaPASEF MS2 frames describe per-precursor
        // selection (driven by a precursor table), not fixed windows, so we do
        // not try to reconstruct them here — emit nothing and leave MS1 alone.
        let acquisition = table.acquisition;
        if acquisition != Acquisition::DiaPasef {
            log::warn!(
                "Bruker MS2 hill detection supports diaPASEF only; '{}' is {:?} — \
                 emitting no MS2 spectra (MS1 output is unaffected)",
                path.display(),
                acquisition,
            );
            return Ok(Vec::new());
        }

        let frame_reader =
            FrameReader::new(path_str).map_err(|e| KothError::TdfError(e.to_string()))?;
        let metadata =
            MetadataReader::new(path_str).map_err(|e| KothError::TdfError(e.to_string()))?;
        let mz_converter = metadata.mz_converter;
        let ims_converter = metadata.im_converter;
        let mobility = calibrated_mobility(path, file)?;

        let filter_params = filter_params(file);
        let watershed_params = watershed_params(file);
        let noise_sigma = file.bruker_noise_sigma;

        // One Vec<Spectrum> per MS2 frame, in ascending frame (== RT) order so
        // the flattened result is deterministic (rayon collect preserves the
        // input index order).
        let ms2: Vec<(usize, &FrameRow)> = table
            .frames
            .iter()
            .enumerate()
            .filter(|(_, f)| f.ms_level == MsLevel::Ms2)
            .collect();
        let per_frame: Vec<Vec<Spectrum>> = ms2
            .into_par_iter()
            .map(|(i, row)| -> Vec<Spectrum> {
                let segments = window_segments(table.quadrupole_settings(row));
                if segments.is_empty() {
                    return Vec::new();
                }
                let frame = match frame_reader.get(i) {
                    Ok(f) => f,
                    Err(e) => {
                        log::error!("Error parsing Bruker MS2 frame: {:?}", e);
                        return Vec::new();
                    }
                };
                let flat = FlatFrame::from_frame(&frame);
                let rt_minutes = row.rt_in_seconds / 60.0;

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

    /// Frame metadata straight from `analysis.tdf`, read the way timsrust 0.4.2's
    /// `FrameReader` did: rows in `Frames.Id` order (0.4.2's table order, and the
    /// order `dnoise::tsr::FrameReader::get` takes), `MsMsType` 0 = MS1, 8/9 =
    /// MS2, acquisition ddaPASEF if any frame is type 8, else diaPASEF if any is
    /// type 9, and diaPASEF isolation windows grouped by `WindowGroup` with their
    /// rows sorted by `ScanNumBegin`. Unreadable cells default to zero, as 0.4.2's
    /// `parse_default` did.
    pub(crate) mod frame_table {
        use std::path::Path;
        use std::sync::Arc;

        use rusqlite::types::FromSql;
        use rusqlite::{Connection, OpenFlags, Row};

        use crate::error::KothError;

        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub(crate) enum MsLevel {
            Ms1,
            Ms2,
            Unknown,
        }

        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub(crate) enum Acquisition {
            DdaPasef,
            DiaPasef,
            Unknown,
        }

        /// One diaPASEF window group: parallel vectors, segment `i` covers scans
        /// `[scan_starts[i], scan_ends[i])` isolated at `isolation_mz[i]` with
        /// total width `isolation_width[i]` (timsrust 0.4's layout).
        // `index` and `collision_energy` are kept for timsrust 0.4's layout and
        // for debugging; the reader uses the scan ranges and windows.
        #[allow(dead_code)]
        #[derive(Debug, Clone, Default, PartialEq)]
        pub(crate) struct QuadrupoleSettings {
            pub index: usize,
            pub scan_starts: Vec<usize>,
            pub scan_ends: Vec<usize>,
            pub isolation_mz: Vec<f64>,
            pub isolation_width: Vec<f64>,
            pub collision_energy: Vec<f64>,
        }

        impl QuadrupoleSettings {
            pub fn len(&self) -> usize {
                self.isolation_mz.len()
            }
        }

        #[derive(Debug, Clone)]
        pub(crate) struct FrameRow {
            /// `Frames.Id`.
            pub index: usize,
            pub rt_in_seconds: f64,
            pub ms_level: MsLevel,
            /// diaPASEF window group (1-based) of an MS2 frame, else `None`.
            pub window_group: Option<usize>,
        }

        #[derive(Debug, Clone)]
        pub(crate) struct FrameTable {
            pub frames: Vec<FrameRow>,
            pub acquisition: Acquisition,
            groups: Vec<Arc<QuadrupoleSettings>>,
            empty: Arc<QuadrupoleSettings>,
        }

        fn get<T: Default + FromSql>(row: &Row, i: usize) -> T {
            row.get(i).unwrap_or_default()
        }

        impl FrameTable {
            pub fn read(path: &Path) -> Result<Self, KothError> {
                let tdf = path.join("analysis.tdf");
                let err =
                    |e: rusqlite::Error| KothError::TdfError(format!("{}: {e}", tdf.display()));
                let con = Connection::open_with_flags(&tdf, OpenFlags::SQLITE_OPEN_READ_ONLY)
                    .map_err(err)?;
                let mut stmt = con
                    .prepare("SELECT Id, MsMsType, Time FROM Frames ORDER BY Id")
                    .map_err(err)?;
                let raw: Vec<(usize, u8, f64)> = stmt
                    .query_map([], |r| {
                        Ok((get::<i64>(r, 0) as usize, get::<u8>(r, 1), get::<f64>(r, 2)))
                    })
                    .map_err(err)?
                    .collect::<Result<_, _>>()
                    .map_err(err)?;
                let acquisition = if raw.iter().any(|f| f.1 == 8) {
                    Acquisition::DdaPasef
                } else if raw.iter().any(|f| f.1 == 9) {
                    Acquisition::DiaPasef
                } else {
                    Acquisition::Unknown
                };

                let mut groups = Vec::new();
                let mut frame_group = std::collections::HashMap::new();
                if acquisition == Acquisition::DiaPasef {
                    let mut stmt = con
                        .prepare("SELECT Frame, WindowGroup FROM DiaFrameMsMsInfo")
                        .map_err(err)?;
                    for row in stmt
                        .query_map([], |r| {
                            Ok((get::<i64>(r, 0) as usize, get::<i64>(r, 1) as usize))
                        })
                        .map_err(err)?
                    {
                        let (frame, group) = row.map_err(err)?;
                        frame_group.insert(frame, group);
                    }
                    let mut stmt = con
                        .prepare(
                            "SELECT WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, \
                             IsolationWidth, CollisionEnergy FROM DiaFrameMsMsWindows",
                        )
                        .map_err(err)?;
                    let rows: Vec<(usize, usize, usize, f64, f64, f64)> = stmt
                        .query_map([], |r| {
                            Ok((
                                get::<i64>(r, 0) as usize,
                                get::<i64>(r, 1) as usize,
                                get::<i64>(r, 2) as usize,
                                get::<f64>(r, 3),
                                get::<f64>(r, 4),
                                get::<f64>(r, 5),
                            ))
                        })
                        .map_err(err)?
                        .collect::<Result<_, _>>()
                        .map_err(err)?;
                    let n_groups = rows.iter().map(|r| r.0).max().unwrap_or(0);
                    let mut settings: Vec<QuadrupoleSettings> = (0..n_groups)
                        .map(|g| QuadrupoleSettings {
                            index: g + 1,
                            ..Default::default()
                        })
                        .collect();
                    for (group, start, end, mz, width, ce) in rows {
                        if group == 0 {
                            continue;
                        }
                        let q = &mut settings[group - 1];
                        q.scan_starts.push(start);
                        q.scan_ends.push(end);
                        q.isolation_mz.push(mz);
                        q.isolation_width.push(width);
                        q.collision_energy.push(ce);
                    }
                    // Stable sort of each group's rows by scan start, as 0.4.2's
                    // `argsort` (a stable `sort_by_key`) did.
                    for q in &mut settings {
                        let mut order: Vec<usize> = (0..q.scan_starts.len()).collect();
                        order.sort_by_key(|&i| q.scan_starts[i]);
                        let pick = |v: &Vec<f64>| order.iter().map(|&i| v[i]).collect::<Vec<_>>();
                        let pick_u =
                            |v: &Vec<usize>| order.iter().map(|&i| v[i]).collect::<Vec<_>>();
                        *q = QuadrupoleSettings {
                            index: q.index,
                            scan_starts: pick_u(&q.scan_starts),
                            scan_ends: pick_u(&q.scan_ends),
                            isolation_mz: pick(&q.isolation_mz),
                            isolation_width: pick(&q.isolation_width),
                            collision_energy: pick(&q.collision_energy),
                        };
                    }
                    groups = settings.into_iter().map(Arc::new).collect();
                }

                let frames = raw
                    .into_iter()
                    .map(|(index, msms_type, rt)| {
                        let ms_level = match msms_type {
                            0 => MsLevel::Ms1,
                            8 | 9 => MsLevel::Ms2,
                            _ => MsLevel::Unknown,
                        };
                        let window_group = (acquisition == Acquisition::DiaPasef
                            && ms_level == MsLevel::Ms2)
                            .then(|| frame_group.get(&index).copied())
                            .flatten();
                        FrameRow {
                            index,
                            rt_in_seconds: rt,
                            ms_level,
                            window_group,
                        }
                    })
                    .collect();
                Ok(Self {
                    frames,
                    acquisition,
                    groups,
                    empty: Arc::new(QuadrupoleSettings::default()),
                })
            }

            /// The isolation windows of `row` (empty unless a diaPASEF MS2 frame).
            pub fn quadrupole_settings(&self, row: &FrameRow) -> &QuadrupoleSettings {
                row.window_group
                    .and_then(|g| g.checked_sub(1))
                    .and_then(|g| self.groups.get(g))
                    .unwrap_or(&self.empty)
            }
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn qs(
            scan_starts: Vec<usize>,
            scan_ends: Vec<usize>,
            isolation_mz: Vec<f64>,
            isolation_width: Vec<f64>,
        ) -> QuadrupoleSettings {
            let n = isolation_mz.len();
            QuadrupoleSettings {
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
            let settings = QuadrupoleSettings::default();
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
