/// Bruker .d folder reader using timsrust.
///
/// Ported from ms1search_rust/src/tdf_reader.rs with adaptations for
/// the koth_ff data model (f64 instead of f32, Spectrum/Peak types).
#[cfg(feature = "tdf")]
pub mod inner {
    use std::cmp::Ordering;
    use std::path::Path;

    use rayon::prelude::*;
    use timsrust::converters::{ConvertableDomain, Scan2ImConverter, Tof2MzConverter};
    use timsrust::readers::FrameReader;

    use crate::error::KothError;
    use crate::models::{Peak, Spectrum};

    const MAX_PEAKS: usize = 10_000;

    pub fn read_bruker(path: &Path, mz_ppm: f64, im_pct: f64) -> Result<Vec<Spectrum>, KothError> {
        let path_str = path.to_str().ok_or_else(|| KothError::UnsupportedFormat("non-UTF8 path".into()))?;

        let frame_reader = FrameReader::new(path_str)
            .map_err(|e| KothError::TdfError(e.to_string()))?;
        let metadata = timsrust::readers::MetadataReader::new(path_str)
            .map_err(|e| KothError::TdfError(e.to_string()))?;

        let mz_converter = metadata.mz_converter;
        let ims_converter = metadata.im_converter;

        let mz_ppm_f32 = mz_ppm as f32;
        let im_pct_f32 = im_pct as f32;

        let mut spectra: Vec<Spectrum> = frame_reader
            .parallel_filter(|f| f.ms_level == timsrust::MSLevel::MS1)
            .filter_map(|frame_result| match frame_result {
                Ok(frame) => {
                    let mut buffer = PeakBuffer::with_capacity(2 * MAX_PEAKS);
                    buffer.with_frame(&frame, &ims_converter, &mz_converter);
                    let (mzs, (intensities, mobilities)) =
                        buffer.fastcentroid_frame(mz_ppm_f32, im_pct_f32);

                    let retention_time = frame.rt_in_seconds as f64 / 60.0;
                    let scan_index = frame.index;

                    let mut peaks: Vec<Peak> = mzs
                        .into_iter()
                        .zip(intensities.into_iter())
                        .zip(mobilities.into_iter())
                        .map(|((mz, intensity), im)| Peak {
                            mz,
                            intensity,
                            ion_mobility: im,
                        })
                        .collect();
                    peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());

                    Some(Spectrum {
                        scan_index,
                        retention_time,
                        peaks,
                    })
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

        // Sort by retention time and re-index
        spectra.sort_by(|a, b| a.retention_time.partial_cmp(&b.retention_time).unwrap());
        for (i, s) in spectra.iter_mut().enumerate() {
            s.scan_index = i;
        }

        log::info!("Read {} MS1 spectra from Bruker .d", spectra.len());
        Ok(spectra)
    }

    #[derive(Clone, Copy)]
    struct ImsPeak {
        mz: f32,
        intensity: f32,
        im: f32,
    }

    #[derive(Clone)]
    struct PeakBuffer {
        peaks: Vec<ImsPeak>,
        order: Vec<usize>,
        agg_buff: Vec<ImsPeak>,
    }

    impl PeakBuffer {
        fn with_capacity(capacity: usize) -> Self {
            Self {
                peaks: Vec::with_capacity(capacity),
                order: Vec::with_capacity(capacity),
                agg_buff: Vec::with_capacity(MAX_PEAKS),
            }
        }

        fn clear(&mut self) {
            self.peaks.clear();
            self.order.clear();
            self.agg_buff.clear();
        }

        fn len(&self) -> usize {
            self.peaks.len()
        }

        fn expand_to_capacity(&mut self, capacity: usize) {
            if capacity <= self.len() {
                return;
            }
            let diff = capacity - self.len();
            let diff = diff.max(self.len() / 5);
            self.peaks.reserve(diff);
            self.order.reserve(diff);
            self.agg_buff.reserve(capacity);
        }

        fn with_frame(
            &mut self,
            frame: &timsrust::Frame,
            ims_converter: &Scan2ImConverter,
            mz_converter: &Tof2MzConverter,
        ) {
            let expect_len = frame.tof_indices.len();
            self.expand_to_capacity(expect_len);

            let mz_iter = frame
                .tof_indices
                .iter()
                .map(|&x| mz_converter.convert(x as f64) as f32);
            let int_iter = frame.intensities.iter().map(|&x| x as f32);
            let im_iter = Self::expand_mobility_iter(&frame.scan_offsets, ims_converter);

            self.peaks.extend(
                mz_iter
                    .zip(int_iter)
                    .zip(im_iter)
                    .map(|((mz, intensity), im)| ImsPeak { mz, intensity, im }),
            );

            self.peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));

            self.order.extend(0..self.len());
            self.order.sort_unstable_by(|&a, &b| {
                self.peaks[b]
                    .intensity
                    .partial_cmp(&self.peaks[a].intensity)
                    .unwrap_or(Ordering::Equal)
            });
        }

        /// Expand the scan_offsets run-length encoding to per-peak ion mobility values.
        fn expand_mobility_iter<'a>(
            scan_offsets: &'a [usize],
            ims_converter: &'a Scan2ImConverter,
        ) -> impl Iterator<Item = f32> + 'a {
            scan_offsets
                .windows(2)
                .enumerate()
                .filter_map(|(i, w)| {
                    let num = w[1] - w[0];
                    if num == 0 {
                        return None;
                    }
                    let im = ims_converter.convert(i as f64) as f32;
                    Some((im, w[0], w[1]))
                })
                .flat_map(|(im, lo, hi)| (lo..hi).map(move |_| im))
        }

        /// Greedy intensity-ordered centroiding within mz_ppm and im_pct tolerances.
        fn fastcentroid_frame(
            &mut self,
            mz_tol_ppm: f32,
            im_tol_pct: f32,
        ) -> (Vec<f32>, (Vec<f32>, Vec<f32>)) {
            debug_assert!(
                self.peaks.windows(2).all(|x| x[0].mz <= x[1].mz),
                "peaks must be mz-sorted"
            );

            let utol = mz_tol_ppm / 1e6;
            let im_tol = im_tol_pct / 100.0;
            let mut global_included = 0;

            for &idx in &self.order {
                if self.peaks[idx].intensity <= 0.0 {
                    continue;
                }
                if self.agg_buff.len() >= MAX_PEAKS {
                    break;
                }

                let mz = self.peaks[idx].mz;
                let im = self.peaks[idx].im;
                let da_tol = mz * utol;
                let abs_im_tol = im * im_tol;

                let ss_start = self.peaks.partition_point(|&x| x.mz < mz - da_tol);
                let ss_end = self.peaks.partition_point(|&x| x.mz <= mz + da_tol);

                let mut curr_intensity = 0.0f32;
                let mut count = 0usize;
                for i in ss_start..ss_end {
                    let p = &mut self.peaks[i];
                    if p.intensity > 0.0
                        && p.im >= im - abs_im_tol
                        && p.im <= im + abs_im_tol
                    {
                        curr_intensity += p.intensity;
                        p.intensity = -1.0;
                        count += 1;
                    }
                }

                if count == 0 {
                    continue;
                }

                self.agg_buff.push(ImsPeak {
                    mz,
                    intensity: curr_intensity,
                    im,
                });
                global_included += count;

                if global_included >= self.len() {
                    break;
                }
            }

            self.agg_buff.sort_unstable_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));

            self.agg_buff
                .drain(..)
                .map(|x| (x.mz, (x.intensity, x.im)))
                .unzip()
        }
    }
}

#[cfg(feature = "tdf")]
pub use inner::read_bruker;
