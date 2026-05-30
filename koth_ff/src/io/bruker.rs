/// Bruker .d folder reader using timsrust.
///
/// Two-stage extraction:
///   1. Raw peaks are flattened from the frame in TOF-index / scan-index
///      space, then optionally smoothed along each axis (rolling average).
///      Index-space smoothing is cheap and amplifies real signal (which
///      clusters along both axes) while dispersing noise. Conversion to
///      physical units (m/z, 1/K0) happens once, after smoothing.
///   2. Greedy intensity-ordered centroiding, followed by a per-centroid
///      satellite-suppression pass that kills sub-threshold "shadow" peaks
///      sitting next to real peaks in a wider m/z window.
#[cfg(feature = "tdf")]
pub mod inner {
    use std::cmp::Ordering;
    use std::collections::{BinaryHeap, HashMap};
    use std::hash::{BuildHasher, Hasher};
    use std::path::Path;
    use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

    use rayon::prelude::*;
    use timsrust::converters::{ConvertableDomain, Scan2ImConverter, Tof2MzConverter};
    use timsrust::readers::FrameReader;

    /// Hasher for `u64` composite keys built from packed `(primary, secondary)`
    /// integer coordinates. SipHash is overkill here (and ~10× too slow); the
    /// keys are well-spread across the u64 space but their low bits can have
    /// poor entropy (e.g. when one axis has only ~900 distinct values), so a
    /// SplitMix64 finalizer is applied to scramble bits into bucket positions.
    #[derive(Default, Clone, Copy)]
    struct U64Hasher(u64);

    impl Hasher for U64Hasher {
        fn finish(&self) -> u64 {
            self.0
        }
        fn write(&mut self, _: &[u8]) {
            unreachable!("U64Hasher is only used with u64 keys")
        }
        fn write_u64(&mut self, n: u64) {
            // SplitMix64 finalizer
            let mut x = n;
            x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            x ^= x >> 31;
            self.0 = x;
        }
    }

    #[derive(Default, Clone, Copy)]
    struct BuildU64Hasher;

    impl BuildHasher for BuildU64Hasher {
        type Hasher = U64Hasher;
        fn build_hasher(&self) -> U64Hasher {
            U64Hasher(0)
        }
    }

    type KeyMap = HashMap<u64, f32, BuildU64Hasher>;

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

        let mz_ppm = file.bruker_mz_ppm as f32;
        let im_pct = file.bruker_im_pct as f32;
        let min_subpeaks = file.bruker_min_subpeaks;
        let im_smoothing = file.bruker_im_smoothing_window as u32;
        let mz_smoothing = file.bruker_mz_smoothing_window as u32;
        let sat_window = file.bruker_satellite_window_da as f32;
        let sat_end_frac = file.bruker_satellite_end_fraction as f32;
        let noise_sigma = file.bruker_noise_sigma;

        let processed = AtomicUsize::new(0);
        const PROGRESS_INTERVAL: usize = 100;

        let mut spectra: Vec<Spectrum> = frame_reader
            .parallel_filter(|f| f.ms_level == timsrust::MSLevel::MS1)
            .filter_map(|frame_result| match frame_result {
                Ok(frame) => {
                    let raw_peak_count = frame.tof_indices.len();

                    let cloud = RawCloud::from_frame(&frame);
                    let cloud = smooth_scan_axis(cloud, im_smoothing);
                    let cloud = smooth_tof_axis(cloud, mz_smoothing);

                    let mut peaks = cloud.into_ims_peaks(&mz_converter, &ims_converter);
                    let n_centroids;
                    let out_peaks: Vec<Peak> = if peaks.is_empty() {
                        n_centroids = 0;
                        Vec::new()
                    } else {
                        peaks.sort_by(|a, b| {
                            a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal)
                        });
                        let heap = build_intensity_heap(&peaks);

                        let centroids = centroid_with_satellite(
                            &mut peaks,
                            heap,
                            mz_ppm,
                            im_pct,
                            min_subpeaks,
                            sat_window,
                            sat_end_frac,
                            MAX_PEAKS,
                        );
                        n_centroids = centroids.len();

                        let mut v: Vec<Peak> = centroids
                            .into_iter()
                            .map(|p| Peak {
                                mz: p.mz,
                                intensity: p.intensity,
                                ion_mobility: p.im,
                            })
                            .collect();
                        v.sort_by(|a, b| {
                            a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal)
                        });
                        v
                    };

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
                                "Bruker frame {}: {} raw -> {} centroids -> {} post-MAD",
                                n,
                                raw_peak_count,
                                n_centroids,
                                n_after_noise,
                            );
                        } else {
                            log::info!(
                                "Bruker frame {}: {} raw peaks -> {} centroids",
                                n,
                                raw_peak_count,
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

    #[derive(Clone, Copy, Debug)]
    struct ImsPeak {
        mz: f32,
        intensity: f32,
        im: f32,
    }

    /// Flat peak cloud in integer TOF/scan-index space, before m/z and 1/K0 conversion.
    struct RawCloud {
        scan_idx: Vec<u32>,
        tof_idx: Vec<u32>,
        intensity: Vec<f32>,
    }

    impl RawCloud {
        fn from_frame(frame: &timsrust::Frame) -> Self {
            let n = frame.tof_indices.len();
            let mut scan_idx = Vec::with_capacity(n);
            let mut tof_idx = Vec::with_capacity(n);
            let mut intensity = Vec::with_capacity(n);
            for (s, w) in frame.scan_offsets.windows(2).enumerate() {
                let (lo, hi) = (w[0], w[1]);
                for i in lo..hi {
                    scan_idx.push(s as u32);
                    tof_idx.push(frame.tof_indices[i]);
                    intensity.push(frame.intensities[i] as f32);
                }
            }
            Self { scan_idx, tof_idx, intensity }
        }

        fn len(&self) -> usize {
            self.intensity.len()
        }

        fn into_ims_peaks(
            self,
            mz_converter: &Tof2MzConverter,
            ims_converter: &Scan2ImConverter,
        ) -> Vec<ImsPeak> {
            let n = self.len();
            let mut peaks = Vec::with_capacity(n);
            for i in 0..n {
                peaks.push(ImsPeak {
                    mz: mz_converter.convert(self.tof_idx[i] as f64) as f32,
                    intensity: self.intensity[i],
                    im: ims_converter.convert(self.scan_idx[i] as f64) as f32,
                });
            }
            peaks
        }
    }

    /// Centered rolling average along the scan-index axis.
    /// Full window = 2 * half_window + 1. Replicas with shifted scan < 0 are
    /// dropped; positive shifts extend beyond the original max scan
    /// (the IM converter handles values slightly outside the observed range).
    fn smooth_scan_axis(cloud: RawCloud, half_window: u32) -> RawCloud {
        if half_window == 0 || cloud.len() == 0 {
            return cloud;
        }
        let window = 2 * half_window + 1;
        let inv_window = 1.0_f32 / window as f32;
        let max_tof = cloud.tof_idx.iter().copied().max().unwrap_or(0);
        let stride = max_tof as u64 + 1;

        let n = cloud.len();
        let mut map: KeyMap = HashMap::with_capacity_and_hasher(n, BuildU64Hasher);
        let hw = half_window as i64;
        for i in 0..n {
            let s = cloud.scan_idx[i] as i64;
            let t = cloud.tof_idx[i] as u64;
            let w = cloud.intensity[i] * inv_window;
            for shift in -hw..=hw {
                let ss = s + shift;
                if ss < 0 {
                    continue;
                }
                let key = (ss as u64) * stride + t;
                *map.entry(key).or_insert(0.0) += w;
            }
        }
        cloud_from_map(map, stride, KeyAxis::Scan)
    }

    /// Centered rolling average along the TOF-index axis.
    /// Full window = 2 * half_window + 1. Negative-shift replicas are dropped;
    /// positive shifts extend beyond the original max TOF (no upper bound).
    fn smooth_tof_axis(cloud: RawCloud, half_window: u32) -> RawCloud {
        if half_window == 0 || cloud.len() == 0 {
            return cloud;
        }
        let window = 2 * half_window + 1;
        let inv_window = 1.0_f32 / window as f32;
        let max_scan = cloud.scan_idx.iter().copied().max().unwrap_or(0);
        let stride = max_scan as u64 + 1;

        let n = cloud.len();
        let mut map: KeyMap = HashMap::with_capacity_and_hasher(n, BuildU64Hasher);
        let hw = half_window as i64;
        for i in 0..n {
            let t = cloud.tof_idx[i] as i64;
            let s = cloud.scan_idx[i] as u64;
            let w = cloud.intensity[i] * inv_window;
            for shift in -hw..=hw {
                let tt = t + shift;
                if tt < 0 {
                    continue;
                }
                let key = (tt as u64) * stride + s;
                *map.entry(key).or_insert(0.0) += w;
            }
        }
        cloud_from_map(map, stride, KeyAxis::Tof)
    }

    #[derive(Clone, Copy)]
    enum KeyAxis {
        /// key = scan * stride + tof  (stride = max_tof + 1)
        Scan,
        /// key = tof * stride + scan  (stride = max_scan + 1)
        Tof,
    }

    /// Drain a deduped key→intensity map into a flat `RawCloud`. Order is not
    /// preserved (HashMap iteration is unspecified) — downstream stages re-sort
    /// by m/z anyway, so this is fine.
    fn cloud_from_map(map: KeyMap, stride: u64, axis: KeyAxis) -> RawCloud {
        let cap = map.len();
        let mut scan_idx: Vec<u32> = Vec::with_capacity(cap);
        let mut tof_idx: Vec<u32> = Vec::with_capacity(cap);
        let mut intensity: Vec<f32> = Vec::with_capacity(cap);

        for (key, sum) in map {
            let (s, t) = match axis {
                KeyAxis::Scan => ((key / stride) as u32, (key % stride) as u32),
                KeyAxis::Tof => ((key % stride) as u32, (key / stride) as u32),
            };
            scan_idx.push(s);
            tof_idx.push(t);
            intensity.push(sum);
        }

        RawCloud { scan_idx, tof_idx, intensity }
    }

    /// Build a max-heap over `(intensity_bits, idx)` so the highest-intensity
    /// peak pops first. Uses `f32::to_bits()`: for non-negative finite floats,
    /// IEEE 754 layout guarantees `a.to_bits() < b.to_bits()` iff `a < b`, so a
    /// max-heap on bits is a max-heap on intensity. All post-smoothing
    /// intensities are >= 0, so this is safe here. Heapify is O(N) via
    /// `BinaryHeap::from(Vec)` — much cheaper than the prior full O(N log N)
    /// sort, especially when we only need the top ~K=10k anchors.
    fn build_intensity_heap(peaks: &[ImsPeak]) -> BinaryHeap<(u32, u32)> {
        let entries: Vec<(u32, u32)> = peaks
            .iter()
            .enumerate()
            .map(|(i, p)| (p.intensity.to_bits(), i as u32))
            .collect();
        BinaryHeap::from(entries)
    }

    /// Greedy intensity-ordered centroiding within `mz_ppm` / `im_pct` tolerances,
    /// followed by a satellite-suppression pass in a wider `sat_window_da` window
    /// using a linear ramp threshold (anchor_raw_intensity at d=0, anchor *
    /// `sat_end_fraction` at d=window).
    ///
    /// `peaks` must be sorted ascending by m/z. The slice is mutated: consumed and
    /// suppressed peaks have their intensity set to -1.0 (sentinel). Anchors are
    /// pulled from `heap` (max on original intensity); the heap's intensity is a
    /// snapshot — stale-popped indices (already consumed/suppressed) are detected
    /// via `peaks[idx].intensity <= 0.0` and skipped.
    fn centroid_with_satellite(
        peaks: &mut [ImsPeak],
        mut heap: BinaryHeap<(u32, u32)>,
        mz_tol_ppm: f32,
        im_tol_pct: f32,
        min_subpeaks: usize,
        sat_window_da: f32,
        sat_end_fraction: f32,
        max_peaks: usize,
    ) -> Vec<ImsPeak> {
        debug_assert!(
            peaks.windows(2).all(|x| x[0].mz <= x[1].mz),
            "peaks must be mz-sorted"
        );

        let utol = mz_tol_ppm / 1e6;
        let im_tol = im_tol_pct / 100.0;
        let satellite_on = sat_window_da > 0.0;

        let mut centroids: Vec<ImsPeak> = Vec::new();
        let mut global_included = 0usize;

        while let Some((_, idx_u32)) = heap.pop() {
            let idx = idx_u32 as usize;
            if peaks[idx].intensity <= 0.0 {
                continue;
            }
            if centroids.len() >= max_peaks {
                break;
            }

            let mz_anchor = peaks[idx].mz;
            let im_anchor = peaks[idx].im;
            let anchor_raw = peaks[idx].intensity;

            let da_tol = mz_anchor * utol;
            let abs_im_tol = im_anchor * im_tol;

            let merge_start = peaks.partition_point(|p| p.mz < mz_anchor - da_tol);
            let merge_end = peaks.partition_point(|p| p.mz <= mz_anchor + da_tol);

            let mut sum_intensity = 0.0_f32;
            let mut sum_mz_weighted = 0.0_f32;
            let mut sum_im_weighted = 0.0_f32;
            let mut count = 0usize;

            for i in merge_start..merge_end {
                let p = &mut peaks[i];
                if p.intensity > 0.0
                    && p.im >= im_anchor - abs_im_tol
                    && p.im <= im_anchor + abs_im_tol
                {
                    sum_intensity += p.intensity;
                    sum_mz_weighted += p.mz * p.intensity;
                    sum_im_weighted += p.im * p.intensity;
                    p.intensity = -1.0;
                    count += 1;
                }
            }

            if count < min_subpeaks {
                global_included += count;
                if global_included >= peaks.len() {
                    break;
                }
                continue;
            }

            let centroid_mz = sum_mz_weighted / sum_intensity;
            let centroid_im = sum_im_weighted / sum_intensity;

            centroids.push(ImsPeak {
                mz: centroid_mz,
                intensity: sum_intensity,
                im: centroid_im,
            });
            global_included += count;

            if satellite_on {
                let sat_start = peaks.partition_point(|p| p.mz < mz_anchor - sat_window_da);
                let sat_end = peaks.partition_point(|p| p.mz <= mz_anchor + sat_window_da);
                let ramp_span = 1.0 - sat_end_fraction;
                for i in sat_start..sat_end {
                    let p = &mut peaks[i];
                    if p.intensity <= 0.0 {
                        continue;
                    }
                    if p.im < im_anchor - abs_im_tol || p.im > im_anchor + abs_im_tol {
                        continue;
                    }
                    let d = (p.mz - mz_anchor).abs();
                    let threshold = anchor_raw * (1.0 - (d / sat_window_da) * ramp_span);
                    if p.intensity < threshold {
                        p.intensity = -1.0;
                    }
                }
            }

            if global_included >= peaks.len() {
                break;
            }
        }

        centroids
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn mk(mz: f32, intensity: f32, im: f32) -> ImsPeak {
            ImsPeak { mz, intensity, im }
        }

        #[test]
        fn smooth_axis_isolated_peak_spreads_to_window() {
            // Single peak at scan=5, tof=100, intensity=10, half_window=1 → 3 cells
            let cloud = RawCloud {
                scan_idx: vec![5],
                tof_idx: vec![100],
                intensity: vec![10.0],
            };
            let out = smooth_scan_axis(cloud, 1);
            // Expect 3 cells: scan ∈ {4, 5, 6}, all tof=100, intensity = 10/3 each
            assert_eq!(out.len(), 3);
            for i in 0..3 {
                assert_eq!(out.tof_idx[i], 100);
                assert!((out.intensity[i] - 10.0 / 3.0).abs() < 1e-5);
            }
            let mut scans: Vec<u32> = out.scan_idx.clone();
            scans.sort();
            assert_eq!(scans, vec![4, 5, 6]);
        }

        #[test]
        fn smooth_axis_clamps_at_zero() {
            // Peak at scan=0, half_window=2 → shifts -2, -1 are dropped, only 0,1,2 remain
            let cloud = RawCloud {
                scan_idx: vec![0],
                tof_idx: vec![10],
                intensity: vec![5.0],
            };
            let out = smooth_scan_axis(cloud, 2);
            assert_eq!(out.len(), 3);
            let mut scans: Vec<u32> = out.scan_idx.clone();
            scans.sort();
            assert_eq!(scans, vec![0, 1, 2]);
            // window divisor is still 5 (2*2+1), so each replica has 5.0/5 = 1.0
            for &v in &out.intensity {
                assert!((v - 1.0).abs() < 1e-5);
            }
        }

        #[test]
        fn smooth_axis_adjacent_peaks_reinforce() {
            // Two peaks at scan=5, scan=6, tof=100 — overlap at scans 5, 6 after window=1
            let cloud = RawCloud {
                scan_idx: vec![5, 6],
                tof_idx: vec![100, 100],
                intensity: vec![9.0, 9.0],
            };
            let out = smooth_scan_axis(cloud, 1);
            // scan 4: 9/3, scan 5: 9/3 + 9/3 = 6, scan 6: 6, scan 7: 3
            let mut paired: Vec<(u32, f32)> =
                out.scan_idx.iter().zip(out.intensity.iter()).map(|(&s, &i)| (s, i)).collect();
            paired.sort_by_key(|&(s, _)| s);
            assert_eq!(paired.len(), 4);
            assert_eq!(paired[0].0, 4);
            assert!((paired[0].1 - 3.0).abs() < 1e-4);
            assert_eq!(paired[1].0, 5);
            assert!((paired[1].1 - 6.0).abs() < 1e-4);
            assert_eq!(paired[2].0, 6);
            assert!((paired[2].1 - 6.0).abs() < 1e-4);
            assert_eq!(paired[3].0, 7);
            assert!((paired[3].1 - 3.0).abs() < 1e-4);
        }

        #[test]
        fn centroid_weighted_mean() {
            // Two peaks within ppm tolerance, intensities 1 and 3 → weighted mz toward the heavy peak
            let mut peaks = vec![mk(500.000, 1.0, 1.0), mk(500.004, 3.0, 1.0)];
            // 8 ppm at 500 = 0.004 Da → both within window
            let heap = build_intensity_heap(&peaks);
            let cs = centroid_with_satellite(&mut peaks, heap, 8.0, 5.0, 1, 0.0, 0.0, 100);
            assert_eq!(cs.len(), 1);
            // weighted mz: (500.000*1 + 500.004*3) / 4 = 500.003
            assert!((cs[0].mz - 500.003).abs() < 1e-4);
            assert!((cs[0].intensity - 4.0).abs() < 1e-5);
        }

        #[test]
        fn satellite_filter_suppresses_shadow() {
            // Anchor at 500.0 intensity 100; satellite at 500.05 (within 0.15 Da) intensity 10.
            // Linear ramp from 100 (d=0) to 0 (d=0.15). At d=0.05, threshold = 100 * (1 - 0.05/0.15)
            // = 66.67. Satellite at 10 < 66.67 → suppressed.
            let mut peaks = vec![mk(500.0, 100.0, 1.0), mk(500.05, 10.0, 1.0)];
            let heap = build_intensity_heap(&peaks);
            let cs = centroid_with_satellite(&mut peaks, heap, 1.0, 5.0, 1, 0.15, 0.0, 100);
            // Only one centroid: the anchor. The satellite is killed and never seeds.
            assert_eq!(cs.len(), 1);
            assert!((cs[0].mz - 500.0).abs() < 1e-4);
        }

        #[test]
        fn satellite_filter_keeps_real_peak_above_ramp() {
            // Same setup but satellite intensity 80 > 66.67 threshold → survives.
            let mut peaks = vec![mk(500.0, 100.0, 1.0), mk(500.05, 80.0, 1.0)];
            let heap = build_intensity_heap(&peaks);
            let cs = centroid_with_satellite(&mut peaks, heap, 1.0, 5.0, 1, 0.15, 0.0, 100);
            assert_eq!(cs.len(), 2);
        }
    }
}

#[cfg(feature = "tdf")]
pub use inner::read_bruker;
