//! Per-scan anomaly normalisation applied before hill detection.
//!
//! Computes a per-scan reference quantity (TIC or median peak intensity),
//! takes a centred rolling median over a wide window, and scales each
//! scan's peak intensities by `clamp(local_median / scan_reference,
//! min_scale, max_scale)`. This pulls back peak intensities in scans
//! whose reference dipped sharply relative to the local elution context
//! (e.g., brief ESI dropouts the instrument's AGC didn't fully
//! compensate) without touching scans whose reference matches local
//! elution.
//!
//! Reference mode `tic` (sum of peak intensities) is dominated by a few
//! intense peaks — real elution apices spike the TIC and get
//! systematically scaled down by the normaliser, flattening apex shape.
//! Reference mode `median` (median peak intensity per scan) is robust
//! to a handful of intense peaks: an apex elution adds tall peaks but
//! barely moves the median, while a dropout drops the median for every
//! peak in the scan. `median` is the preferred mode for Orbitrap DDA.

use crate::models::Spectrum;

/// Which per-scan quantity to use as the rolling reference.
#[derive(Debug, Clone, Copy)]
pub enum RefMode {
    /// Sum of all peak intensities (TIC).
    Tic,
    /// Median peak intensity per scan.
    Median,
}

impl RefMode {
    pub fn parse(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "tic" => RefMode::Tic,
            _ => RefMode::Median,
        }
    }
}

/// Apply rolling-median TIC anomaly normalisation to a slice of MS1 spectra
/// in place. Non-MS1 spectra are not modified.
///
/// `window`: full centred window size in scans (e.g., 100 = ±50 scans).
///           Must be > 0; `0` means "off" and the caller should not call.
/// `min_scale, max_scale`: clamping bounds on the per-scan multiplier.
///                         Typical (0.5, 3.0).
pub fn normalize_tic_in_place(
    spectra: &mut [Spectrum],
    window: usize,
    min_scale: f64,
    max_scale: f64,
    mode: RefMode,
) {
    if window == 0 || spectra.len() < 2 {
        return;
    }

    // Indices of MS1 spectra and per-scan reference quantities.
    let ms1_idx: Vec<usize> = spectra
        .iter()
        .enumerate()
        .filter_map(|(i, s)| if s.ms_level == 1 { Some(i) } else { None })
        .collect();
    if ms1_idx.len() < 2 {
        return;
    }

    let mut tics: Vec<f64> = Vec::with_capacity(ms1_idx.len());
    let mut buf: Vec<f64> = Vec::new();
    for &i in &ms1_idx {
        let val = match mode {
            RefMode::Tic => spectra[i].peaks.iter().map(|p| p.intensity as f64).sum(),
            RefMode::Median => {
                buf.clear();
                buf.extend(spectra[i].peaks.iter().map(|p| p.intensity as f64));
                if buf.is_empty() {
                    0.0
                } else {
                    buf.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                    let m = buf.len();
                    if m % 2 == 1 {
                        buf[m / 2]
                    } else {
                        0.5 * (buf[m / 2 - 1] + buf[m / 2])
                    }
                }
            }
        };
        tics.push(val);
    }

    let half = window / 2;
    let n = tics.len();

    // Rolling median over centred window. Reuses a sortable scratch Vec.
    let mut scratch: Vec<f64> = Vec::with_capacity(window + 1);
    let mut log_scales: Vec<f64> = Vec::with_capacity(n);
    let mut n_clamped_low = 0usize;
    let mut n_clamped_high = 0usize;
    let mut sum_abs_log_scale = 0.0f64;

    for i in 0..n {
        let lo = i.saturating_sub(half);
        let hi = (i + half + 1).min(n);
        scratch.clear();
        scratch.extend_from_slice(&tics[lo..hi]);
        scratch.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let m = scratch.len();
        let median = if m == 0 {
            tics[i].max(1.0)
        } else if m % 2 == 1 {
            scratch[m / 2]
        } else {
            0.5 * (scratch[m / 2 - 1] + scratch[m / 2])
        };

        let raw_scale = if tics[i] > 0.0 {
            median / tics[i]
        } else {
            max_scale
        };
        let clamped = raw_scale.clamp(min_scale, max_scale);
        if raw_scale < min_scale {
            n_clamped_low += 1;
        } else if raw_scale > max_scale {
            n_clamped_high += 1;
        }
        log_scales.push(clamped.ln());
        sum_abs_log_scale += clamped.ln().abs();

        if (clamped - 1.0).abs() > 1e-9 {
            let sp_idx = ms1_idx[i];
            let s = clamped as f32;
            for p in &mut spectra[sp_idx].peaks {
                p.intensity *= s;
            }
        }
    }

    let mean_abs_log = sum_abs_log_scale / n as f64;
    let mean_distort = (mean_abs_log.exp() - 1.0) * 100.0;
    let mode_name = match mode {
        RefMode::Tic => "tic",
        RefMode::Median => "median",
    };
    log::info!(
        "TIC normalisation: mode={} window={} scans, ({}/{}) clamped low/high, \
         mean |Δlog scale| = {:.4} (~{:.1}% mean intensity distortion)",
        mode_name,
        window,
        n_clamped_low,
        n_clamped_high,
        mean_abs_log,
        mean_distort,
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::{Peak, Spectrum};

    fn spec(rt: f64, intensities: &[f32]) -> Spectrum {
        Spectrum {
            scan_index: 0,
            ms_level: 1,
            retention_time: rt,
            peaks: intensities
                .iter()
                .enumerate()
                .map(|(i, &v)| Peak {
                    mz: 100.0 + i as f32,
                    intensity: v,
                    ion_mobility: 0.0,
                })
                .collect(),
            isolation_window: None,
        }
    }

    #[test]
    fn off_when_window_is_zero() {
        let mut s = vec![spec(0.0, &[1.0, 2.0]), spec(1.0, &[1.0, 2.0])];
        let orig = s.iter().map(|sp| sp.peaks[0].intensity).collect::<Vec<_>>();
        normalize_tic_in_place(&mut s, 0, 0.5, 3.0, RefMode::Tic);
        let after = s.iter().map(|sp| sp.peaks[0].intensity).collect::<Vec<_>>();
        assert_eq!(orig, after);
    }

    #[test]
    fn boosts_anomalously_low_scan_tic_mode() {
        let mut s: Vec<Spectrum> = (0..21)
            .map(|i| {
                let t = if i == 10 { 10.0 } else { 100.0 };
                spec(i as f64, &[t])
            })
            .collect();
        normalize_tic_in_place(&mut s, 11, 0.5, 3.0, RefMode::Tic);
        assert!((s[10].peaks[0].intensity - 30.0).abs() < 1e-3);
        assert!((s[5].peaks[0].intensity - 100.0).abs() < 1e-3);
    }

    #[test]
    fn median_mode_ignores_apex_but_corrects_dropout() {
        // 21 scans. Most have peaks [10, 20, 30] (median 20). Scan 10 has
        // a giant apex peak: [10, 20, 30, 10000] (median ~25, similar to
        // baseline) — TIC mode would scale apex DOWN here; median mode
        // should leave it alone. Scan 15 is a uniform dropout: [1, 2, 3]
        // (median 2, 10× lower) — both modes should boost it.
        let mut s: Vec<Spectrum> = (0..21)
            .map(|i| {
                if i == 10 {
                    spec(i as f64, &[10.0, 20.0, 30.0, 10000.0])
                } else if i == 15 {
                    spec(i as f64, &[1.0, 2.0, 3.0])
                } else {
                    spec(i as f64, &[10.0, 20.0, 30.0])
                }
            })
            .collect();
        normalize_tic_in_place(&mut s, 11, 0.5, 3.0, RefMode::Median);
        // Apex scan median is 25 vs window median 20 → scale ≈ 0.8.
        // Apex peak 10000 should be barely changed (within ~25%).
        assert!(
            s[10].peaks[3].intensity > 7000.0,
            "apex peak should not be dragged down too much in median mode (got {})",
            s[10].peaks[3].intensity
        );
        // Dropout scan median is 2 vs window median 20 → ratio 10×, clamped to 3.0.
        // Peak 3.0 should become 9.0.
        assert!(
            (s[15].peaks[2].intensity - 9.0).abs() < 0.5,
            "dropout peak should be boosted ~3× (got {})",
            s[15].peaks[2].intensity
        );
    }
}
