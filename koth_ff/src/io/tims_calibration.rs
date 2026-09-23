//! Bruker timsTOF scan-to-1/K0 conversion using the acquisition calibration
//! stored in each run's `analysis.tdf`.
//!
//! timsrust 0.4 converts a TIMS scan index to inverse reduced mobility with a
//! straight line between `OneOverK0AcqRangeUpper` (scan 0) and
//! `OneOverK0AcqRangeLower` (the largest `NumScans`). Bruker's own software (the
//! timsdata SDK, and every tool built on it or on its mzML export) instead uses
//! the per-run `TimsCalibration` table, which differs from the line by up to
//! ~0.03 1/K0 on a 0.64–1.45 acquisition range. This module is a native port of
//! that calibration, so koth reports the same 1/K0 as the vendor library without
//! linking it.
//!
//! # The ModelType 2 calibration
//!
//! A `TimsCalibration` row carries `ModelType` and ten coefficients `C0..C9`.
//! For `ModelType = 2` (the only type observed in any timsTOF file available to
//! us) the SDK's `tims_scannum_to_oneoverk0` is, for a scan number `s`:
//!
//! ```text
//! V(s)      = C2 + (C2 - C3) / C1 * (C4 + C0 - s)     ramp voltage at scan s
//! g(V)      = V / (C6 * V + C7)                        1/K0 = 1 / (C6 + C7 / V)
//! 1/K0(s)   = g(V)                                     for C8 <= V <= C9
//!           = g(C8) + g'(C8) * (V - C8)                for V < C8
//!           = g(C9) + g'(C9) * (V - C9)                for V > C9
//! g'(V)     = C7 / (C6 * V + C7)^2
//! ```
//!
//! Reading of the coefficients: `C1` is the last scan number, `C2`/`C3` the ramp
//! start/end voltages, `C4` a scan offset (it equals 3600 / scan period in µs on
//! every file examined), `C6`/`C7` the mobility-voltage relation, and `C8`/`C9`
//! the voltage range outside which the SDK extrapolates linearly (C¹-continuous)
//! instead of evaluating the hyperbola. `C0` and `C5` were 1 in every file, so the
//! data cannot distinguish `C4 + C0` from `C4 + C5` or `C4 + 1` in the offset;
//! `C0` is used.
//!
//! The model was derived by fitting the SDK's output and then verified against
//! `libtimsdata.so` (SHA-256 `c1fbe908…3151`, as shipped with AlphaPept) on 304
//! local `.d` runs covering 25 distinct calibration rows (including a run with two
//! rows), 933 run/frame checks in all, at fractional scan steps over the
//! acquisition range and far outside it: the maximum absolute difference was
//! 8.9e-16 1/K0 inside the scan range and 1.4e-14 outside it (floating-point
//! rounding). Any other `ModelType` is refused rather than approximated.
//!
//! The same conversion could live in timsrust's `Scan2ImConverter`; it is kept
//! here as a small, self-contained module so it can be upstreamed unchanged.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

/// Which 1/K0 scale a Bruker reader reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MobilityScale {
    /// Bruker acquisition calibration (`TimsCalibration`), identical to the
    /// timsdata SDK opened in its default (non-recalibrated) state.
    #[default]
    Calibrated,
    /// timsrust's straight line between the acquisition-range bounds. This was
    /// koth's only behaviour up to 0.9.0.
    Linear,
}

impl MobilityScale {
    /// Stable label written into `report.json` so a downstream join can refuse
    /// to mix scales.
    pub fn label(self) -> &'static str {
        match self {
            MobilityScale::Calibrated => crate::input::SCALE_CALIBRATED,
            MobilityScale::Linear => crate::input::SCALE_LINEAR,
        }
    }
}

/// The only calibration model type implemented (and the only one observed).
pub const SUPPORTED_MODEL_TYPE: i64 = 2;

/// One `TimsCalibration` row of `ModelType` 2.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimsCalibrationModel {
    c: [f64; 10],
}

impl TimsCalibrationModel {
    /// Build from `ModelType` and `C0..C9`. Errors on an unsupported model type
    /// or on coefficients that would divide by zero.
    pub fn new(model_type: i64, c: [f64; 10]) -> Result<Self, String> {
        if model_type != SUPPORTED_MODEL_TYPE {
            return Err(format!(
                "TimsCalibration ModelType {model_type} is not supported (only {SUPPORTED_MODEL_TYPE}); \
                 set `bruker_mobility_scale = \"linear\"` to use the uncalibrated scale"
            ));
        }
        if !c.iter().all(|v| v.is_finite()) || c[1] == 0.0 || c[7] == 0.0 || c[8] > c[9] {
            return Err(format!("invalid TimsCalibration coefficients {c:?}"));
        }
        Ok(Self { c })
    }

    /// `C0..C9` as stored.
    pub fn coefficients(&self) -> [f64; 10] {
        self.c
    }

    #[inline]
    fn g(&self, v: f64) -> f64 {
        v / (self.c[6] * v + self.c[7])
    }

    #[inline]
    fn dg(&self, v: f64) -> f64 {
        let d = self.c[6] * v + self.c[7];
        self.c[7] / (d * d)
    }

    /// Ramp voltage at scan number `scan` (fractional scans allowed).
    #[inline]
    pub fn voltage(&self, scan: f64) -> f64 {
        let c = &self.c;
        c[2] + (c[2] - c[3]) / c[1] * (c[4] + c[0] - scan)
    }

    /// Inverse reduced mobility (1/K0, V·s/cm²) at scan number `scan`, exactly
    /// as the timsdata SDK's `tims_scannum_to_oneoverk0` computes it.
    #[inline]
    pub fn scan_to_inv_k0(&self, scan: f64) -> f64 {
        let v = self.voltage(scan);
        let (lo, hi) = (self.c[8], self.c[9]);
        if v < lo {
            self.g(lo) + self.dg(lo) * (v - lo)
        } else if v > hi {
            self.g(hi) + self.dg(hi) * (v - hi)
        } else {
            self.g(v)
        }
    }
}

/// timsrust 0.4's linear scale, reproduced so the `linear` option and the
/// calibrated path share one code path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearScale {
    upper: f64,
    slope: f64,
}

impl LinearScale {
    /// `im_lower`/`im_upper` are `OneOverK0AcqRange{Lower,Upper}`; `max_scans`
    /// is the largest `Frames.NumScans`.
    pub fn new(im_lower: f64, im_upper: f64, max_scans: u32) -> Self {
        Self {
            upper: im_upper,
            slope: (im_lower - im_upper) / max_scans as f64,
        }
    }

    #[inline]
    pub fn scan_to_inv_k0(&self, scan: f64) -> f64 {
        self.upper + self.slope * scan
    }
}

/// Per-run scan → 1/K0 conversion honouring each frame's calibration row.
///
/// Integer scans (the only ones the readers produce) are served from a
/// precomputed table per calibration row, so the per-peak cost is one lookup.
#[derive(Debug, Clone)]
pub struct ScanToMobility {
    scale: MobilityScale,
    /// Calibration row id per frame id (`Frames.Id` → `Frames.TimsCalibration`).
    frame_calibration: HashMap<usize, i64>,
    /// Row id → (model, table over scans `0..=max_scans`).
    models: HashMap<i64, (TimsCalibrationModel, Vec<f64>)>,
    linear: LinearScale,
    linear_table: Vec<f64>,
}

impl ScanToMobility {
    /// Build from already-read metadata. `frame_calibration` may be empty only
    /// when `scale` is [`MobilityScale::Linear`].
    pub fn new(
        scale: MobilityScale,
        linear: LinearScale,
        max_scans: u32,
        rows: Vec<(i64, TimsCalibrationModel)>,
        frame_calibration: HashMap<usize, i64>,
    ) -> Result<Self, String> {
        let table = |f: &dyn Fn(f64) -> f64| (0..=max_scans).map(|s| f(s as f64)).collect();
        let linear_table = table(&|s| linear.scan_to_inv_k0(s));
        let mut models = HashMap::new();
        if scale == MobilityScale::Calibrated {
            for (id, m) in rows {
                models.insert(id, (m, table(&|s| m.scan_to_inv_k0(s))));
            }
            if let Some((frame, id)) = frame_calibration
                .iter()
                .find(|(_, id)| !models.contains_key(*id))
            {
                return Err(format!(
                    "frame {frame} references TimsCalibration row {id}, which does not exist"
                ));
            }
            if frame_calibration.is_empty() {
                return Err("no frame carries a TimsCalibration id".into());
            }
        }
        Ok(Self {
            scale,
            frame_calibration,
            models,
            linear,
            linear_table,
        })
    }

    pub fn scale(&self) -> MobilityScale {
        self.scale
    }

    /// 1/K0 of integer scan `scan` in frame `frame_id` (`Frames.Id`).
    #[inline]
    pub fn convert(&self, frame_id: usize, scan: u32) -> f64 {
        match self.scale {
            MobilityScale::Linear => self
                .linear_table
                .get(scan as usize)
                .copied()
                .unwrap_or_else(|| self.linear.scan_to_inv_k0(scan as f64)),
            MobilityScale::Calibrated => {
                let id = self.frame_calibration[&frame_id];
                let (model, table) = &self.models[&id];
                table
                    .get(scan as usize)
                    .copied()
                    .unwrap_or_else(|| model.scan_to_inv_k0(scan as f64))
            }
        }
    }

    /// Distinct calibration rows used by the run (empty for the linear scale).
    pub fn calibration_rows(&self) -> Vec<(i64, TimsCalibrationModel)> {
        let mut v: Vec<_> = self.models.iter().map(|(id, (m, _))| (*id, *m)).collect();
        v.sort_by_key(|(id, _)| *id);
        v
    }
}

/// Read the conversion for a `.d` folder straight from its `analysis.tdf`.
#[cfg(feature = "tdf")]
pub fn load(path: &std::path::Path, scale: MobilityScale) -> Result<ScanToMobility, String> {
    use rusqlite::{Connection, OpenFlags};

    let tdf = path.join("analysis.tdf");
    let con = Connection::open_with_flags(&tdf, OpenFlags::SQLITE_OPEN_READ_ONLY)
        .map_err(|e| format!("{}: {e}", tdf.display()))?;
    let meta = |key: &str| -> Result<f64, String> {
        let v: String = con
            .query_row(
                "SELECT Value FROM GlobalMetadata WHERE Key = ?1",
                [key],
                |r| r.get(0),
            )
            .map_err(|e| format!("GlobalMetadata {key}: {e}"))?;
        v.trim()
            .parse()
            .map_err(|e| format!("GlobalMetadata {key} = {v:?}: {e}"))
    };
    let lower = meta("OneOverK0AcqRangeLower")?;
    let upper = meta("OneOverK0AcqRangeUpper")?;
    let max_scans: u32 = con
        .query_row("SELECT MAX(NumScans) FROM Frames", [], |r| r.get(0))
        .map_err(|e| format!("Frames.NumScans: {e}"))?;
    let linear = LinearScale::new(lower, upper, max_scans);

    let mut rows = Vec::new();
    let mut frame_calibration = HashMap::new();
    if scale == MobilityScale::Calibrated {
        let err = |e: rusqlite::Error| {
            format!(
                "cannot read TimsCalibration ({e}); set `bruker_mobility_scale = \"linear\"` \
                 to use the uncalibrated scale"
            )
        };
        let mut stmt = con
            .prepare(
                "SELECT Id, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 FROM TimsCalibration",
            )
            .map_err(err)?;
        let parsed = stmt
            .query_map([], |r| {
                let mut c = [0.0; 10];
                for (i, v) in c.iter_mut().enumerate() {
                    *v = r.get(i + 2)?;
                }
                Ok((r.get::<_, i64>(0)?, r.get::<_, i64>(1)?, c))
            })
            .map_err(err)?;
        for row in parsed {
            let (id, model_type, c) = row.map_err(err)?;
            rows.push((id, TimsCalibrationModel::new(model_type, c)?));
        }
        let mut stmt = con
            .prepare("SELECT Id, TimsCalibration FROM Frames")
            .map_err(err)?;
        let frames = stmt
            .query_map([], |r| {
                Ok((r.get::<_, i64>(0)? as usize, r.get::<_, i64>(1)?))
            })
            .map_err(err)?;
        for f in frames {
            let (frame, id) = f.map_err(err)?;
            frame_calibration.insert(frame, id);
        }
    }
    ScanToMobility::new(scale, linear, max_scans, rows, frame_calibration)
        .map_err(|e| format!("{}: {e}", tdf.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Deserialize)]
    struct Fixture {
        cases: Vec<Case>,
    }

    #[derive(Deserialize)]
    struct Case {
        run: String,
        calibration_id: i64,
        model_type: i64,
        coefficients: Vec<f64>,
        num_scans: u32,
        one_over_k0_acq_range: [f64; 2],
        scans: Vec<f64>,
        sdk_one_over_k0: Vec<f64>,
    }

    fn fixture() -> Fixture {
        let text = include_str!("../../tests/data/tims_calibration_sdk.json");
        serde_json::from_str(text).unwrap()
    }

    fn model(case: &Case) -> TimsCalibrationModel {
        TimsCalibrationModel::new(
            case.model_type,
            case.coefficients.clone().try_into().unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn matches_sdk_on_real_runs() {
        let fx = fixture();
        // Two runs of the benchmark cohort plus three other instruments/methods,
        // one of which has two calibration rows.
        assert!(fx.cases.len() >= 5);
        assert!(fx.cases.iter().filter(|c| c.run.contains("15min")).count() >= 2);
        let mut worst: f64 = 0.0;
        for case in &fx.cases {
            let m = model(case);
            for (&s, &want) in case.scans.iter().zip(&case.sdk_one_over_k0) {
                let err = (m.scan_to_inv_k0(s) - want).abs();
                assert!(
                    err <= 1e-12,
                    "{} cal {} scan {s}: got {} want {want}",
                    case.run,
                    case.calibration_id,
                    m.scan_to_inv_k0(s)
                );
                worst = worst.max(err);
            }
        }
        assert!(worst <= 1e-12);
    }

    #[test]
    fn table_lookup_equals_direct_evaluation_and_frames_route_to_their_row() {
        let fx = fixture();
        let two: Vec<&Case> = fx.cases.iter().filter(|c| c.run.contains("IRT")).collect();
        assert_eq!(two.len(), 2, "fixture keeps the two-calibration run");
        let rows: Vec<_> = two.iter().map(|c| (c.calibration_id, model(c))).collect();
        let frames = HashMap::from([(10, two[0].calibration_id), (11, two[1].calibration_id)]);
        let [lo, hi] = two[0].one_over_k0_acq_range;
        let lin = LinearScale::new(lo, hi, two[0].num_scans);
        let conv = ScanToMobility::new(
            MobilityScale::Calibrated,
            lin,
            two[0].num_scans,
            rows.clone(),
            frames,
        )
        .unwrap();
        for scan in [0u32, 1, 500, two[0].num_scans, two[0].num_scans + 7] {
            assert_eq!(
                conv.convert(10, scan),
                rows[0].1.scan_to_inv_k0(scan as f64)
            );
            assert_eq!(
                conv.convert(11, scan),
                rows[1].1.scan_to_inv_k0(scan as f64)
            );
        }
        assert_ne!(conv.convert(10, 500), conv.convert(11, 500));
    }

    #[cfg(feature = "tdf")]
    #[test]
    fn linear_scale_reproduces_timsrust() {
        use timsrust::converters::{ConvertableDomain, Scan2ImConverter};
        let fx = fixture();
        let case = &fx.cases[0];
        let [lo, hi] = case.one_over_k0_acq_range;
        let ours = LinearScale::new(lo, hi, case.num_scans);
        let theirs = Scan2ImConverter::from_boundaries(lo, hi, case.num_scans);
        let conv = ScanToMobility::new(
            MobilityScale::Linear,
            ours,
            case.num_scans,
            vec![],
            HashMap::new(),
        )
        .unwrap();
        for scan in 0..=case.num_scans + 3 {
            assert_eq!(conv.convert(1, scan), theirs.convert(scan as f64));
        }
    }

    #[test]
    fn calibrated_and_linear_differ_on_the_cohort_calibration() {
        // The reason this module exists: on the benchmark calibration the two
        // scales disagree by ~0.032 1/K0 at scan 0.
        let fx = fixture();
        let case = &fx.cases[0];
        let [lo, hi] = case.one_over_k0_acq_range;
        let d = model(case).scan_to_inv_k0(0.0)
            - LinearScale::new(lo, hi, case.num_scans).scan_to_inv_k0(0.0);
        assert!((d - 0.0318).abs() < 1e-3, "{d}");
    }

    #[test]
    fn refuses_unknown_model_type_and_missing_rows() {
        let c = [
            1.0, 935.0, 239.3, 103.0, 33.6, 1.0, -0.027, 171.4, 16.8, 1732.7,
        ];
        assert!(TimsCalibrationModel::new(1, c).is_err());
        let m = TimsCalibrationModel::new(2, c).unwrap();
        let lin = LinearScale::new(0.6, 1.4, 935);
        let bad = HashMap::from([(1usize, 9i64)]);
        assert!(
            ScanToMobility::new(MobilityScale::Calibrated, lin, 935, vec![(1, m)], bad).is_err()
        );
        assert!(ScanToMobility::new(
            MobilityScale::Calibrated,
            lin,
            935,
            vec![(1, m)],
            HashMap::new()
        )
        .is_err());
    }

    #[test]
    fn scale_labels_are_stable() {
        assert_eq!(MobilityScale::default(), MobilityScale::Calibrated);
        assert_eq!(
            MobilityScale::Calibrated.label(),
            "bruker-acquisition-calibrated-1/K0"
        );
        assert_eq!(MobilityScale::Linear.label(), "timsrust-linear-1/K0");
    }
}
