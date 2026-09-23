//! Bruker timsTOF scan-to-1/K0 conversion using the acquisition calibration
//! stored in each run's `analysis.tdf`.
//!
//! timsrust 0.4 converts a TIMS scan index to inverse reduced mobility with a
//! straight line between `OneOverK0AcqRangeUpper` (scan 0) and
//! `OneOverK0AcqRangeLower` (the largest `NumScans`). Bruker's own software (the
//! timsdata SDK, and every tool built on it or on its mzML export) instead uses
//! the per-run `TimsCalibration` table, which differs from the line by up to
//! ~0.03 1/K0 on a 0.64–1.45 acquisition range.
//!
//! The calibration model itself (`TimsCalibration` `ModelType` 2, a native port
//! of the SDK's `tims_scannum_to_oneoverk0`) lives in [`dnoise::mobility`], so
//! koth's reported 1/K0 and dnoise's MS1 gates use one implementation. It was
//! verified against `libtimsdata.so` (SHA-256 `c1fbe908…3151`, as shipped with
//! AlphaPept) on 304 local `.d` runs covering 25 distinct calibration rows: the
//! maximum absolute difference was 8.9e-16 1/K0 inside the scan range and
//! 1.4e-14 outside it. The SDK reference values in
//! `tests/data/tims_calibration_sdk.json` are checked below. Any other
//! `ModelType` is refused rather than approximated.
//!
//! This module adds what koth needs on top: the config-facing [`MobilityScale`]
//! (available without the `tdf` feature), per-frame routing for runs with more
//! than one calibration row, per-scan lookup tables, and the coefficients that
//! `report.json` records.

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

/// The one place koth's config enum maps onto dnoise's.
#[cfg(feature = "tdf")]
impl From<MobilityScale> for dnoise::MobilityScale {
    fn from(s: MobilityScale) -> Self {
        match s {
            MobilityScale::Calibrated => dnoise::MobilityScale::Calibrated,
            MobilityScale::Linear => dnoise::MobilityScale::Linear,
        }
    }
}

#[cfg(feature = "tdf")]
pub use dnoise::mobility::{TimsCalibrationModel, SUPPORTED_MODEL_TYPE};

#[cfg(feature = "tdf")]
pub use tdf::*;

#[cfg(feature = "tdf")]
mod tdf {
    use std::collections::HashMap;

    use dnoise::mobility::ScanToMobility as RowConverter;
    use timsrust::converters::Scan2ImConverter;

    use super::{MobilityScale, TimsCalibrationModel};

    /// One scale (a calibration row, or the linear line) with its per-scan table.
    #[derive(Debug, Clone)]
    struct Table {
        conv: RowConverter,
        /// 1/K0 at integer scans `0..=max_scans`.
        values: Vec<f64>,
    }

    impl Table {
        fn new(conv: RowConverter, max_scans: u32) -> Self {
            let values = (0..=max_scans).map(|s| conv.convert(s as f64)).collect();
            Self { conv, values }
        }

        #[inline]
        fn get(&self, scan: u32) -> f64 {
            self.values
                .get(scan as usize)
                .copied()
                .unwrap_or_else(|| self.conv.convert(scan as f64))
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
        /// Row id → table (calibrated scale only).
        rows: HashMap<i64, Table>,
        /// Row id → `C0..C9` as stored, for `report.json`.
        coefficients: Vec<(i64, [f64; 10])>,
        linear: Table,
    }

    impl ScanToMobility {
        /// Build from already-read metadata. `rows` are `(Id, ModelType, C0..C9)`
        /// of `TimsCalibration`; they and `frame_calibration` may be empty only
        /// when `scale` is [`MobilityScale::Linear`]. `linear` is timsrust's
        /// converter for the run.
        pub fn new(
            scale: MobilityScale,
            linear: Scan2ImConverter,
            max_scans: u32,
            rows: Vec<(i64, i64, [f64; 10])>,
            frame_calibration: HashMap<usize, i64>,
        ) -> Result<Self, String> {
            let mut tables = HashMap::new();
            let mut coefficients = Vec::new();
            if scale == MobilityScale::Calibrated {
                for (id, model_type, c) in rows {
                    let model = TimsCalibrationModel::new(model_type, c)?;
                    tables.insert(id, Table::new(RowConverter::Calibrated(model), max_scans));
                    coefficients.push((id, c));
                }
                coefficients.sort_by_key(|(id, _)| *id);
                if let Some((frame, id)) = frame_calibration
                    .iter()
                    .find(|(_, id)| !tables.contains_key(*id))
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
                rows: tables,
                coefficients,
                linear: Table::new(RowConverter::Linear(linear), max_scans),
            })
        }

        pub fn scale(&self) -> MobilityScale {
            self.scale
        }

        /// 1/K0 of integer scan `scan` in frame `frame_id` (`Frames.Id`).
        #[inline]
        pub fn convert(&self, frame_id: usize, scan: u32) -> f64 {
            match self.scale {
                MobilityScale::Linear => self.linear.get(scan),
                MobilityScale::Calibrated => {
                    self.rows[&self.frame_calibration[&frame_id]].get(scan)
                }
            }
        }

        /// `(Id, C0..C9)` of the calibration rows used by the run, by id (empty
        /// for the linear scale).
        pub fn calibration_rows(&self) -> &[(i64, [f64; 10])] {
            &self.coefficients
        }
    }

    /// Read the conversion for a `.d` folder straight from its `analysis.tdf`.
    pub fn load(path: &std::path::Path, scale: MobilityScale) -> Result<ScanToMobility, String> {
        use rusqlite::{Connection, OpenFlags};

        let tdf = path.join("analysis.tdf");
        let at = |e: String| format!("{}: {e}", tdf.display());
        // The readers' own linear converter, so `linear` matches them exactly.
        let linear = timsrust::readers::MetadataReader::new(&tdf)
            .map_err(|e| at(e.to_string()))?
            .im_converter;
        let con = Connection::open_with_flags(&tdf, OpenFlags::SQLITE_OPEN_READ_ONLY)
            .map_err(|e| at(e.to_string()))?;
        let max_scans: u32 = con
            .query_row("SELECT MAX(NumScans) FROM Frames", [], |r| r.get(0))
            .map_err(|e| at(format!("Frames.NumScans: {e}")))?;

        let mut rows = Vec::new();
        let mut frame_calibration = HashMap::new();
        if scale == MobilityScale::Calibrated {
            let err = |e: rusqlite::Error| {
                at(format!(
                    "cannot read TimsCalibration ({e}); set `bruker_mobility_scale = \"linear\"` \
                     to use the uncalibrated scale"
                ))
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
                rows.push(row.map_err(err)?);
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
        ScanToMobility::new(scale, linear, max_scans, rows, frame_calibration).map_err(at)
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use serde::Deserialize;
        use timsrust::converters::ConvertableDomain;

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

        fn coeffs(case: &Case) -> [f64; 10] {
            case.coefficients.clone().try_into().unwrap()
        }

        fn model(case: &Case) -> TimsCalibrationModel {
            TimsCalibrationModel::new(case.model_type, coeffs(case)).unwrap()
        }

        fn linear(case: &Case) -> Scan2ImConverter {
            let [lo, hi] = case.one_over_k0_acq_range;
            Scan2ImConverter::from_boundaries(lo, hi, case.num_scans)
        }

        /// A one-row calibrated converter for `case`, every frame on that row.
        fn converter(case: &Case) -> ScanToMobility {
            ScanToMobility::new(
                MobilityScale::Calibrated,
                linear(case),
                case.num_scans,
                vec![(case.calibration_id, case.model_type, coeffs(case))],
                HashMap::from([(1, case.calibration_id)]),
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
            for case in &fx.cases {
                let m = model(case);
                let conv = converter(case);
                for (&s, &want) in case.scans.iter().zip(&case.sdk_one_over_k0) {
                    let got = m.scan_to_inv_k0(s);
                    assert!(
                        (got - want).abs() <= 1e-12,
                        "{} cal {} scan {s}: got {got} want {want}",
                        case.run,
                        case.calibration_id,
                    );
                    // koth's per-frame table serves the same values at integer scans.
                    if s >= 0.0 && s.fract() == 0.0 {
                        let t = conv.convert(1, s as u32);
                        assert!((t - want).abs() <= 1e-12, "{} table scan {s}", case.run);
                    }
                }
            }
        }

        #[test]
        fn table_lookup_equals_direct_evaluation_and_frames_route_to_their_row() {
            let fx = fixture();
            let two: Vec<&Case> = fx.cases.iter().filter(|c| c.run.contains("IRT")).collect();
            assert_eq!(two.len(), 2, "fixture keeps the two-calibration run");
            let rows: Vec<_> = two
                .iter()
                .map(|c| (c.calibration_id, c.model_type, coeffs(c)))
                .collect();
            let frames = HashMap::from([(10, two[0].calibration_id), (11, two[1].calibration_id)]);
            let conv = ScanToMobility::new(
                MobilityScale::Calibrated,
                linear(two[0]),
                two[0].num_scans,
                rows,
                frames,
            )
            .unwrap();
            for scan in [0u32, 1, 500, two[0].num_scans, two[0].num_scans + 7] {
                assert_eq!(
                    conv.convert(10, scan),
                    model(two[0]).scan_to_inv_k0(scan as f64)
                );
                assert_eq!(
                    conv.convert(11, scan),
                    model(two[1]).scan_to_inv_k0(scan as f64)
                );
            }
            assert_ne!(conv.convert(10, 500), conv.convert(11, 500));
            let ids: Vec<i64> = conv.calibration_rows().iter().map(|(id, _)| *id).collect();
            assert!(ids.windows(2).all(|w| w[0] < w[1]) && ids.len() == 2);
        }

        #[test]
        fn linear_scale_reproduces_timsrust() {
            let fx = fixture();
            let case = &fx.cases[0];
            let theirs = linear(case);
            let conv = ScanToMobility::new(
                MobilityScale::Linear,
                theirs,
                case.num_scans,
                vec![],
                HashMap::new(),
            )
            .unwrap();
            for scan in 0..=case.num_scans + 3 {
                assert_eq!(conv.convert(1, scan), theirs.convert(scan as f64));
            }
            assert!(conv.calibration_rows().is_empty());
        }

        #[test]
        fn calibrated_and_linear_differ_on_the_cohort_calibration() {
            // On the benchmark calibration the two scales disagree by ~0.032
            // 1/K0 at scan 0.
            let fx = fixture();
            let case = &fx.cases[0];
            let d = model(case).scan_to_inv_k0(0.0) - linear(case).convert(0.0);
            assert!((d - 0.0318).abs() < 1e-3, "{d}");
        }

        #[test]
        fn refuses_unknown_model_type_and_missing_rows() {
            let c = [
                1.0, 935.0, 239.3, 103.0, 33.6, 1.0, -0.027, 171.4, 16.8, 1732.7,
            ];
            let lin = Scan2ImConverter::from_boundaries(0.6, 1.4, 935);
            let one = |mt| vec![(1i64, mt, c)];
            let frames = || HashMap::from([(1usize, 1i64)]);
            assert!(
                ScanToMobility::new(MobilityScale::Calibrated, lin, 935, one(1), frames()).is_err()
            );
            assert!(
                ScanToMobility::new(MobilityScale::Calibrated, lin, 935, one(2), frames()).is_ok()
            );
            let bad = HashMap::from([(1usize, 9i64)]);
            assert!(ScanToMobility::new(MobilityScale::Calibrated, lin, 935, one(2), bad).is_err());
            assert!(ScanToMobility::new(
                MobilityScale::Calibrated,
                lin,
                935,
                one(2),
                HashMap::new()
            )
            .is_err());
        }

        #[test]
        fn config_scale_maps_onto_dnoise() {
            assert_eq!(
                dnoise::MobilityScale::from(MobilityScale::Calibrated),
                dnoise::MobilityScale::Calibrated
            );
            assert_eq!(
                dnoise::MobilityScale::from(MobilityScale::Linear),
                dnoise::MobilityScale::Linear
            );
            assert_eq!(
                dnoise::MobilityScale::from(MobilityScale::default()),
                dnoise::MobilityScale::default()
            );
        }
    }
}

#[cfg(test)]
mod label_tests {
    use super::*;

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
