/// Readers for hill and feature files produced by `koth_ff`.
///
/// Supports both TSV (default) and Parquet output formats.
/// Auto-detection is based on the file extension.
use std::path::Path;
use std::sync::Arc;

use serde::Deserialize;

use crate::error::KothError;
use crate::models::{Feature, Hill, IsolationWindow, ScoredFeature};

const C13_NEUTRON: f64 = 1.003_354_835;

// ──────────────────────────────────────────────────────────────────────────────
// Public API
// ──────────────────────────────────────────────────────────────────────────────

/// Read hills from a file, detecting format from the extension (.tsv / .parquet).
pub fn read_hills(path: &Path) -> Result<Vec<Hill>, KothError> {
    match path.extension().and_then(|e| e.to_str()) {
        Some("parquet") => read_hills_parquet(path),
        _ => read_hills_tsv(path),
    }
}

/// Read scored features from a file, detecting format from extension.
pub fn read_features(path: &Path) -> Result<Vec<ScoredFeature>, KothError> {
    match path.extension().and_then(|e| e.to_str()) {
        Some("parquet") => read_features_parquet(path),
        _ => read_features_tsv(path),
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// TSV readers
// ──────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
struct HillRow {
    #[serde(default, deserialize_with = "de_opt_u64")]
    hill_id: Option<u64>,
    mz: f64,
    mz_std: f64,
    #[serde(default)]
    mz_se: f64,
    rt: f64,
    rt_start: f64,
    rt_end: f64,
    rt_width: f64,
    im: f64,
    im_std: f64,
    scan_start: usize,
    scan_apex: usize,
    scan_end: usize,
    n_scans: usize,
    skipped_scans: usize,
    intensity_sum: f64,
    intensity_max: f64,
    hill_score: f64,
    intensity_profile: String,
    #[serde(default, deserialize_with = "de_opt_float")]
    iso_target_mz: Option<f64>,
    #[serde(default, deserialize_with = "de_opt_float")]
    iso_lower_mz: Option<f64>,
    #[serde(default, deserialize_with = "de_opt_float")]
    iso_upper_mz: Option<f64>,
}

fn de_opt_float<'de, D: serde::Deserializer<'de>>(de: D) -> Result<Option<f64>, D::Error> {
    let s: Option<String> = Option::deserialize(de)?;
    match s {
        Some(s) if !s.is_empty() => s.parse::<f64>().map(Some).map_err(serde::de::Error::custom),
        _ => Ok(None),
    }
}

fn de_opt_u64<'de, D: serde::Deserializer<'de>>(de: D) -> Result<Option<u64>, D::Error> {
    let s: Option<String> = Option::deserialize(de)?;
    match s {
        Some(s) if !s.is_empty() => s.parse::<u64>().map(Some).map_err(serde::de::Error::custom),
        _ => Ok(None),
    }
}

pub fn read_hills_tsv(path: &Path) -> Result<Vec<Hill>, KothError> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let mut hills = Vec::new();
    for (row_idx, result) in rdr.deserialize::<HillRow>().enumerate() {
        let row = result?;
        let profile: Vec<f32> = serde_json::from_str(&row.intensity_profile)?;
        let isolation_window = match (row.iso_target_mz, row.iso_lower_mz, row.iso_upper_mz) {
            (Some(t), Some(l), Some(u)) => Some(IsolationWindow {
                target: t,
                lower: l,
                upper: u,
            }),
            _ => None,
        };
        hills.push(Hill {
            hill_id: row.hill_id.unwrap_or(row_idx as u64),
            mz: row.mz,
            mz_std: row.mz_std,
            mz_se: row.mz_se,
            rt: row.rt,
            rt_start: row.rt_start,
            rt_end: row.rt_end,
            rt_width: row.rt_width,
            im: row.im,
            im_std: row.im_std,
            scan_start: row.scan_start,
            scan_apex: row.scan_apex,
            scan_end: row.scan_end,
            n_scans: row.n_scans,
            skipped_scans: row.skipped_scans,
            intensity_sum: row.intensity_sum,
            intensity_max: row.intensity_max,
            hill_score: row.hill_score,
            intensity_profile: Arc::from(profile.as_slice()),
            isolation_window,
            faims_cv: None,
        });
    }
    Ok(hills)
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct FeatureRow {
    #[serde(rename = "massCalib")]
    mass_calib: String,
    mz: f64,
    #[serde(rename = "rtApex")]
    rt_apex: f64,
    #[serde(rename = "rtStart")]
    rt_start: f64,
    #[serde(rename = "rtEnd")]
    rt_end: f64,
    #[serde(rename = "FAIMS", default, deserialize_with = "de_opt_float")]
    faims_cv: Option<f64>,
    #[serde(rename = "intensityApex")]
    intensity_apex: f64,
    #[serde(rename = "intensitySum")]
    intensity_sum: f64,
    charge: u8,
    #[serde(rename = "nIsotopes")]
    n_isotopes: usize,
    #[serde(rename = "nScans")]
    n_scans: usize,
    #[serde(deserialize_with = "de_optional_float")]
    im: f64,
    cosine_score: f64,
    ppm_error: f64,
    neutron_offset: i8,
    isotope_score: f64,
    combined_score: f64,
    theoretical_pattern: String,
    // remaining columns are not needed for alignment/LFQ
    isotope_profile: String,
    elution_profile: String,
}

pub fn read_features_tsv(path: &Path) -> Result<Vec<ScoredFeature>, KothError> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let mut features = Vec::new();
    for result in rdr.deserialize::<FeatureRow>() {
        let row = result?;
        if row.charge == 0 {
            continue;
        }
        let theoretical_pattern: Vec<f64> = serde_json::from_str(&row.theoretical_pattern)?;
        features.push(scored_feature_from_row(row, theoretical_pattern));
    }
    Ok(features)
}

fn scored_feature_from_row(row: FeatureRow, theoretical_pattern: Vec<f64>) -> ScoredFeature {
    // The stored `mz` is the calibrated monoisotopic mz (neutron-offset corrected).
    // To reconstruct the raw hill mz: hill.mz = mz + offset * C13/charge
    let raw_mz = row.mz + row.neutron_offset as f64 * C13_NEUTRON / row.charge as f64;

    let hill = Hill {
        hill_id: 0,
        mz: raw_mz,
        mz_std: 0.0,
        mz_se: 0.0,
        rt: row.rt_apex,
        rt_start: row.rt_start,
        rt_end: row.rt_end,
        rt_width: row.rt_end - row.rt_start,
        im: row.im,
        im_std: 0.0,
        scan_start: 0,
        scan_apex: 0,
        scan_end: 0,
        n_scans: row.n_scans,
        skipped_scans: 0,
        intensity_sum: row.intensity_sum,
        intensity_max: row.intensity_apex,
        hill_score: 1.0,
        intensity_profile: Arc::from(&[] as &[f32]),
        isolation_window: None,
        faims_cv: row.faims_cv.map(|v| v as f32),
    };

    ScoredFeature {
        feature: Feature {
            hills: vec![hill],
            charge: row.charge,
            cosine_score: row.cosine_score,
            ppm_error: row.ppm_error,
        },
        neutron_offset: row.neutron_offset,
        isotope_score: row.isotope_score,
        cosine_score: row.cosine_score,
        combined_score: row.combined_score,
        theoretical_pattern,
    }
}

fn de_optional_float<'de, D: serde::Deserializer<'de>>(de: D) -> Result<f64, D::Error> {
    let s = String::deserialize(de)?;
    if s.is_empty() {
        Ok(0.0)
    } else {
        s.parse::<f64>().map_err(serde::de::Error::custom)
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Parquet readers
// ──────────────────────────────────────────────────────────────────────────────

pub fn read_hills_parquet(path: &Path) -> Result<Vec<Hill>, KothError> {
    use arrow::array::{Array, Float64Array, Int64Array, StringArray};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = std::fs::File::open(path)?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| KothError::ParquetError(e.to_string()))?
        .build()
        .map_err(|e| KothError::ParquetError(e.to_string()))?;

    let mut hills = Vec::new();

    for batch in reader {
        let batch = batch.map_err(|e| KothError::ParquetError(e.to_string()))?;

        macro_rules! f64_col {
            ($name:expr) => {
                batch
                    .column_by_name($name)
                    .ok_or_else(|| KothError::ParquetError(format!("missing column: {}", $name)))?
                    .as_any()
                    .downcast_ref::<Float64Array>()
                    .ok_or_else(|| KothError::ParquetError(format!("wrong type: {}", $name)))?
            };
        }
        macro_rules! i64_col {
            ($name:expr) => {
                batch
                    .column_by_name($name)
                    .ok_or_else(|| KothError::ParquetError(format!("missing column: {}", $name)))?
                    .as_any()
                    .downcast_ref::<Int64Array>()
                    .ok_or_else(|| KothError::ParquetError(format!("wrong type: {}", $name)))?
            };
        }
        macro_rules! str_col {
            ($name:expr) => {
                batch
                    .column_by_name($name)
                    .ok_or_else(|| KothError::ParquetError(format!("missing column: {}", $name)))?
                    .as_any()
                    .downcast_ref::<StringArray>()
                    .ok_or_else(|| KothError::ParquetError(format!("wrong type: {}", $name)))?
            };
        }

        let mz_col = f64_col!("mz");
        let mz_std_col = f64_col!("mz_std");
        // mz_se was added later; tolerate older parquets by defaulting to 0.
        let mz_se_col = batch
            .column_by_name("mz_se")
            .and_then(|c| c.as_any().downcast_ref::<Float64Array>());
        let rt_col = f64_col!("rt");
        let rt_start_col = f64_col!("rt_start");
        let rt_end_col = f64_col!("rt_end");
        let rt_width_col = f64_col!("rt_width");
        let im_col = f64_col!("im");
        let im_std_col = f64_col!("im_std");
        let scan_start_col = i64_col!("scan_start");
        let scan_apex_col = i64_col!("scan_apex");
        let scan_end_col = i64_col!("scan_end");
        let n_scans_col = i64_col!("n_scans");
        let skipped_col = i64_col!("skipped_scans");
        let int_sum_col = f64_col!("intensity_sum");
        let int_max_col = f64_col!("intensity_max");
        let score_col = f64_col!("hill_score");
        let profile_col = str_col!("intensity_profile");

        let iso_target_col = batch
            .column_by_name("iso_target_mz")
            .and_then(|c| c.as_any().downcast_ref::<Float64Array>());
        let iso_lower_col = batch
            .column_by_name("iso_lower_mz")
            .and_then(|c| c.as_any().downcast_ref::<Float64Array>());
        let iso_upper_col = batch
            .column_by_name("iso_upper_mz")
            .and_then(|c| c.as_any().downcast_ref::<Float64Array>());
        let hill_id_col = batch
            .column_by_name("hill_id")
            .and_then(|c| c.as_any().downcast_ref::<arrow::array::UInt64Array>());

        let row_offset = hills.len() as u64;
        for i in 0..batch.num_rows() {
            let profile: Vec<f32> = serde_json::from_str(profile_col.value(i))?;
            let isolation_window = match (iso_target_col, iso_lower_col, iso_upper_col) {
                (Some(t), Some(l), Some(u)) if !t.is_null(i) && !l.is_null(i) && !u.is_null(i) => {
                    Some(IsolationWindow {
                        target: t.value(i),
                        lower: l.value(i),
                        upper: u.value(i),
                    })
                }
                _ => None,
            };
            let hill_id = hill_id_col
                .filter(|c| !c.is_null(i))
                .map(|c| c.value(i))
                .unwrap_or(row_offset + i as u64);
            hills.push(Hill {
                hill_id,
                mz: mz_col.value(i),
                mz_std: mz_std_col.value(i),
                mz_se: mz_se_col.map(|c| c.value(i)).unwrap_or(0.0),
                rt: rt_col.value(i),
                rt_start: rt_start_col.value(i),
                rt_end: rt_end_col.value(i),
                rt_width: rt_width_col.value(i),
                im: im_col.value(i),
                im_std: im_std_col.value(i),
                scan_start: scan_start_col.value(i) as usize,
                scan_apex: scan_apex_col.value(i) as usize,
                scan_end: scan_end_col.value(i) as usize,
                n_scans: n_scans_col.value(i) as usize,
                skipped_scans: skipped_col.value(i) as usize,
                intensity_sum: int_sum_col.value(i),
                intensity_max: int_max_col.value(i),
                hill_score: score_col.value(i),
                intensity_profile: Arc::from(profile.as_slice()),
                isolation_window,
                faims_cv: None,
            });
        }
    }

    Ok(hills)
}

pub fn read_features_parquet(path: &Path) -> Result<Vec<ScoredFeature>, KothError> {
    use arrow::array::{Array, Float64Array, Int64Array, Int8Array, StringArray, UInt8Array};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = std::fs::File::open(path)?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| KothError::ParquetError(e.to_string()))?
        .build()
        .map_err(|e| KothError::ParquetError(e.to_string()))?;

    let mut features = Vec::new();

    for batch in reader {
        let batch = batch.map_err(|e| KothError::ParquetError(e.to_string()))?;

        macro_rules! f64_col {
            ($name:expr) => {
                batch
                    .column_by_name($name)
                    .ok_or_else(|| KothError::ParquetError(format!("missing column: {}", $name)))?
                    .as_any()
                    .downcast_ref::<Float64Array>()
                    .ok_or_else(|| KothError::ParquetError(format!("wrong type: {}", $name)))?
            };
        }

        let mz_col = f64_col!("mz");
        let rt_apex_col = f64_col!("rtApex");
        let rt_start_col = f64_col!("rtStart");
        let rt_end_col = f64_col!("rtEnd");
        let int_apex_col = f64_col!("intensityApex");
        let int_sum_col = f64_col!("intensitySum");
        let im_col = f64_col!("im");
        let cosine_col = f64_col!("cosine_score");
        let ppm_col = f64_col!("ppm_error");
        let isotope_col = f64_col!("isotope_score");
        let combined_col = f64_col!("combined_score");
        let faims_col = batch
            .column_by_name("FAIMS")
            .and_then(|c| c.as_any().downcast_ref::<Float64Array>());

        let charge_col = batch
            .column_by_name("charge")
            .ok_or_else(|| KothError::ParquetError("missing column: charge".into()))?
            .as_any()
            .downcast_ref::<UInt8Array>()
            .ok_or_else(|| KothError::ParquetError("wrong type: charge".into()))?;

        let n_scans_col = batch
            .column_by_name("nScans")
            .ok_or_else(|| KothError::ParquetError("missing column: nScans".into()))?
            .as_any()
            .downcast_ref::<Int64Array>()
            .ok_or_else(|| KothError::ParquetError("wrong type: nScans".into()))?;

        let offset_col = batch
            .column_by_name("neutron_offset")
            .ok_or_else(|| KothError::ParquetError("missing column: neutron_offset".into()))?
            .as_any()
            .downcast_ref::<Int8Array>()
            .ok_or_else(|| KothError::ParquetError("wrong type: neutron_offset".into()))?;

        let theo_col = batch
            .column_by_name("theoretical_pattern")
            .ok_or_else(|| KothError::ParquetError("missing column: theoretical_pattern".into()))?
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| KothError::ParquetError("wrong type: theoretical_pattern".into()))?;

        for i in 0..batch.num_rows() {
            let charge = charge_col.value(i);
            if charge == 0 {
                continue;
            }
            let neutron_offset = offset_col.value(i);
            let mz = mz_col.value(i);
            let raw_mz = mz + neutron_offset as f64 * C13_NEUTRON / charge as f64;
            let theoretical_pattern: Vec<f64> = serde_json::from_str(theo_col.value(i))?;

            let hill = Hill {
                hill_id: 0,
                mz: raw_mz,
                mz_std: 0.0,
                mz_se: 0.0,
                rt: rt_apex_col.value(i),
                rt_start: rt_start_col.value(i),
                rt_end: rt_end_col.value(i),
                rt_width: rt_end_col.value(i) - rt_start_col.value(i),
                im: im_col.value(i),
                im_std: 0.0,
                scan_start: 0,
                scan_apex: 0,
                scan_end: 0,
                n_scans: n_scans_col.value(i) as usize,
                skipped_scans: 0,
                intensity_sum: int_sum_col.value(i),
                intensity_max: int_apex_col.value(i),
                hill_score: 1.0,
                intensity_profile: Arc::from(&[] as &[f32]),
                isolation_window: None,
                faims_cv: faims_col
                    .filter(|c| !c.is_null(i))
                    .map(|c| c.value(i) as f32),
            };

            features.push(ScoredFeature {
                feature: Feature {
                    hills: vec![hill],
                    charge,
                    cosine_score: cosine_col.value(i),
                    ppm_error: ppm_col.value(i),
                },
                neutron_offset,
                isotope_score: isotope_col.value(i),
                cosine_score: cosine_col.value(i),
                combined_score: combined_col.value(i),
                theoretical_pattern,
            });
        }
    }

    Ok(features)
}

// ──────────────────────────────────────────────────────────────────────────────
// Batch discovery
// ──────────────────────────────────────────────────────────────────────────────

/// A discovered run: the name plus paths to its hills and features files.
pub struct RunPaths {
    pub name: String,
    pub hills_path: std::path::PathBuf,
    pub features_path: std::path::PathBuf,
}

/// Scan `batch_dir` for subdirectories that contain both a hills file and a
/// features file.  Returns one `RunPaths` per discovered run, sorted by name.
pub fn discover_runs(batch_dir: &Path) -> Result<Vec<RunPaths>, KothError> {
    let mut runs = Vec::new();

    let entries = std::fs::read_dir(batch_dir)?;
    for entry in entries {
        let entry = entry?;
        let meta = entry.metadata()?;
        if !meta.is_dir() {
            continue;
        }

        let run_dir = entry.path();
        let name = run_dir
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .to_string();

        let hills_path = find_file(&run_dir, "hills")?;
        let features_path = find_file(&run_dir, "features")?;

        match (hills_path, features_path) {
            (Some(h), Some(f)) => runs.push(RunPaths {
                name,
                hills_path: h,
                features_path: f,
            }),
            _ => {
                log::debug!(
                    "Skipping '{}': no hills+features pair found",
                    run_dir.display()
                );
            }
        }
    }

    runs.sort_by(|a, b| a.name.cmp(&b.name));
    Ok(runs)
}

/// Find `<stem>.parquet` or `<stem>.tsv` in a directory (parquet preferred).
fn find_file(dir: &Path, stem: &str) -> Result<Option<std::path::PathBuf>, KothError> {
    for ext in &["parquet", "tsv"] {
        let p = dir.join(format!("{}.{}", stem, ext));
        if p.exists() {
            return Ok(Some(p));
        }
    }
    Ok(None)
}

// ──────────────────────────────────────────────────────────────────────────────
// Ion-mobility scale guard for cross-run consumers
// ──────────────────────────────────────────────────────────────────────────────

pub use koth_core::config::{SCALE_CALIBRATED, SCALE_LINEAR};
/// `report.json` label for input koth did not convert (mzML, Thermo).
pub const SCALE_AS_READ: &str = "as-read-from-input";

/// What one run's outputs say about their 1/K0 scale.
#[derive(Debug, Clone, PartialEq)]
pub struct RunMobility {
    pub name: String,
    /// `mobility.scale` from the run's `report.json`; `None` when the report or
    /// the field is absent (koth 0.9.0 and earlier, or a hand-assembled batch).
    pub recorded: Option<String>,
    /// Whether any feature carries a nonzero ion mobility.
    pub has_mobility: bool,
}

/// The scale a run's `im` values are on, for the purpose of mixing runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EffectiveScale {
    /// No mobility to compare (mzML/Thermo input, or no IM values at all).
    NotApplicable,
    Calibrated,
    Linear {
        recorded: bool,
    },
}

impl RunMobility {
    fn effective(&self) -> Result<EffectiveScale, KothError> {
        match self.recorded.as_deref() {
            Some(SCALE_CALIBRATED) => Ok(EffectiveScale::Calibrated),
            Some(SCALE_LINEAR) => Ok(EffectiveScale::Linear { recorded: true }),
            Some(SCALE_AS_READ) => Ok(EffectiveScale::NotApplicable),
            Some(other) => Err(KothError::ConfigError(format!(
                "run '{}': unknown ion-mobility scale '{other}' in report.json",
                self.name
            ))),
            // Written before the scale was recorded: timsTOF output from koth
            // 0.9.0 or earlier, which only had the linear scale.
            None if self.has_mobility => Ok(EffectiveScale::Linear { recorded: false }),
            None => Ok(EffectiveScale::NotApplicable),
        }
    }
}

/// Read `mobility.scale` from the `report.json` beside a run's features file.
/// A missing report or field is `None`; an unreadable report is an error.
pub fn recorded_mobility_scale(features_path: &Path) -> Result<Option<String>, KothError> {
    let Some(dir) = features_path.parent() else {
        return Ok(None);
    };
    let report = dir.join("report.json");
    if !report.exists() {
        return Ok(None);
    }
    let text = std::fs::read_to_string(&report)?;
    let value: serde_json::Value = serde_json::from_str(&text)
        .map_err(|e| KothError::ConfigError(format!("cannot parse {}: {e}", report.display())))?;
    Ok(value
        .get("mobility")
        .and_then(|m| m.get("scale"))
        .and_then(|s| s.as_str())
        .map(str::to_owned))
}

/// Describe one run's scale from its paths and loaded features.
pub fn run_mobility(
    paths: &RunPaths,
    features: &[ScoredFeature],
) -> Result<RunMobility, KothError> {
    Ok(RunMobility {
        name: paths.name.clone(),
        recorded: recorded_mobility_scale(&paths.features_path)?,
        has_mobility: features.iter().any(|f| f.feature.im_apex() != 0.0),
    })
}

/// Refuse a batch whose timsTOF runs report 1/K0 on different scales.
///
/// Calibrated and linear 1/K0 differ by up to ~0.03 (about 2%), more than the
/// alignment and extraction IM windows, so mixing them silently mis-matches
/// features across runs. Runs without ion mobility are ignored. A timsTOF run
/// with no recorded scale is treated as linear (koth ≤ 0.9.0).
pub fn check_mobility_scales(runs: &[RunMobility]) -> Result<(), KothError> {
    let mut calibrated = Vec::new();
    let mut linear = Vec::new();
    for run in runs {
        match run.effective()? {
            EffectiveScale::NotApplicable => {}
            EffectiveScale::Calibrated => {
                calibrated.push(format!("{} ({SCALE_CALIBRATED})", run.name))
            }
            EffectiveScale::Linear { recorded: true } => {
                linear.push(format!("{} ({SCALE_LINEAR})", run.name))
            }
            EffectiveScale::Linear { recorded: false } => linear.push(format!(
                "{} (no scale recorded: koth 0.9.0 or earlier, linear 1/K0)",
                run.name
            )),
        }
    }
    if calibrated.is_empty() || linear.is_empty() {
        return Ok(());
    }
    Err(KothError::ConfigError(format!(
        "runs report ion mobility on different 1/K0 scales and cannot be aligned together.\n  \
         calibrated: {}\n  linear: {}\n\
         Re-run feature detection on the linear runs with koth 0.10.0 or later (calibrated \
         is the default), or run every input with `[file] bruker_mobility_scale = \"linear\"`.",
        calibrated.join(", "),
        linear.join(", ")
    )))
}

/// [`run_mobility`] for every run, then [`check_mobility_scales`].
pub fn check_batch_mobility_scales<'a>(
    runs: impl IntoIterator<Item = (&'a RunPaths, &'a [ScoredFeature])>,
) -> Result<(), KothError> {
    let described = runs
        .into_iter()
        .map(|(p, f)| run_mobility(p, f))
        .collect::<Result<Vec<_>, _>>()?;
    check_mobility_scales(&described)
}

#[cfg(test)]
mod mobility_scale_tests {
    use super::*;

    fn run(name: &str, recorded: Option<&str>, has_mobility: bool) -> RunMobility {
        RunMobility {
            name: name.into(),
            recorded: recorded.map(str::to_owned),
            has_mobility,
        }
    }

    #[test]
    fn all_calibrated_passes() {
        let runs = [
            run("a", Some(SCALE_CALIBRATED), true),
            run("b", Some(SCALE_CALIBRATED), true),
        ];
        assert!(check_mobility_scales(&runs).is_ok());
    }

    #[test]
    fn all_linear_passes_recorded_or_not() {
        let runs = [
            run("a", Some(SCALE_LINEAR), true),
            run("b", None, true),
            run("c", None, true),
        ];
        assert!(check_mobility_scales(&runs).is_ok());
    }

    #[test]
    fn mixed_scales_fail_and_name_the_runs() {
        let runs = [
            run("new_a", Some(SCALE_CALIBRATED), true),
            run("lin_b", Some(SCALE_LINEAR), true),
        ];
        let msg = check_mobility_scales(&runs).unwrap_err().to_string();
        assert!(
            msg.contains("new_a (bruker-acquisition-calibrated-1/K0)"),
            "{msg}"
        );
        assert!(msg.contains("lin_b (timsrust-linear-1/K0)"), "{msg}");
        assert!(msg.contains("bruker_mobility_scale = \"linear\""), "{msg}");
    }

    #[test]
    fn missing_metadata_is_linear_and_conflicts_with_calibrated() {
        let runs = [
            run("new_a", Some(SCALE_CALIBRATED), true),
            run("old_b", None, true),
        ];
        let msg = check_mobility_scales(&runs).unwrap_err().to_string();
        assert!(msg.contains("old_b (no scale recorded"), "{msg}");
        assert!(msg.contains("Re-run"), "{msg}");
    }

    #[test]
    fn runs_without_mobility_are_ignored() {
        let runs = [
            run("tims", Some(SCALE_CALIBRATED), true),
            run("mzml", Some(SCALE_AS_READ), true),
            run("orbitrap_old", None, false),
            run("tims_lin_no_im", Some(SCALE_AS_READ), false),
        ];
        assert!(check_mobility_scales(&runs).is_ok());
        let only_old = [
            run("orbi", None, false),
            run("tims", Some(SCALE_LINEAR), true),
        ];
        assert!(check_mobility_scales(&only_old).is_ok());
    }

    #[test]
    fn unknown_label_is_an_error() {
        let runs = [run("a", Some("furlongs"), true)];
        assert!(check_mobility_scales(&runs).is_err());
    }

    #[test]
    fn labels_match_the_writer() {
        use crate::io::tims_calibration::MobilityScale;
        assert_eq!(MobilityScale::Calibrated.label(), SCALE_CALIBRATED);
        assert_eq!(MobilityScale::Linear.label(), SCALE_LINEAR);
        assert_eq!(
            crate::output::MobilityReport::as_read().scale,
            SCALE_AS_READ
        );
    }

    #[test]
    fn reads_scale_from_report_json_and_tolerates_absence() {
        let dir =
            std::env::temp_dir().join(format!("koth_mobility_scale_test_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let features = dir.join("features.parquet");
        assert_eq!(recorded_mobility_scale(&features).unwrap(), None);
        std::fs::write(dir.join("report.json"), r#"{"hills":{},"features":{}}"#).unwrap();
        assert_eq!(recorded_mobility_scale(&features).unwrap(), None);
        std::fs::write(
            dir.join("report.json"),
            format!(r#"{{"mobility":{{"scale":"{SCALE_CALIBRATED}"}}}}"#),
        )
        .unwrap();
        assert_eq!(
            recorded_mobility_scale(&features).unwrap().as_deref(),
            Some(SCALE_CALIBRATED)
        );
        std::fs::write(dir.join("report.json"), "not json").unwrap();
        assert!(recorded_mobility_scale(&features).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
