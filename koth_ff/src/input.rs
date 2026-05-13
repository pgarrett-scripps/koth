/// Readers for hill and feature files produced by `koth_ff`.
///
/// Supports both TSV (default) and Parquet output formats.
/// Auto-detection is based on the file extension.
use std::path::Path;
use std::sync::Arc;

use serde::Deserialize;

use crate::error::KothError;
use crate::models::{Feature, Hill, ScoredFeature};

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
    mz: f64,
    mz_std: f64,
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
}

pub fn read_hills_tsv(path: &Path) -> Result<Vec<Hill>, KothError> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let mut hills = Vec::new();
    for result in rdr.deserialize::<HillRow>() {
        let row = result?;
        let profile: Vec<f32> = serde_json::from_str(&row.intensity_profile)?;
        hills.push(Hill {
            mz: row.mz,
            mz_std: row.mz_std,
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
    cosine_similarity: f64,
    ppm_error: f64,
    neutron_offset: i8,
    score: f64,
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
        mz: raw_mz,
        mz_std: 0.0,
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
    };

    ScoredFeature {
        feature: Feature {
            hills: vec![hill],
            charge: row.charge,
            cosine_similarity: row.cosine_similarity,
            ppm_error: row.ppm_error,
        },
        neutron_offset: row.neutron_offset,
        score: row.score,
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
    use arrow::array::{Float64Array, Int64Array, StringArray};
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

        for i in 0..batch.num_rows() {
            let profile: Vec<f32> = serde_json::from_str(profile_col.value(i))?;
            hills.push(Hill {
                mz: mz_col.value(i),
                mz_std: mz_std_col.value(i),
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
            });
        }
    }

    Ok(hills)
}

pub fn read_features_parquet(path: &Path) -> Result<Vec<ScoredFeature>, KothError> {
    use arrow::array::{Float64Array, Int8Array, Int64Array, UInt8Array, StringArray};
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
        let cosine_col = f64_col!("cosine_similarity");
        let ppm_col = f64_col!("ppm_error");
        let score_col = f64_col!("score");

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
                mz: raw_mz,
                mz_std: 0.0,
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
            };

            features.push(ScoredFeature {
                feature: Feature {
                    hills: vec![hill],
                    charge,
                    cosine_similarity: cosine_col.value(i),
                    ppm_error: ppm_col.value(i),
                },
                neutron_offset,
                score: score_col.value(i),
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
