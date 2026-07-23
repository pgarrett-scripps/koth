pub mod align_report;
pub mod report;
pub use align_report::{build_align_report, write_align_report, AlignReport, Timing as AlignTiming};
pub use report::{build_features_report, build_hills_report, write_report, RunReport};

use std::path::Path;

use crate::error::KothError;
use crate::models::{Hill, ScoredFeature};

/// Hill output columns, in emission order. Shared by the TSV header and the
/// Parquet schema so the two backends can't drift apart. The Parquet element
/// types are supplied alongside at the schema-construction site (all base
/// columns are non-nullable).
const HILL_COLUMNS: &[&str] = &[
    "hill_id",
    "mz",
    "mz_std",
    "mz_se",
    "rt",
    "rt_start",
    "rt_end",
    "rt_width",
    "im",
    "im_std",
    "scan_start",
    "scan_apex",
    "scan_end",
    "n_scans",
    "skipped_scans",
    "intensity_sum",
    "intensity_max",
    "hill_score",
    "intensity_profile",
];

/// Extra hill columns appended when the isolation window is included (MS2).
const HILL_ISO_COLUMNS: &[&str] = &["iso_target_mz", "iso_lower_mz", "iso_upper_mz"];

/// Scored-feature output columns, in emission order. Shared by the TSV header
/// and the Parquet schema. `massCalib` is nullable in Parquet; the remaining
/// columns are non-nullable (types/nullability supplied at the schema site).
const FEATURE_COLUMNS: &[&str] = &[
    "massCalib",
    "mz",
    "rtApex",
    "rtStart",
    "rtEnd",
    "intensityApex",
    "intensitySum",
    "charge",
    "nIsotopes",
    "nScans",
    "im",
    "cosine_score",
    "ppm_error",
    "neutron_offset",
    "isotope_score",
    "combined_score",
    "theoretical_pattern",
    "isotope_profile",
    "elution_profile",
    "hill_ids",
    "intensityApexParab",
    "intensityScattered5",
    "intensityConsec5",
];

/// Alternative quant estimators derived from a feature's per-scan elution
/// profile (the total-across-isotopes intensity at each scan — the same source
/// as `intensityApex` = max and `intensitySum` = sum). Returns
/// `(apex_parabolic, scattered5, consec5)`:
///   * `apex_parabolic` — peak height refined by 3-point parabolic
///     interpolation around the max scan (standard MS centroiding). Removes
///     scan-grid-alignment jitter from the raw apex; falls back to the raw
///     apex when there is no concave-down interior maximum.
///   * `scattered5` — sum of the 5 highest scans (or all scans if fewer). A
///     robust peak-top area that ignores noisy flanks/tails and dodges dips.
///   * `consec5` — max sum over contiguous windows of length min(5, n). The
///     contiguous-window analogue of `scattered5`.
/// These accompany the conventional `intensityApex` / `intensitySum` so the
/// LFQ layer can compare quant readouts (see SI estimator comparison).
fn quant_variants(profile: &[f64]) -> (f64, f64, f64) {
    if profile.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let n = profile.len();

    // argmax
    let mut i_max = 0;
    for (i, &v) in profile.iter().enumerate() {
        if v > profile[i_max] {
            i_max = i;
        }
    }

    // parabolic-refined apex (vertex of the parabola through the 3 points
    // around the max). Only when the fit is concave-down (denom < 0).
    let apex_parab = if i_max > 0 && i_max + 1 < n {
        let (y0, y1, y2) = (profile[i_max - 1], profile[i_max], profile[i_max + 1]);
        let denom = y0 - 2.0 * y1 + y2;
        if denom < 0.0 {
            y1 - 0.125 * (y2 - y0).powi(2) / denom
        } else {
            y1
        }
    } else {
        profile[i_max]
    };

    // scattered top-5: sum of the 5 largest scans.
    let k = n.min(5);
    let mut sorted = profile.to_vec();
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let scattered5: f64 = sorted[..k].iter().sum();

    // consecutive-5: max sum over contiguous windows of length min(5, n),
    // via a sliding window.
    let w = n.min(5);
    let mut window: f64 = profile[..w].iter().sum();
    let mut consec5 = window;
    for j in w..n {
        window += profile[j] - profile[j - w];
        if window > consec5 {
            consec5 = window;
        }
    }

    (apex_parab, scattered5, consec5)
}

/// Write hills to a TSV file.
///
/// Column order matches the Python zenith_feature_finder hills.tsv exactly.
/// The `intensity_profile` column is serialized as a JSON array.
/// When `include_isolation_window` is true, three extra columns are appended:
/// `iso_target_mz`, `iso_lower_mz`, `iso_upper_mz`.
pub fn write_hills_tsv(hills: &[Hill], path: &Path) -> Result<(), KothError> {
    write_hills_tsv_inner(hills, path, false)
}

/// Write MS2 hills (one row per hill, with isolation window columns).
pub fn write_ms2_hills_tsv(hills: &[Hill], path: &Path) -> Result<(), KothError> {
    write_hills_tsv_inner(hills, path, true)
}

fn write_hills_tsv_inner(
    hills: &[Hill],
    path: &Path,
    include_isolation_window: bool,
) -> Result<(), KothError> {
    // Sort by intensity_sum descending (matches Python output)
    let mut indices: Vec<usize> = (0..hills.len()).collect();
    indices.sort_by(|&a, &b| {
        hills[b]
            .intensity_sum
            .partial_cmp(&hills[a].intensity_sum)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)?;

    // Header
    let mut header: Vec<&str> = HILL_COLUMNS.to_vec();
    if include_isolation_window {
        header.extend_from_slice(HILL_ISO_COLUMNS);
    }
    wtr.write_record(&header)?;

    for &i in &indices {
        let h = &hills[i];
        let profile_json = serde_json::to_string(h.intensity_profile.as_ref())?;
        let mut row = vec![
            h.hill_id.to_string(),
            format!("{:.6}", h.mz),
            format!("{:.6}", h.mz_std),
            format!("{:.6e}", h.mz_se),
            format!("{:.6}", h.rt),
            format!("{:.6}", h.rt_start),
            format!("{:.6}", h.rt_end),
            format!("{:.6}", h.rt_width),
            format!("{:.6}", h.im),
            format!("{:.6}", h.im_std),
            h.scan_start.to_string(),
            h.scan_apex.to_string(),
            h.scan_end.to_string(),
            h.n_scans.to_string(),
            h.skipped_scans.to_string(),
            format!("{:.5e}", h.intensity_sum),
            format!("{:.5e}", h.intensity_max),
            format!("{:.6}", h.hill_score),
            profile_json,
        ];
        if include_isolation_window {
            match h.isolation_window {
                Some(iw) => {
                    row.push(format!("{:.6}", iw.target));
                    row.push(format!("{:.6}", iw.lower));
                    row.push(format!("{:.6}", iw.upper));
                }
                None => {
                    row.push(String::new());
                    row.push(String::new());
                    row.push(String::new());
                }
            }
        }
        wtr.write_record(&row)?;
    }

    wtr.flush()?;
    log::info!("Wrote {} hills to {}", hills.len(), path.display());
    Ok(())
}

/// Write scored features to a TSV file.
///
/// Column order matches the Python zenith_feature_finder features.tsv exactly.
/// List columns (elution_profile, isotope_profile, etc.) are JSON arrays.
pub fn write_features_tsv(features: &[ScoredFeature], path: &Path) -> Result<(), KothError> {
    // Sort by intensitySum descending, exclude charge=0
    let mut scored: Vec<&ScoredFeature> = features
        .iter()
        .filter(|f| f.feature.charge > 0)
        .collect();
    scored.sort_by(|a, b| {
        b.feature
            .total_intensity()
            .partial_cmp(&a.feature.total_intensity())
            .unwrap()
    });

    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)?;

    // Header
    wtr.write_record(FEATURE_COLUMNS)?;

    for sf in &scored {
        let f = &sf.feature;

        let mass_calib = sf
            .monoisotopic_neutral_mass()
            .map(|m| format!("{:.6}", m))
            .unwrap_or_default();

        let mono_mz = sf.monoisotopic_mz();
        let rt_apex = f.rt_apex();
        let rt_start = f.rt_start();
        let rt_end = f.rt_end();
        let intensity_apex = f.total_intensity_at_apex();
        let intensity_sum = f.total_intensity();
        let n_isotopes = f.hills.len();
        let n_scans = f.n_scans_total();
        let im = f.im_apex();

        let (min_scan, max_scan, elution_vec) = f.elution_profile();
        let _ = (min_scan, max_scan);
        let (apex_parab, scattered5, consec5) = quant_variants(&elution_vec);
        let isotope_profile = f.isotope_profile_apex();

        let theo_json = serde_json::to_string(&sf.theoretical_pattern)?;
        let isotope_json = serde_json::to_string(&isotope_profile)?;
        let elution_json = serde_json::to_string(&elution_vec)?;
        let hill_ids: Vec<u64> = f.hills.iter().map(|h| h.hill_id).collect();
        let hill_ids_json = serde_json::to_string(&hill_ids)?;

        let im_str = if im != 0.0 {
            format!("{:.6}", im)
        } else {
            String::new()
        };

        wtr.write_record([
            mass_calib,
            format!("{:.6}", mono_mz),
            format!("{:.6}", rt_apex),
            format!("{:.6}", rt_start),
            format!("{:.6}", rt_end),
            format!("{:.5e}", intensity_apex),
            format!("{:.5e}", intensity_sum),
            f.charge.to_string(),
            n_isotopes.to_string(),
            n_scans.to_string(),
            im_str,
            format!("{:.6}", f.cosine_score),
            format!("{:.6}", f.ppm_error),
            sf.neutron_offset.to_string(),
            format!("{:.6}", sf.isotope_score),
            format!("{:.6}", sf.combined_score),
            theo_json,
            isotope_json,
            elution_json,
            hill_ids_json,
            format!("{:.5e}", apex_parab),
            format!("{:.5e}", scattered5),
            format!("{:.5e}", consec5),
        ])?;
    }

    wtr.flush()?;
    log::info!(
        "Wrote {} features to {}",
        scored.len(),
        path.display()
    );
    Ok(())
}

/// Write hills to a Parquet file.
///
/// Schema mirrors hills.tsv; `intensity_profile` is stored as a JSON string.
pub fn write_hills_parquet(hills: &[Hill], path: &Path) -> Result<(), KothError> {
    write_hills_parquet_inner(hills, path, false)
}

/// Write MS2 hills to a Parquet file (extra `iso_*_mz` columns).
pub fn write_ms2_hills_parquet(hills: &[Hill], path: &Path) -> Result<(), KothError> {
    write_hills_parquet_inner(hills, path, true)
}

fn write_hills_parquet_inner(
    hills: &[Hill],
    path: &Path,
    include_isolation_window: bool,
) -> Result<(), KothError> {
    use std::sync::Arc;
    use arrow::array::{ArrayRef, Float64Array, Int64Array, StringArray, UInt64Array};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;

    let mut indices: Vec<usize> = (0..hills.len()).collect();
    indices.sort_by(|&a, &b| {
        hills[b].intensity_sum.partial_cmp(&hills[a].intensity_sum).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Parquet element types, one per `HILL_COLUMNS` entry (same order). Names
    // come from the shared const so the TSV header and Parquet schema stay in
    // lock-step; all base columns are non-nullable.
    let hill_types: [DataType; 19] = [
        DataType::UInt64,  // hill_id
        DataType::Float64, // mz
        DataType::Float64, // mz_std
        DataType::Float64, // mz_se
        DataType::Float64, // rt
        DataType::Float64, // rt_start
        DataType::Float64, // rt_end
        DataType::Float64, // rt_width
        DataType::Float64, // im
        DataType::Float64, // im_std
        DataType::Int64,   // scan_start
        DataType::Int64,   // scan_apex
        DataType::Int64,   // scan_end
        DataType::Int64,   // n_scans
        DataType::Int64,   // skipped_scans
        DataType::Float64, // intensity_sum
        DataType::Float64, // intensity_max
        DataType::Float64, // hill_score
        DataType::Utf8,    // intensity_profile
    ];
    let mut fields: Vec<Field> = HILL_COLUMNS
        .iter()
        .zip(hill_types)
        .map(|(name, dt)| Field::new(*name, dt, false))
        .collect();
    if include_isolation_window {
        for name in HILL_ISO_COLUMNS {
            fields.push(Field::new(*name, DataType::Float64, true));
        }
    }
    let schema = Arc::new(Schema::new(fields));

    let profiles: Result<Vec<String>, serde_json::Error> = indices.iter()
        .map(|&i| serde_json::to_string(hills[i].intensity_profile.as_ref()))
        .collect();
    let profiles = profiles?;

    let mut columns: Vec<ArrayRef> = vec![
        Arc::new(indices.iter().map(|&i| hills[i].hill_id).collect::<UInt64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].mz).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].mz_std).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].mz_se).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].rt).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].rt_start).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].rt_end).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].rt_width).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].im).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].im_std).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].scan_start as i64).collect::<Int64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].scan_apex as i64).collect::<Int64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].scan_end as i64).collect::<Int64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].n_scans as i64).collect::<Int64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].skipped_scans as i64).collect::<Int64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].intensity_sum).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].intensity_max).collect::<Float64Array>()),
        Arc::new(indices.iter().map(|&i| hills[i].hill_score).collect::<Float64Array>()),
        Arc::new(profiles.iter().map(|s| Some(s.as_str())).collect::<StringArray>()),
    ];
    if include_isolation_window {
        let target: Float64Array = indices
            .iter()
            .map(|&i| hills[i].isolation_window.map(|iw| iw.target))
            .collect();
        let lower: Float64Array = indices
            .iter()
            .map(|&i| hills[i].isolation_window.map(|iw| iw.lower))
            .collect();
        let upper: Float64Array = indices
            .iter()
            .map(|&i| hills[i].isolation_window.map(|iw| iw.upper))
            .collect();
        columns.push(Arc::new(target));
        columns.push(Arc::new(lower));
        columns.push(Arc::new(upper));
    }

    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| KothError::ParquetError(e.to_string()))?;

    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, None)
        .map_err(|e| KothError::ParquetError(e.to_string()))?;
    writer.write(&batch).map_err(|e| KothError::ParquetError(e.to_string()))?;
    writer.close().map_err(|e| KothError::ParquetError(e.to_string()))?;

    log::info!("Wrote {} hills to {}", hills.len(), path.display());
    Ok(())
}

/// Write scored features to a Parquet file.
///
/// Schema mirrors features.tsv; list columns are stored as JSON strings.
pub fn write_features_parquet(features: &[ScoredFeature], path: &Path) -> Result<(), KothError> {
    use std::sync::Arc;
    use arrow::array::{ArrayRef, Float64Array, Int8Array, Int64Array, UInt8Array, StringArray};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;

    let mut scored: Vec<&ScoredFeature> = features
        .iter()
        .filter(|f| f.feature.charge > 0)
        .collect();
    scored.sort_by(|a, b| {
        b.feature.total_intensity().partial_cmp(&a.feature.total_intensity()).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Parquet element type + nullability, one per `FEATURE_COLUMNS` entry (same
    // order). Names come from the shared const so the TSV header and Parquet
    // schema stay in lock-step.
    let feature_types: [(DataType, bool); 23] = [
        (DataType::Float64, true),  // massCalib
        (DataType::Float64, false), // mz
        (DataType::Float64, false), // rtApex
        (DataType::Float64, false), // rtStart
        (DataType::Float64, false), // rtEnd
        (DataType::Float64, false), // intensityApex
        (DataType::Float64, false), // intensitySum
        (DataType::UInt8, false),   // charge
        (DataType::Int64, false),   // nIsotopes
        (DataType::Int64, false),   // nScans
        (DataType::Float64, false), // im
        (DataType::Float64, false), // cosine_score
        (DataType::Float64, false), // ppm_error
        (DataType::Int8, false),    // neutron_offset
        (DataType::Float64, false), // isotope_score
        (DataType::Float64, false), // combined_score
        (DataType::Utf8, false),    // theoretical_pattern
        (DataType::Utf8, false),    // isotope_profile
        (DataType::Utf8, false),    // elution_profile
        (DataType::Utf8, false),    // hill_ids
        (DataType::Float64, false), // intensityApexParab
        (DataType::Float64, false), // intensityScattered5
        (DataType::Float64, false), // intensityConsec5
    ];
    let schema = Arc::new(Schema::new(
        FEATURE_COLUMNS
            .iter()
            .zip(feature_types)
            .map(|(name, (dt, nullable))| Field::new(*name, dt, nullable))
            .collect::<Vec<Field>>(),
    ));

    let n = scored.len();
    let mut mass_calib       = Vec::<Option<f64>>::with_capacity(n);
    let mut mz               = Vec::<f64>::with_capacity(n);
    let mut rt_apex          = Vec::<f64>::with_capacity(n);
    let mut rt_start         = Vec::<f64>::with_capacity(n);
    let mut rt_end           = Vec::<f64>::with_capacity(n);
    let mut int_apex         = Vec::<f64>::with_capacity(n);
    let mut int_sum          = Vec::<f64>::with_capacity(n);
    let mut charge           = Vec::<u8>::with_capacity(n);
    let mut n_isotopes       = Vec::<i64>::with_capacity(n);
    let mut n_scans          = Vec::<i64>::with_capacity(n);
    let mut im               = Vec::<f64>::with_capacity(n);
    let mut cosine_score     = Vec::<f64>::with_capacity(n);
    let mut ppm_err          = Vec::<f64>::with_capacity(n);
    let mut neutron_offset   = Vec::<i8>::with_capacity(n);
    let mut isotope_score    = Vec::<f64>::with_capacity(n);
    let mut combined_score   = Vec::<f64>::with_capacity(n);
    let mut theo_json        = Vec::<String>::with_capacity(n);
    let mut iso_json         = Vec::<String>::with_capacity(n);
    let mut elut_json        = Vec::<String>::with_capacity(n);
    let mut hill_ids_json    = Vec::<String>::with_capacity(n);
    let mut apex_parab_col   = Vec::<f64>::with_capacity(n);
    let mut scattered5_col   = Vec::<f64>::with_capacity(n);
    let mut consec5_col      = Vec::<f64>::with_capacity(n);

    for sf in &scored {
        let f = &sf.feature;
        mass_calib.push(sf.monoisotopic_neutral_mass());
        mz.push(sf.monoisotopic_mz());
        rt_apex.push(f.rt_apex());
        rt_start.push(f.rt_start());
        rt_end.push(f.rt_end());
        int_apex.push(f.total_intensity_at_apex());
        int_sum.push(f.total_intensity());
        charge.push(f.charge);
        n_isotopes.push(f.hills.len() as i64);
        n_scans.push(f.n_scans_total() as i64);
        im.push(f.im_apex());
        cosine_score.push(f.cosine_score);
        ppm_err.push(f.ppm_error);
        neutron_offset.push(sf.neutron_offset);
        isotope_score.push(sf.isotope_score);
        combined_score.push(sf.combined_score);

        let (_, _, elution_vec) = f.elution_profile();
        let iso_profile = f.isotope_profile_apex();
        let (apex_parab, scattered5, consec5) = quant_variants(&elution_vec);
        apex_parab_col.push(apex_parab);
        scattered5_col.push(scattered5);
        consec5_col.push(consec5);

        theo_json.push(serde_json::to_string(&sf.theoretical_pattern)?);
        iso_json.push(serde_json::to_string(&iso_profile)?);
        elut_json.push(serde_json::to_string(&elution_vec)?);

        let hill_ids: Vec<u64> = f.hills.iter().map(|h| h.hill_id).collect();
        hill_ids_json.push(serde_json::to_string(&hill_ids)?);
    }

    let columns: Vec<ArrayRef> = vec![
        Arc::new(mass_calib.into_iter().collect::<Float64Array>()),
        Arc::new(mz.into_iter().collect::<Float64Array>()),
        Arc::new(rt_apex.into_iter().collect::<Float64Array>()),
        Arc::new(rt_start.into_iter().collect::<Float64Array>()),
        Arc::new(rt_end.into_iter().collect::<Float64Array>()),
        Arc::new(int_apex.into_iter().collect::<Float64Array>()),
        Arc::new(int_sum.into_iter().collect::<Float64Array>()),
        Arc::new(charge.into_iter().collect::<UInt8Array>()),
        Arc::new(n_isotopes.into_iter().collect::<Int64Array>()),
        Arc::new(n_scans.into_iter().collect::<Int64Array>()),
        Arc::new(im.into_iter().collect::<Float64Array>()),
        Arc::new(cosine_score.into_iter().collect::<Float64Array>()),
        Arc::new(ppm_err.into_iter().collect::<Float64Array>()),
        Arc::new(neutron_offset.into_iter().collect::<Int8Array>()),
        Arc::new(isotope_score.into_iter().collect::<Float64Array>()),
        Arc::new(combined_score.into_iter().collect::<Float64Array>()),
        Arc::new(theo_json.iter().map(|s| Some(s.as_str())).collect::<StringArray>()),
        Arc::new(iso_json.iter().map(|s| Some(s.as_str())).collect::<StringArray>()),
        Arc::new(elut_json.iter().map(|s| Some(s.as_str())).collect::<StringArray>()),
        Arc::new(hill_ids_json.iter().map(|s| Some(s.as_str())).collect::<StringArray>()),
        Arc::new(apex_parab_col.into_iter().collect::<Float64Array>()),
        Arc::new(scattered5_col.into_iter().collect::<Float64Array>()),
        Arc::new(consec5_col.into_iter().collect::<Float64Array>()),
    ];

    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| KothError::ParquetError(e.to_string()))?;

    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, None)
        .map_err(|e| KothError::ParquetError(e.to_string()))?;
    writer.write(&batch).map_err(|e| KothError::ParquetError(e.to_string()))?;
    writer.close().map_err(|e| KothError::ParquetError(e.to_string()))?;

    log::info!("Wrote {} features to {}", scored.len(), path.display());
    Ok(())
}
