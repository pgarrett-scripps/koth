pub mod align_report;
pub mod report;
pub use align_report::{build_align_report, write_align_report, AlignReport, Timing as AlignTiming};
pub use report::{build_features_report, build_hills_report, write_report, RunReport};

use std::path::Path;

use crate::error::KothError;
use crate::models::{Hill, ScoredFeature};

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
            .unwrap()
    });

    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)?;

    // Header
    let mut header: Vec<&str> = vec![
        "hill_id",
        "mz",
        "mz_std",
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
    if include_isolation_window {
        header.extend_from_slice(&["iso_target_mz", "iso_lower_mz", "iso_upper_mz"]);
    }
    wtr.write_record(&header)?;

    for &i in &indices {
        let h = &hills[i];
        let profile_json = serde_json::to_string(h.intensity_profile.as_ref())?;
        let mut row = vec![
            h.hill_id.to_string(),
            format!("{:.6}", h.mz),
            format!("{:.6}", h.mz_std),
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
    wtr.write_record([
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
    ])?;

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
        hills[b].intensity_sum.partial_cmp(&hills[a].intensity_sum).unwrap()
    });

    let mut fields: Vec<Field> = vec![
        Field::new("hill_id",         DataType::UInt64,  false),
        Field::new("mz",              DataType::Float64, false),
        Field::new("mz_std",          DataType::Float64, false),
        Field::new("rt",              DataType::Float64, false),
        Field::new("rt_start",        DataType::Float64, false),
        Field::new("rt_end",          DataType::Float64, false),
        Field::new("rt_width",        DataType::Float64, false),
        Field::new("im",              DataType::Float64, false),
        Field::new("im_std",          DataType::Float64, false),
        Field::new("scan_start",      DataType::Int64,   false),
        Field::new("scan_apex",       DataType::Int64,   false),
        Field::new("scan_end",        DataType::Int64,   false),
        Field::new("n_scans",         DataType::Int64,   false),
        Field::new("skipped_scans",   DataType::Int64,   false),
        Field::new("intensity_sum",   DataType::Float64, false),
        Field::new("intensity_max",   DataType::Float64, false),
        Field::new("hill_score",      DataType::Float64, false),
        Field::new("intensity_profile", DataType::Utf8,  false),
    ];
    if include_isolation_window {
        fields.push(Field::new("iso_target_mz", DataType::Float64, true));
        fields.push(Field::new("iso_lower_mz",  DataType::Float64, true));
        fields.push(Field::new("iso_upper_mz",  DataType::Float64, true));
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
        b.feature.total_intensity().partial_cmp(&a.feature.total_intensity()).unwrap()
    });

    let schema = Arc::new(Schema::new(vec![
        Field::new("massCalib",                DataType::Float64, true),
        Field::new("mz",                       DataType::Float64, false),
        Field::new("rtApex",                   DataType::Float64, false),
        Field::new("rtStart",                  DataType::Float64, false),
        Field::new("rtEnd",                    DataType::Float64, false),
        Field::new("intensityApex",            DataType::Float64, false),
        Field::new("intensitySum",             DataType::Float64, false),
        Field::new("charge",                   DataType::UInt8,   false),
        Field::new("nIsotopes",                DataType::Int64,   false),
        Field::new("nScans",                   DataType::Int64,   false),
        Field::new("im",                       DataType::Float64, false),
        Field::new("cosine_score",             DataType::Float64, false),
        Field::new("ppm_error",                DataType::Float64, false),
        Field::new("neutron_offset",           DataType::Int8,    false),
        Field::new("isotope_score",            DataType::Float64, false),
        Field::new("combined_score",           DataType::Float64, false),
        Field::new("theoretical_pattern",      DataType::Utf8,    false),
        Field::new("isotope_profile",          DataType::Utf8,    false),
        Field::new("elution_profile",          DataType::Utf8,    false),
        Field::new("hill_ids",                 DataType::Utf8,    false),
    ]));

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
