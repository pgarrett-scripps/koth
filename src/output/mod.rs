use std::path::Path;

use crate::error::KothError;
use crate::models::{Hill, ScoredFeature};

/// Write hills to a TSV file.
///
/// Column order matches the Python zenith_feature_finder hills.tsv exactly.
/// The `intensity_profile` column is serialized as a JSON array.
pub fn write_hills_tsv(hills: &[Hill], path: &Path) -> Result<(), KothError> {
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
    wtr.write_record([
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
        "intensity_profile",
    ])?;

    for &i in &indices {
        let h = &hills[i];
        let profile_json = serde_json::to_string(h.intensity_profile.as_ref())?;
        wtr.write_record([
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
            profile_json,
        ])?;
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
        "cosine_similarity",
        "ppm_error",
        "neutron_offset",
        "score",
        "theoretical_pattern",
        "isotope_profile",
        "elution_profile",
        "mono_hills_scan_lists",
        "mono_hills_intensity_list",
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
        let mono_scans: Vec<usize> = f.mono_scan_list();
        let mono_intensities: Vec<f64> = f.mono_intensity_list().iter().map(|&x| x as f64).collect();

        let theo_json = serde_json::to_string(&sf.theoretical_pattern)?;
        let isotope_json = serde_json::to_string(&isotope_profile)?;
        let elution_json = serde_json::to_string(&elution_vec)?;
        let mono_scans_json = serde_json::to_string(&mono_scans)?;
        let mono_int_json = serde_json::to_string(&mono_intensities)?;

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
            format!("{:.6}", f.cosine_similarity),
            format!("{:.6}", f.ppm_error),
            sf.neutron_offset.to_string(),
            format!("{:.6}", sf.score),
            theo_json,
            isotope_json,
            elution_json,
            mono_scans_json,
            mono_int_json,
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
