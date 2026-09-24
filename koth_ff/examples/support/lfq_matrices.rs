use koth_ms::lfq::IntensityMatrix;
use std::path::Path;

pub fn matrices(m: &IntensityMatrix, folder: &Path) -> anyhow::Result<()> {
    let meta = [
        "massCalib",
        "mz",
        "charge",
        "rtApex",
        "im",
        "combined_score",
        "seed_run",
        "n_contributing_runs",
    ];
    for (kind, name) in [
        "consensus_features.tsv",
        "intensity_matrix.tsv",
        "qvalue_matrix.tsv",
        "decoy_intensity_matrix.tsv",
    ]
    .iter()
    .enumerate()
    {
        let mut out = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(folder.join(name))?;
        let mut header: Vec<String> = meta.iter().map(|s| s.to_string()).collect();
        if kind == 0 {
            header.extend(["n_runs_detected", "group_score", "group_qvalue"].map(String::from));
        } else {
            header.extend(m.run_names.clone());
        }
        out.write_record(header)?;
        for (i, g) in m.consensus.iter().enumerate() {
            let mut row = vec![
                format!("{:.6}", g.neutral_mass),
                format!("{:.6}", g.ref_mz),
                g.charge.to_string(),
                format!("{:.6}", g.ref_rt),
                if g.ref_im == 0.0 {
                    String::new()
                } else {
                    format!("{:.6}", g.ref_im)
                },
                format!("{:.6}", g.seed_combined_score),
                m.run_names[g.seed_run_idx].clone(),
                g.n_contributing_runs.to_string(),
            ];
            if kind == 0 {
                row.extend([
                    (0..m.n_runs)
                        .filter(|&r| m.intensity(i, r) > 0.0)
                        .count()
                        .to_string(),
                    format!("{:.6}", g.group_score),
                    format!("{:.6}", g.group_qvalue),
                ]);
            } else {
                for r in 0..m.n_runs {
                    let v = match kind {
                        1 => m.intensity(i, r),
                        2 => m.q_value(i, r),
                        _ => m.decoy_intensity(i, r),
                    };
                    row.push(if kind == 2 {
                        format!("{v:.4}")
                    } else if v > 0.0 {
                        format!("{v:.5e}")
                    } else {
                        "0".into()
                    });
                }
            }
            out.write_record(row)?;
        }
        out.flush()?;
    }
    let mut audit = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(folder.join("ownership.tsv"))?;
    audit.write_record([
        "feature_idx",
        "run_idx",
        "is_decoy",
        "intensity",
        "ownership_status",
        "owned_samples",
        "excluded_samples",
        "competing_feature",
        "apex_rt",
        "bhattacharyya",
        "coelution",
        "preceding_signal_fraction",
    ])?;
    for e in &m.entries {
        audit.write_record([
            e.feature_idx.to_string(),
            e.run_idx.to_string(),
            e.is_decoy.to_string(),
            e.intensity.to_string(),
            e.ownership_status.to_string(),
            e.owned_samples.to_string(),
            e.excluded_samples.to_string(),
            e.competing_feature
                .map(|i| i.to_string())
                .unwrap_or_default(),
            e.apex_rt.to_string(),
            e.spectral_bhattacharyya.to_string(),
            e.coelution.to_string(),
            e.preceding_signal_fraction.to_string(),
        ])?;
    }
    audit.flush()?;
    Ok(())
}
