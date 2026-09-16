//! Versioned long evidence export. Original descriptors come from the source
//! table, not reconstructed single-hill placeholders used during alignment.
use crate::{
    alignment::{AlignmentResult, RunInput},
    config::AlignConfig,
    lfq::IntensityMatrix,
};
use anyhow::{ensure, Context, Result};
use sha2::{Digest, Sha256};
use std::{
    collections::HashMap,
    fs::File,
    io::Read,
    path::{Path, PathBuf},
};

pub fn sha256(path: &Path) -> Result<String> {
    let mut f = File::open(path)?;
    let mut h = Sha256::new();
    let mut buf = [0; 65536];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(format!("{:x}", h.finalize()))
}

const OBS: &[&str] = &[
    "run_id",
    "observation_id",
    "consensus_id",
    "is_primary",
    "is_seed",
    "exclusion_reason",
    "feature_mass",
    "feature_mz",
    "feature_charge",
    "feature_rt",
    "feature_im",
    "feature_combined_score",
    "feature_n_isotopes",
    "feature_n_scans",
    "feature_cosine_score",
    "feature_isotope_score",
    "feature_ppm_error",
    "feature_neutron_offset",
    "feature_rt_width",
    "feature_isoerror",
    "feature_isoerror2",
    "feature_intensity_apex",
    "feature_intensity_sum",
    "aligned_mass",
    "aligned_rt",
    "aligned_im",
];

/// Stream TSV or Parquet rows without retaining elution profiles in memory.
pub(crate) fn source_rows(
    path: &Path,
    mut visit: impl FnMut(&csv::StringRecord, &csv::StringRecord) -> Result<()>,
) -> Result<()> {
    if path.extension().and_then(|s| s.to_str()) == Some("parquet") {
        use arrow::{array::Array, util::display::array_value_to_string};
        let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(
            File::open(path)?,
        )?
        .build()?;
        for batch in reader {
            let batch = batch?;
            let header = csv::StringRecord::from(
                batch
                    .schema()
                    .fields()
                    .iter()
                    .map(|f| f.name().as_str())
                    .collect::<Vec<_>>(),
            );
            for i in 0..batch.num_rows() {
                let row = csv::StringRecord::from(
                    batch
                        .columns()
                        .iter()
                        .map(|a| {
                            if a.is_null(i) {
                                Ok(String::new())
                            } else {
                                array_value_to_string(a.as_ref(), i)
                            }
                        })
                        .collect::<std::result::Result<Vec<_>, _>>()?,
                );
                visit(&header, &row)?;
            }
        }
    } else {
        let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
        let header = reader.headers()?.clone();
        for row in reader.records() {
            visit(&header, &row?)?;
        }
    }
    Ok(())
}

pub fn write_bundle(
    matrix: &IntensityMatrix,
    runs: &[RunInput],
    alignment: &AlignmentResult,
    sources: &[PathBuf],
    config: &AlignConfig,
    out: &Path,
) -> Result<()> {
    ensure!(
        sources.len() == runs.len(),
        "feature source/run count mismatch"
    );
    let manifest_path = out.join("matrix_manifest.json");
    // Never leave an old complete manifest pointing at partially replaced files.
    if manifest_path.exists() {
        std::fs::remove_file(&manifest_path)?;
    }
    let obs_path = out.join("feature_observations.tsv");
    let cells_path = out.join("lfq_matrix.long.tsv");
    let mut membership = HashMap::new();
    for (c, g) in matrix.consensus.iter().enumerate() {
        for &(run, obs) in &g.members {
            ensure!(
                membership.insert((run, obs as usize), c).is_none(),
                "observation belongs to multiple groups"
            );
        }
    }
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&obs_path)?;
    writer.write_record(OBS)?;
    let mut run_manifest = Vec::new();
    // The alignment reader skips charge-zero rows. Map its indices back to
    // source row IDs, retaining those skipped rows in the observation table.
    let mut source_ids: Vec<Vec<usize>> = Vec::new();
    for (run, path) in sources.iter().enumerate() {
        let hash = sha256(path)?;
        let mut row_id = 0;
        let mut ids = Vec::new();
        source_rows(path, |header, row| {
            let get = |names: &[&str]| -> String {
                names
                    .iter()
                    .find_map(|n| {
                        header
                            .iter()
                            .position(|h| h == *n)
                            .map(|i| row.get(i).unwrap_or("").to_owned())
                    })
                    .unwrap_or_default()
            };
            let charge: u8 = get(&["charge"]).parse()?;
            let fi = ids.len();
            let group = if charge == 0 {
                None
            } else {
                membership.get(&(run, fi)).copied()
            };
            let primary =
                group.is_some_and(|c| matrix.consensus[c].per_run_feature[run] == Some(fi as u32));
            let seed = group.is_some_and(|c| {
                matrix.consensus[c].seed_run_idx == run
                    && matrix.consensus[c].seed_feature_idx == fi as u32
            });
            let mz: f64 = get(&["mz"]).parse()?;
            let rt: f64 = get(&["rtApex"]).parse()?;
            let im = get(&["im"]);
            let im_num = if im.is_empty() { 0.0 } else { im.parse()? };
            let (am, ar, ai) = if run == alignment.reference_idx {
                (mz, rt, im_num)
            } else {
                let a = &alignment.alignments[&runs[run].name];
                (
                    a.correct_mz(mz, rt),
                    a.warp_rt(rt),
                    a.correct_im(im_num, rt),
                )
            };
            let reason = if charge == 0 {
                "unknown_charge"
            } else if group.is_some() {
                ""
            } else if !get(&["combined_score"]).parse::<f64>()?.is_finite() {
                "invalid_feature_score"
            } else {
                "group_confidence_or_singleton"
            };
            let mut record = vec![
                run.to_string(),
                row_id.to_string(),
                group.map(|x| x.to_string()).unwrap_or_default(),
                primary.to_string(),
                seed.to_string(),
                reason.into(),
            ];
            for names in [
                vec!["massCalib"],
                vec!["mz"],
                vec!["charge"],
                vec!["rtApex"],
                vec!["im"],
                vec!["combined_score"],
                vec!["nIsotopes", "n_isotopes"],
                vec!["nScans", "n_scans"],
                vec!["cosine_score"],
                vec!["isotope_score"],
                vec!["ppm_error"],
                vec!["neutron_offset"],
                vec!["rt_width"],
                vec!["isoerror"],
                vec!["isoerror2"],
                vec!["intensityApex"],
                vec!["intensitySum"],
            ] {
                record.push(get(&names));
            }
            record.extend([
                (am * charge as f64 - charge as f64 * 1.007_276_466_621).to_string(),
                ar.to_string(),
                if im.is_empty() {
                    String::new()
                } else {
                    ai.to_string()
                },
            ]);
            writer.write_record(record)?;
            if charge != 0 {
                ids.push(row_id);
            }
            row_id += 1;
            Ok(())
        })
        .with_context(|| format!("exporting original observations from {}", path.display()))?;
        ensure!(
            ids.len() == runs[run].features.len(),
            "source/loaded feature count mismatch"
        );
        ensure!(
            hash == sha256(path)?,
            "feature source changed during export"
        );
        run_manifest.push(serde_json::json!({"run_id":run,"name":runs[run].name,"source":path,"sha256":hash,"n_observations":row_id}));
        source_ids.push(ids);
    }
    writer.flush()?;
    let mut cells = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&cells_path)?;
    cells.write_record([
        "consensus_id",
        "run_id",
        "extraction_kind",
        "observation_id",
        "seed_run_id",
        "seed_observation_id",
        "mass",
        "mz",
        "charge",
        "rt",
        "im",
        "combined_score",
        "n_contributing_runs",
        "group_score",
        "group_qvalue",
        "lfq_intensity",
        "lfq_q_value",
        "lfq_is_mbr",
        "lfq_hybrid_score",
        "lfq_isotope_pattern",
        "lfq_coelution",
        "lfq_n_isotopes_found",
        "lfq_expected_rt",
        "lfq_apex_rt",
        "lfq_expected_mz",
        "lfq_observed_mz",
        "lfq_status",
        "lfq_ownership_status",
        "lfq_owned_samples",
        "lfq_excluded_samples",
        "lfq_competing_consensus_id",
        "lfq_preceding_signal_fraction",
    ])?;
    let finite = |x: f64| {
        if x.is_finite() {
            x.to_string()
        } else {
            String::new()
        }
    };
    for e in &matrix.entries {
        let c = e.feature_idx;
        let r = e.run_idx;
        let g = &matrix.consensus[c];
        let original = if e.is_decoy {
            None
        } else {
            g.per_run_feature[r].map(|i| source_ids[r][i as usize])
        };
        let intensity = if e.is_decoy {
            matrix.decoy_intensity(c, r)
        } else {
            matrix.intensity(c, r)
        };
        cells.write_record(vec![
            c.to_string(),
            r.to_string(),
            if e.is_decoy { "decoy" } else { "target" }.into(),
            original.map(|x| x.to_string()).unwrap_or_default(),
            g.seed_run_idx.to_string(),
            if matrix.search_guidance.is_some() {
                String::new()
            } else {
                source_ids[g.seed_run_idx][g.seed_feature_idx as usize].to_string()
            },
            g.neutral_mass.to_string(),
            g.ref_mz.to_string(),
            g.charge.to_string(),
            g.ref_rt.to_string(),
            if g.ref_im == 0.0 {
                String::new()
            } else {
                g.ref_im.to_string()
            },
            finite(g.seed_combined_score),
            g.n_contributing_runs.to_string(),
            finite(g.group_score),
            finite(g.group_qvalue),
            intensity.to_string(),
            if e.is_decoy || !config.lfq.run_tdc {
                String::new()
            } else {
                matrix.q_value(c, r).to_string()
            },
            e.is_mbr.to_string(),
            finite(e.hybrid_score as f64),
            finite(e.spectral_bhattacharyya as f64),
            finite(e.coelution as f64),
            e.n_isotopes_found.to_string(),
            finite(e.expected_rt),
            finite(e.apex_rt),
            finite(e.expected_mz),
            finite(e.observed_mz),
            if intensity > 0.0 {
                "signal"
            } else {
                "no_signal"
            }
            .into(),
            e.ownership_status.into(),
            e.owned_samples.to_string(),
            e.excluded_samples.to_string(),
            e.competing_feature
                .map(|i| i.to_string())
                .unwrap_or_default(),
            e.preceding_signal_fraction.to_string(),
        ])?;
    }
    cells.flush()?;
    let manifest = serde_json::json!({"schema_version":if matrix.search_guidance.is_some() { 4 } else { 3 },"search_guided":matrix.search_guidance.is_some(),"n_consensus":matrix.n_features,"reference_run_id":alignment.reference_idx,
        "runs":run_manifest,"config":config,"units":{"mass":"Da","rt":"minutes","im":"1/K0"},
        "observations":{"path":"feature_observations.tsv","sha256":sha256(&obs_path)?},
        "cells":{"path":"lfq_matrix.long.tsv","sha256":sha256(&cells_path)?}});
    let temp = out.join("matrix_manifest.json.tmp");
    std::fs::write(&temp, serde_json::to_vec_pretty(&manifest)?)?;
    std::fs::rename(temp, manifest_path)?;
    Ok(())
}
