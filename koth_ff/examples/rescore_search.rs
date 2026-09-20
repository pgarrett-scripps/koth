//! Diagnostic replay of exported, rounded LFQ features. Native extraction is
//! required for final validation; exported ppm/RT/features have finite precision.
use anyhow::{Context, Result};
use koth_ff::lfq::{
    rescore::search::{compute, SearchScoring},
    LfqEntry,
};
use std::{
    collections::{BTreeMap, HashMap},
    path::PathBuf,
};

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    anyhow::ensure!(
        args.len() == 3,
        "usage: rescore_search INPUT_SEARCH_FOLDER NEW_OUTPUT.tsv"
    );
    let input = PathBuf::from(&args[1]);
    let output = PathBuf::from(&args[2]);
    anyhow::ensure!(!output.exists(), "refusing to overwrite replay output");
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(input.join("peptide_quant.tsv"))?;
    let headers = reader.headers()?.clone();
    let col = |name: &str| headers.iter().position(|h| h == name).unwrap();
    let (fi, pi, ci, ri, ei) = (
        col("target_id"),
        col("modified_peptide"),
        col("charge"),
        col("run"),
        col("evidence"),
    );
    let mut peptides = Vec::new();
    let mut charges = Vec::new();
    let mut runs = BTreeMap::new();
    let mut direct = HashMap::new();
    for result in reader.records() {
        let row = result?;
        let i: usize = row[fi].parse()?;
        if i >= peptides.len() {
            peptides.resize(i + 1, String::new());
            charges.resize(i + 1, 0);
        }
        peptides[i] = row[pi].to_owned();
        charges[i] = row[ci].parse()?;
        let next = runs.len();
        let r = *runs.entry(row[ri].to_owned()).or_insert(next);
        direct.insert((i, r), &row[ei] == "direct_ms2");
    }
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(input.join("lfq_details.tsv"))?;
    let header = reader.headers()?.clone();
    let positions: HashMap<_, _> = header
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_owned(), i))
        .collect();
    let mut entries = Vec::new();
    let mut exported = HashMap::new();
    for result in reader.records() {
        let row = result?;
        let get = |key: &str| -> &str { &row[positions[key]] };
        let f = |key: &str| get(key).parse::<f64>().unwrap_or(f64::NAN);
        let feature_idx = get("feature_idx").parse::<usize>()?;
        let run_idx = *runs.get(get("run_name")).context("unknown run")?;
        let is_decoy = get("is_decoy") == "true";
        // Recreate the five model inputs from exported residuals, avoiding
        // confusion between reference m/z and run/decoy-specific expected m/z.
        entries.push(LfqEntry {
            feature_idx,
            run_idx,
            is_decoy,
            is_mbr: !direct[&(feature_idx, run_idx)],
            intensity: f("intensity"),
            hybrid_score: f("hybrid_score") as f32,
            spectral_bhattacharyya: f("spectral_bhattacharyya") as f32,
            coelution: f("coelution") as f32,
            rt_score: f("rt_score") as f32,
            int_score: f("int_score") as f32,
            n_isotopes_found: get("n_isotopes_found").parse()?,
            expected_mz: 500.0,
            observed_mz: 500.0 * (1.0 + f("ppm_error") / 1e6),
            expected_rt: 0.0,
            apex_rt: f("rt_diff"),
            peak_width_rt: f("peak_width_rt"),
            expected_im: 1.0,
            observed_im: 1.0 + f("im_delta"),
            owned_samples: 0,
            excluded_samples: 0,
            competing_feature: None,
            ownership_status: "replay",
            preceding_signal_fraction: f("preceding_signal_fraction") as f32,
        });
        if !is_decoy {
            exported.insert((feature_idx, run_idx), f("q_value"));
        }
    }
    let mut maps = Vec::new();
    for mode in [
        SearchScoring::Legacy,
        SearchScoring::PeptideGrouped,
        SearchScoring::ChargeStratified,
        SearchScoring::ChargeSigned,
        SearchScoring::ChargeSignedQuality,
    ] {
        let mut q = HashMap::new();
        for is_direct in [false, true] {
            let subset: Vec<_> = entries
                .iter()
                .filter(|e| direct[&(e.feature_idx, e.run_idx)] == is_direct)
                .cloned()
                .collect();
            q.extend(compute(&subset, &peptides, &charges, mode));
        }
        eprintln!(
            "{mode:?}: {} accepted target cells",
            q.values().filter(|&&v| v <= 0.01).count()
        );
        maps.push(q);
    }
    let mut names = vec![String::new(); runs.len()];
    for (name, i) in runs {
        names[i] = name;
    }
    let mut w = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output)?;
    w.write_record([
        "target_id",
        "modified_peptide",
        "charge",
        "run",
        "evidence",
        "intensity",
        "exported_q",
        "replayed_legacy_q",
        "peptide_grouped_q",
        "charge_stratified_q",
        "charge_signed_q",
        "charge_signed_quality_q",
    ])?;
    for e in entries.iter().filter(|e| !e.is_decoy) {
        let key = (e.feature_idx, e.run_idx);
        w.write_record([
            e.feature_idx.to_string(),
            peptides[e.feature_idx].clone(),
            charges[e.feature_idx].to_string(),
            names[e.run_idx].clone(),
            if direct[&key] { "direct_ms2" } else { "mbr" }.into(),
            e.intensity.to_string(),
            exported[&key].to_string(),
            maps[0][&key].to_string(),
            maps[1][&key].to_string(),
            maps[2][&key].to_string(),
            maps[3][&key].to_string(),
            maps[4][&key].to_string(),
        ])?;
    }
    w.flush()?;
    Ok(())
}
