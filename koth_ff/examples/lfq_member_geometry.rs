//! Inspect aligned original observations in selected cached candidate groups.
use koth_ff::{
    alignment::{align_runs, RunInput},
    config::AlignConfig,
    input::{check_batch_mobility_scales, discover_runs, read_features},
};
use std::{collections::HashSet, path::PathBuf};
fn main() -> anyhow::Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    anyhow::ensure!(
        args.len() == 5,
        "usage: lfq_member_geometry BATCH CONFIG MEMBERS_TSV OUT_TSV ID,ID"
    );
    let cfg = AlignConfig::from_toml(&PathBuf::from(&args[1]))?;
    let paths = discover_runs(&PathBuf::from(&args[0]))?;
    let runs: Vec<_> = paths
        .iter()
        .map(|r| {
            Ok(RunInput {
                name: r.name.clone(),
                features: read_features(&r.features_path)?,
                hills: Vec::new(),
                scan_times: Vec::new(),
                rt_bounds: None,
            })
        })
        .collect::<anyhow::Result<_>>()?;
    check_batch_mobility_scales(
        paths
            .iter()
            .zip(&runs)
            .map(|(p, r)| (p, r.features.as_slice())),
    )?;
    let alignment = align_runs(&runs, &cfg.alignment);
    let ids: HashSet<usize> = args[4]
        .split(',')
        .map(str::parse)
        .collect::<Result<_, _>>()?;
    let mut input = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&args[2])?;
    let mut output = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&args[3])?;
    output.write_record([
        "candidate_id",
        "run_id",
        "feature_idx",
        "mass",
        "rt",
        "im",
        "quality",
        "is_primary",
        "neutron_offset",
    ])?;
    for row in input.records() {
        let row = row?;
        let id: usize = row[0].parse()?;
        if !ids.contains(&id) {
            continue;
        }
        let r: usize = row[1].parse()?;
        let f: usize = row[2].parse()?;
        let feat = &runs[r].features[f];
        let mz = feat.monoisotopic_mz();
        let rt = feat.feature.rt_apex();
        let im = feat.feature.im_apex();
        let (mz, rt, im) = if r == alignment.reference_idx {
            (mz, rt, im)
        } else {
            let a = &alignment.alignments[&runs[r].name];
            (a.correct_mz(mz, rt), a.warp_rt(rt), a.correct_im(im, rt))
        };
        let mass = (mz - 1.007276466621) * feat.feature.charge as f64;
        output.write_record([
            id.to_string(),
            r.to_string(),
            f.to_string(),
            mass.to_string(),
            rt.to_string(),
            im.to_string(),
            feat.combined_score.to_string(),
            row[3].to_string(),
            feat.neutron_offset.to_string(),
        ])?;
    }
    output.flush()?;
    Ok(())
}
