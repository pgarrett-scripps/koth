//! Export all candidate groups (including rejected groups) without loading hills.
//! cargo run --release --example lfq_group_audit -- BATCH CONFIG OUT
use koth_ff::{
    alignment::{align_runs, RunInput},
    config::AlignConfig,
    input::{discover_runs, read_features},
    lfq::consensus::build_consensus_candidates,
};
use std::{io::Write, path::PathBuf};

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().skip(1).map(PathBuf::from).collect();
    anyhow::ensure!(args.len() == 3, "usage: lfq_group_audit BATCH CONFIG OUT");
    let cfg = AlignConfig::from_toml(&args[1])?;
    let paths = discover_runs(&args[0])?;
    anyhow::ensure!(!paths.is_empty(), "no feature runs found");
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
    let alignment = align_runs(&runs, &cfg.alignment);
    let groups = build_consensus_candidates(&runs, &alignment, &cfg.lfq);
    std::fs::create_dir_all(&args[2])?;
    let mut out = std::io::BufWriter::new(std::fs::File::create(args[2].join("candidates.tsv"))?);
    writeln!(out, "candidate_id\tmassCalib\tmz\tcharge\trtApex\tim\tcombined_score\tn_contributing_runs\tgroup_score\tgroup_qvalue\tweak_members\tretained")?;
    for (i, g) in groups.iter().enumerate() {
        let weak = g
            .members
            .iter()
            .filter(|&&(r, f)| runs[r].features[f as usize].combined_score < 0.5)
            .count();
        writeln!(
            out,
            "{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{weak}\t{}",
            g.neutral_mass,
            g.ref_mz,
            g.charge,
            g.ref_rt,
            g.ref_im,
            g.seed_combined_score,
            g.n_contributing_runs,
            g.group_score,
            g.group_qvalue,
            g.n_contributing_runs >= 2 && g.group_qvalue <= cfg.lfq.consensus.max_group_qvalue
        )?;
    }
    out.flush()?;
    Ok(())
}
