//! Independent-permutation controls; intentionally outside the production CLI.
#[path = "support/lfq_matrices.rs"]
mod lfq_matrices;
use clap::Parser;
use koth_ms::{
    alignment::{align_runs, RunInput},
    config::AlignConfig,
    input::{check_batch_mobility_scales, discover_runs, read_features, read_hills},
    lfq::{consensus::audit_permutation_consensus, quantify_consensus},
};
use std::{path::PathBuf, time::Instant};

#[derive(Parser)]
struct Args {
    batch: PathBuf,
    config: PathBuf,
    output: PathBuf,
    #[arg(long)]
    take: Option<usize>,
    /// Extract the pooled-ten accepted set with the unchanged production pipeline.
    #[arg(long)]
    quantify: bool,
}

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let args = Args::parse();
    let cfg = AlignConfig::from_toml(&args.config)?;
    let all = discover_runs(&args.batch)?;
    let n = args.take.unwrap_or(all.len());
    anyhow::ensure!(n >= 2 && n <= all.len(), "take must be 2..=number of runs");
    let paths: Vec<_> = (0..n).map(|i| &all[i * all.len() / n]).collect();
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
            .map(|(p, r)| (*p, r.features.as_slice())),
    )?;
    if let Some(name) = &cfg.alignment.reference_run {
        anyhow::ensure!(
            runs.iter().any(|r| &r.name == name),
            "reference run is not selected"
        );
    }
    std::fs::create_dir_all(&args.output)?;
    std::fs::write(args.output.join("runs.json"),serde_json::to_vec_pretty(&paths.iter().enumerate().map(|(id,p)|
        serde_json::json!({"run_id":id,"name":p.name,"features":p.features_path})).collect::<Vec<_>>())?)?;
    let t = Instant::now();
    let alignment = align_runs(&runs, &cfg.alignment);
    let audit = audit_permutation_consensus(&runs, &alignment, &cfg.lfq);
    let mut out = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(args.output.join("candidates.tsv"))?;
    out.write_record([
        "candidate_id",
        "massCalib",
        "mz",
        "charge",
        "rtApex",
        "im",
        "combined_score",
        "n_contributing_runs",
        "weak_members",
        "average_score",
        "average_first",
        "average_repeat",
        "average_permuted",
    ])?;
    let mut counts = [0usize; 3];
    let (mut intersection, mut union) = (0usize, 0usize);
    for (i, g) in audit.groups.iter().enumerate() {
        let keep = [audit.q_first[i], audit.q_repeat[i], audit.q_pooled[i]]
            .map(|q| g.n_contributing_runs >= 2 && q <= cfg.lfq.consensus.max_group_qvalue);
        for j in 0..3 {
            counts[j] += usize::from(keep[j]);
        }
        intersection += usize::from(keep[0] && keep[1]);
        union += usize::from(keep[0] || keep[1]);
        let weak = g
            .members
            .iter()
            .filter(|&&(r, f)| runs[r].features[f as usize].combined_score < 0.5)
            .count();
        out.write_record([
            i.to_string(),
            g.neutral_mass.to_string(),
            g.ref_mz.to_string(),
            g.charge.to_string(),
            g.ref_rt.to_string(),
            g.ref_im.to_string(),
            g.seed_combined_score.to_string(),
            g.n_contributing_runs.to_string(),
            weak.to_string(),
            g.group_score.to_string(),
            audit.q_first[i].to_string(),
            audit.q_repeat[i].to_string(),
            audit.q_pooled[i].to_string(),
        ])?;
    }
    out.flush()?;
    let summary = serde_json::json!({"method":"independent RT permutation within run/charge/quality stratum","n_runs":n,"reference":alignment.reference_name,"elapsed_sec":t.elapsed().as_secs_f64(),"controls":audit.null_group_counts,"retained_first":counts[0],"retained_repeat":counts[1],"retained_pooled":counts[2],"retained_jaccard":if union>0 {Some(intersection as f64 / union as f64)}else{None},"gate":cfg.lfq.consensus.max_group_qvalue});
    std::fs::write(
        args.output.join("audit.json"),
        serde_json::to_vec_pretty(&summary)?,
    )?;
    println!("{summary}");
    if args.quantify {
        let groups = audit
            .groups
            .into_iter()
            .filter(|g| {
                g.n_contributing_runs >= 2 && g.group_qvalue <= cfg.lfq.consensus.max_group_qvalue
            })
            .collect();
        let matrix = quantify_consensus(&runs, &alignment, &cfg.lfq, groups, |i| {
            read_hills(&paths[i].hills_path).expect("read hills")
        });
        let folder = args.output.join("average");
        std::fs::create_dir_all(&folder)?;
        lfq_matrices::matrices(&matrix, &folder)?;
        std::fs::write(folder.join("align_config.toml"), cfg.to_toml_string()?)?;
        std::fs::write(
            folder.join("variant.json"),
            serde_json::to_vec_pretty(&summary)?,
        )?;
    }
    Ok(())
}
