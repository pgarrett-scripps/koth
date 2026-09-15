//! Controlled scoring/control ablations; intentionally outside the production CLI.
#[path = "support/lfq_matrices.rs"]
mod lfq_matrices;
use clap::Parser;
use koth_ff::{
    alignment::{align_runs, RunInput},
    config::AlignConfig,
    input::{discover_runs, read_features, read_hills},
    lfq::{
        consensus::{audit_consensus, VARIANT_NAMES},
        quantify_consensus,
    },
    output::{build_align_report, write_align_report, AlignTiming},
};
use lfq_matrices::matrices;
use std::{io::Write, path::PathBuf, time::Instant};

#[derive(Parser)]
struct Args {
    batch: PathBuf,
    config: PathBuf,
    output: PathBuf,
    /// Evenly spaced run subset (original names are preserved).
    #[arg(long)]
    take: Option<usize>,
    /// Comma-separated variant indices to quantify. Empty = audit only.
    #[arg(long, default_value = "")]
    quantify: String,
    /// Only the first five global controls; full audit defaults to 25 controls.
    #[arg(long)]
    quick: bool,
    /// Use the research quality-conditioned controls for extraction.
    #[arg(long)]
    conditioned: bool,
}

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let args = Args::parse();
    anyhow::ensure!(
        !(args.quick && args.conditioned),
        "conditioned controls require the full audit"
    );
    let config = AlignConfig::from_toml(&args.config)?;
    let all_paths = discover_runs(&args.batch)?;
    anyhow::ensure!(!all_paths.is_empty(), "no runs");
    let count = args.take.unwrap_or(all_paths.len());
    anyhow::ensure!(
        count >= 2 && count <= all_paths.len(),
        "take must be 2..=number of runs"
    );
    let selected: Vec<_> = (0..count).map(|i| i * all_paths.len() / count).collect();
    let paths: Vec<_> = selected.iter().map(|&i| &all_paths[i]).collect();
    let mut runs = Vec::new();
    for p in &paths {
        runs.push(RunInput {
            name: p.name.clone(),
            features: read_features(&p.features_path)?,
            hills: Vec::new(),
            scan_times: Vec::new(),
        });
    }
    if let Some(name) = &config.alignment.reference_run {
        anyhow::ensure!(
            runs.iter().any(|r| &r.name == name),
            "reference run is not selected"
        );
    }
    std::fs::create_dir_all(&args.output)?;
    std::fs::write(args.output.join("runs.json"),serde_json::to_vec_pretty(&paths.iter().enumerate().map(|(id,p)|
        serde_json::json!({"run_id":id,"name":p.name,"features":p.features_path})).collect::<Vec<_>>())?)?;
    let t = Instant::now();
    let alignment = align_runs(&runs, &config.alignment);
    let alignment_sec = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let audit = audit_consensus(&runs, &alignment, &config.lfq, !args.quick);
    let audit_sec = t.elapsed().as_secs_f64();
    std::fs::write(
        args.output.join("audit.json"),
        serde_json::to_vec_pretty(&serde_json::json!({
        "variants":VARIANT_NAMES,"controls":audit.null_group_counts,"full_controls":!args.quick,"alignment_sec":alignment_sec,"audit_sec":audit_sec,
        "n_runs":runs.len(),"n_groups":audit.groups.len(),"reference":alignment.reference_name}))?,
    )?;
    let mut out = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(args.output.join("candidates.tsv"))?;
    let mut header: Vec<String> = [
        "candidate_id",
        "massCalib",
        "mz",
        "charge",
        "rtApex",
        "im",
        "combined_score",
        "n_contributing_runs",
        "weak_members",
    ]
    .map(String::from)
    .to_vec();
    for v in VARIANT_NAMES {
        for field in [
            "score",
            "global",
            "repeat",
            "density",
            "pooled",
            "conditioned",
        ] {
            header.push(format!("{v}_{field}"));
        }
    }
    out.write_record(header)?;
    let mut members = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(args.output.join("members.tsv"))?;
    members.write_record(["candidate_id", "run_id", "feature_idx", "is_primary"])?;
    for (i, g) in audit.groups.iter().enumerate() {
        let weak = g
            .members
            .iter()
            .filter(|&&(r, f)| runs[r].features[f as usize].combined_score < 0.5)
            .count();
        let mut row = vec![
            i.to_string(),
            g.neutral_mass.to_string(),
            g.ref_mz.to_string(),
            g.charge.to_string(),
            g.ref_rt.to_string(),
            g.ref_im.to_string(),
            g.seed_combined_score.to_string(),
            g.n_contributing_runs.to_string(),
            weak.to_string(),
        ];
        for v in 0..4 {
            row.extend(
                [
                    audit.scores[i][v],
                    audit.q_global[i][v],
                    audit.q_repeat[i][v],
                    audit.q_density[i][v],
                    audit.q_pooled[i][v],
                    audit.q_conditioned[i][v],
                ]
                .map(|x| x.to_string()),
            );
        }
        out.write_record(row)?;
        for &(r, f) in &g.members {
            members.write_record([
                i.to_string(),
                r.to_string(),
                f.to_string(),
                (g.per_run_feature[r] == Some(f)).to_string(),
            ])?;
        }
    }
    out.flush()?;
    members.flush()?;
    let selected_q = if args.conditioned {
        &audit.q_conditioned
    } else {
        &audit.q_global
    };
    for variant in args.quantify.split(',').filter(|s| !s.is_empty()) {
        let v: usize = variant.parse()?;
        anyhow::ensure!(v < 4, "variant must be 0..3");
        let folder = args.output.join(VARIANT_NAMES[v]);
        std::fs::create_dir_all(&folder)?;
        let groups = audit
            .groups
            .iter()
            .enumerate()
            .filter(|(i, g)| {
                g.n_contributing_runs >= 2
                    && selected_q[*i][v] <= config.lfq.consensus.max_group_qvalue
            })
            .map(|(i, g)| {
                let mut g = g.clone();
                g.group_score = audit.scores[i][v];
                g.group_qvalue = selected_q[i][v];
                g
            })
            .collect();
        let t = Instant::now();
        let matrix = quantify_consensus(&runs, &alignment, &config.lfq, groups, |i| {
            read_hills(&paths[i].hills_path).expect("read hills")
        });
        let lfq_sec = t.elapsed().as_secs_f64();
        matrices(&matrix, &folder)?;
        std::fs::write(folder.join("align_config.toml"), config.to_toml_string()?)?;
        let report = build_align_report(
            &runs,
            &alignment,
            &matrix,
            &config,
            AlignTiming {
                alignment_sec,
                lfq_sec,
                total_sec: alignment_sec + audit_sec + lfq_sec,
            },
        );
        write_align_report(&report, &folder.join("align_report.json"))?;
        let mut meta = std::fs::File::create(folder.join("variant.json"))?;
        writeln!(
            meta,
            "{}",
            serde_json::json!({"variant":VARIANT_NAMES[v],"group_gate":config.lfq.consensus.max_group_qvalue,"control": if args.conditioned { "ten quality-conditioned rank shifts" } else { "first five global shifts" }})
        )?;
    }
    Ok(())
}
