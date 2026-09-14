//! Replay feature detection/scoring on cached hills with one or more configs.
//! This avoids repeating raw-file parsing and hill detection during experiments.
use std::path::PathBuf;
use std::time::Instant;

use anyhow::Context;
use clap::Parser;
use koth_ff::{config::KothConfig, input::read_hills, output::write_features_parquet};

#[derive(Parser)]
struct Args {
    hills: PathBuf,
    #[arg(long, required = true)]
    config: Vec<PathBuf>,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value_t = 8)]
    threads: usize,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;
    let start = Instant::now();
    let hills = read_hills(&args.hills).context("read cached hills")?;
    let read_s = start.elapsed().as_secs_f64();
    for config_path in &args.config {
        let label = config_path
            .file_stem()
            .context("config needs a file stem")?;
        let out = args.output.join(label);
        anyhow::ensure!(!out.exists(), "refusing to overwrite {}", out.display());
        let config = KothConfig::from_toml(config_path)?;
        let start = Instant::now();
        let features = koth_ff::run_features(&hills, &config.features, &config.file)?;
        let feature_s = start.elapsed().as_secs_f64();
        let start = Instant::now();
        let scored = koth_ff::run_scoring(
            &features,
            &config.scoring,
            &config.features,
            config.file.polarity,
        );
        let scoring_s = start.elapsed().as_secs_f64();
        std::fs::create_dir_all(&out)?;
        let start = Instant::now();
        write_features_parquet(&scored, &out.join("features.parquet"), config.file.polarity)?;
        let write_s = start.elapsed().as_secs_f64();
        std::fs::write(out.join("config.toml"), config.to_toml_string()?)?;
        let report = serde_json::json!({
            "hills": hills.len(), "features_before_filter": features.len(),
            "features_after_filter": scored.len(), "read_hills_s": read_s,
            "feature_s": feature_s, "scoring_s": scoring_s, "write_s": write_s,
            "input": args.hills, "config": config_path, "threads": args.threads,
        });
        std::fs::write(
            out.join("timing.json"),
            serde_json::to_string_pretty(&report)?,
        )?;
        println!(
            "{}: {} features, detection {:.2}s, scoring {:.2}s, writing {:.2}s",
            label.to_string_lossy(),
            scored.len(),
            feature_s,
            scoring_s,
            write_s
        );
    }
    Ok(())
}
