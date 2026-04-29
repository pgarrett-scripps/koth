use std::path::PathBuf;
use std::time::Instant;

use anyhow::Context;
use clap::Parser;

use koth_ff::{
    config::KothConfig,
    mem::log_mem,
    output::{write_features_tsv, write_hills_tsv},
    run_features, run_hills_streaming, run_scoring,
};

#[derive(Parser, Debug)]
#[command(
    name = "koth_ff",
    about = "High-performance LC-MS feature finder for timsTOF and mzML data",
    version
)]
struct Args {
    /// Input file path (.mzML) or Bruker .d directory
    input: PathBuf,

    /// Output directory (default: current directory)
    #[arg(short, long, default_value = ".")]
    output: PathBuf,

    /// TOML configuration file (uses defaults if not provided)
    #[arg(short, long)]
    config: Option<PathBuf>,

    /// Skip the isotope pattern scoring stage
    #[arg(long)]
    no_scoring: bool,

    /// Log level: error, warn, info, debug, trace (default: info)
    #[arg(long, default_value = "info")]
    log_level: String,

}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let log_filter = format!("koth_ff={}", args.log_level);
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(&log_filter))
        .init();

    let config = match &args.config {
        Some(cfg_path) => KothConfig::from_toml(cfg_path)
            .with_context(|| format!("Failed to load config from {}", cfg_path.display()))?,
        None => KothConfig::default(),
    };

    let file_stem = args
        .input
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output");
    let out_dir = args.output.join(file_stem);
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("Failed to create output directory {}", out_dir.display()))?;

    log::info!("Input: {}", args.input.display());
    log_mem("startup");

    log::info!("Output directory: {}", out_dir.display());

    let total_start = Instant::now();

    // Stage 1: Streaming hill detection (no Vec<Spectrum> held in memory)
    log::info!("Detecting hills (streaming)...");
    let t = Instant::now();
    let hills = run_hills_streaming(&args.input, &config.hills)
        .with_context(|| format!("Failed to detect hills from {}", args.input.display()))?;
    log::info!("[timing] hill detection: {:.2?} ({} hills)", t.elapsed(), hills.len());
    log_mem(&format!("after hill detection ({} hills)", hills.len()));

    let hills_path = out_dir.join("hills.tsv");
    let t = Instant::now();
    write_hills_tsv(&hills, &hills_path)
        .with_context(|| format!("Failed to write hills to {}", hills_path.display()))?;
    log::info!("[timing] write hills.tsv: {:.2?}", t.elapsed());

    // Stage 2: Feature detection — hills borrowed here, features clone Hill structs
    // (intensity_profile is Arc so no data duplication)
    log::info!("Detecting features...");
    let t = Instant::now();
    let features = run_features(&hills, &config.features)
        .context("Feature detection failed")?;
    log::info!("[timing] feature detection: {:.2?} ({} features)", t.elapsed(), features.len());

    // Drop hills now — profile data stays alive via Arc refs inside features
    drop(hills);
    log_mem(&format!("after feature detection ({} features)", features.len()));

    // Stage 3: Scoring (optional)
    let features_path = out_dir.join("features.tsv");
    if args.no_scoring {
        let unscored: Vec<_> = features
            .iter()
            .map(|f| koth_ff::models::ScoredFeature {
                feature: f.clone(),
                neutron_offset: 0,
                score: 0.0,
                theoretical_pattern: Vec::new(),
            })
            .collect();
        let t = Instant::now();
        write_features_tsv(&unscored, &features_path)
            .with_context(|| format!("Failed to write features to {}", features_path.display()))?;
        log::info!("[timing] write features.tsv: {:.2?}", t.elapsed());
    } else {
        log::info!("Scoring features...");
        let t = Instant::now();
        let scored = run_scoring(&features, &config.scoring);
        log::info!("[timing] scoring: {:.2?}", t.elapsed());
        let t = Instant::now();
        write_features_tsv(&scored, &features_path)
            .with_context(|| format!("Failed to write features to {}", features_path.display()))?;
        log::info!("[timing] write features.tsv: {:.2?}", t.elapsed());
    }

    log::info!("[timing] total: {:.2?}", total_start.elapsed());

    let config_path = out_dir.join("config.toml");
    let config_str = config.to_toml_string().context("Failed to serialize config")?;
    std::fs::write(&config_path, config_str)
        .with_context(|| format!("Failed to write config to {}", config_path.display()))?;

    log::info!("Done. Results in {}", out_dir.display());
    Ok(())
}
