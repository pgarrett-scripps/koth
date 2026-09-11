use std::path::PathBuf;
use std::time::Instant;

use anyhow::Context;
use clap::Parser;

use koth_ff::{
    config::{KothConfig, OutputFormat},
    mem::log_mem,
    output::{
        build_features_report, build_hills_report, write_features_parquet, write_features_tsv,
        write_hills_parquet, write_hills_tsv, write_ms2_hills_parquet, write_ms2_hills_tsv,
        write_report, RunReport,
    },
    run_features, run_hills_streaming, run_ms2_hills_streaming, run_scoring,
};

#[derive(Parser, Debug)]
#[command(
    name = "koth_ff",
    about = "High-performance LC-MS feature finder for timsTOF and mzML data",
    version = koth_ff::VERSION
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

    /// Also detect MS2 hills (per isolation window) and write them to
    /// `hills_ms2.{tsv|parquet}`. Overrides the `file.ms2_hills_enabled`
    /// config setting when set.
    #[arg(long)]
    ms2: bool,

    /// Enable ID-free isotope-consistency m/z recalibration (a pass-1 feature
    /// pass learns a per-region m/z-offset surface applied in pass 2, and the
    /// isotope-match tolerance becomes region-adaptive). Overrides
    /// `file.mz_recalibration` when set.
    #[arg(long)]
    recalibrate: bool,

    /// Drop large baseline-like hills (≥40 scans without a clear apex).
    /// Overrides `hills.filter_large_baseline_hills` when set. Targets
    /// column bleed / solvent ions / plasticizer contamination.
    #[arg(long)]
    filter_baseline_hills: bool,

    /// Log level: error, warn, info, debug, trace (default: info)
    #[arg(long, default_value = "info")]
    log_level: String,

    /// Worker threads for the parallel stages. Overrides `file.n_threads`
    /// (default when both are unset: all available cores)
    #[arg(long)]
    threads: Option<usize>,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let log_filter = format!("koth_ff={}", args.log_level);
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(&log_filter)).init();

    let mut config = match &args.config {
        Some(cfg_path) => KothConfig::from_toml(cfg_path)
            .with_context(|| format!("Failed to load config from {}", cfg_path.display()))?,
        None => KothConfig::default(),
    };
    if args.ms2 {
        config.file.ms2_hills_enabled = true;
    }
    if args.recalibrate {
        config.file.mz_recalibration = true;
    }
    if args.filter_baseline_hills {
        config.hills.filter_large_baseline_hills = true;
    }
    if args.threads.is_some() {
        config.file.n_threads = args.threads;
    }
    // Must run before the first rayon use anywhere: build_global is a no-op
    // (an Err) once the pool exists.
    if let Some(n) = config.file.n_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .context("Failed to configure the global thread pool")?;
        log::info!("Thread pool: {n} threads");
    }

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
    let hills = run_hills_streaming(&args.input, &config.hills, &config.file)
        .with_context(|| format!("Failed to detect hills from {}", args.input.display()))?;
    log::info!(
        "[timing] hill detection: {:.2?} ({} hills)",
        t.elapsed(),
        hills.len()
    );
    log_mem(&format!("after hill detection ({} hills)", hills.len()));

    let hills_ext = match config.output.format {
        OutputFormat::Tsv => "tsv",
        OutputFormat::Parquet => "parquet",
    };
    let hills_path = out_dir.join(format!("hills.{hills_ext}"));
    let t = Instant::now();
    match config.output.format {
        OutputFormat::Tsv => write_hills_tsv(&hills, &hills_path),
        OutputFormat::Parquet => write_hills_parquet(&hills, &hills_path),
    }
    .with_context(|| format!("Failed to write hills to {}", hills_path.display()))?;
    log::info!("[timing] write hills.{hills_ext}: {:.2?}", t.elapsed());

    // Optional Stage 1b: MS2 hill detection per isolation window (DIA channels).
    // Independent of MS1 feature finding — we just emit a separate hills_ms2 file.
    if config.file.ms2_hills_enabled {
        log::info!("Detecting MS2 hills (streaming, per isolation window)...");
        let t = Instant::now();
        let ms2_hills = run_ms2_hills_streaming(&args.input, &config.ms2_hills(), &config.file)
            .with_context(|| format!("Failed to detect MS2 hills from {}", args.input.display()))?;
        log::info!(
            "[timing] MS2 hill detection: {:.2?} ({} hills)",
            t.elapsed(),
            ms2_hills.len()
        );
        log_mem(&format!(
            "after MS2 hill detection ({} hills)",
            ms2_hills.len()
        ));

        if !ms2_hills.is_empty() {
            let ms2_path = out_dir.join(format!("hills_ms2.{hills_ext}"));
            let t = Instant::now();
            match config.output.format {
                OutputFormat::Tsv => write_ms2_hills_tsv(&ms2_hills, &ms2_path),
                OutputFormat::Parquet => write_ms2_hills_parquet(&ms2_hills, &ms2_path),
            }
            .with_context(|| format!("Failed to write MS2 hills to {}", ms2_path.display()))?;
            log::info!("[timing] write hills_ms2.{hills_ext}: {:.2?}", t.elapsed());
        } else {
            log::info!("No MS2 hills produced (input may have no MS2 spectra).");
        }
    }

    // Stage 2: Feature detection — hills borrowed here, features clone Hill structs
    // (intensity_profile is Arc so no data duplication)
    log::info!("Detecting features...");
    let t = Instant::now();
    let features =
        run_features(&hills, &config.features, &config.file).context("Feature detection failed")?;
    log::info!(
        "[timing] feature detection: {:.2?} ({} features)",
        t.elapsed(),
        features.len()
    );

    // Build hills summary before we release the hills vec
    let hills_report = build_hills_report(&hills);

    // Drop hills — profile data stays alive via Arc refs inside features
    drop(hills);
    log_mem(&format!(
        "after feature detection ({} features)",
        features.len()
    ));

    // Stage 3: Scoring (optional)
    let scored = if args.no_scoring {
        features
            .iter()
            .map(|f| koth_ff::models::ScoredFeature {
                feature: f.clone(),
                neutron_offset: 0,
                isotope_score: 0.0,
                cosine_score: f.cosine_score,
                combined_score: 0.0,
                theoretical_pattern: Vec::new(),
            })
            .collect::<Vec<_>>()
    } else {
        log::info!("Scoring features...");
        let t = Instant::now();
        let s = run_scoring(
            &features,
            &config.scoring,
            &config.features,
            config.file.polarity,
        );
        log::info!("[timing] scoring: {:.2?}", t.elapsed());
        s
    };

    let features_ext = match config.output.format {
        OutputFormat::Tsv => "tsv",
        OutputFormat::Parquet => "parquet",
    };
    let features_path = out_dir.join(format!("features.{features_ext}"));
    let t = Instant::now();
    match config.output.format {
        OutputFormat::Tsv => write_features_tsv(&scored, &features_path, config.file.polarity),
        OutputFormat::Parquet => {
            write_features_parquet(&scored, &features_path, config.file.polarity)
        }
    }
    .with_context(|| format!("Failed to write features to {}", features_path.display()))?;
    log::info!(
        "[timing] write features.{features_ext}: {:.2?}",
        t.elapsed()
    );

    // Report
    let features_report = build_features_report(&scored);
    let report = RunReport {
        hills: hills_report,
        features: features_report,
    };
    let report_path = out_dir.join("report.json");
    write_report(&report, &report_path).context("Failed to write report")?;

    log::info!("[timing] total: {:.2?}", total_start.elapsed());

    let config_path = out_dir.join("config.toml");
    let config_str = config
        .to_toml_string()
        .context("Failed to serialize config")?;
    std::fs::write(&config_path, config_str)
        .with_context(|| format!("Failed to write config to {}", config_path.display()))?;

    log::info!("Done. Results in {}", out_dir.display());
    Ok(())
}
