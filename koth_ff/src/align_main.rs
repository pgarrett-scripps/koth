use std::io::Write as IoWrite;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::Context;
use clap::Parser;

use koth_ff::{
    alignment::{align_runs, RunInput},
    config::{AlignConfig, OutputFormat},
    input::{discover_runs, read_features, read_hills},
    lfq::{quantify, IntensityMatrix},
    mem::log_mem,
};

#[derive(Parser, Debug)]
#[command(
    name = "koth_align",
    about = "Multi-run alignment and label-free quantification from koth_ff output",
    long_about = "Discovers all per-run subdirectories inside <batch_dir>, aligns retention \
                  time, m/z ppm, and ion mobility across runs, then performs Sage-style LFQ \
                  by extracting ion chromatogram grids from hill data.",
    version
)]
struct Args {
    /// Directory containing one subdirectory per run.
    /// Each subdirectory must contain hills.{tsv,parquet} and features.{tsv,parquet}.
    batch_dir: PathBuf,

    /// Output directory (default: <batch_dir>/align_output)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Alignment + LFQ configuration TOML (uses defaults if not provided)
    #[arg(short, long)]
    config: Option<PathBuf>,

    /// Log level: error, warn, info, debug, trace
    #[arg(long, default_value = "info")]
    log_level: String,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let log_filter = format!("koth_ff={}", args.log_level);
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(&log_filter))
        .init();

    let config = match &args.config {
        Some(p) => AlignConfig::from_toml(p)
            .with_context(|| format!("Failed to load config from {}", p.display()))?,
        None => AlignConfig::default(),
    };

    let out_dir = args
        .output
        .clone()
        .unwrap_or_else(|| args.batch_dir.join("align_output"));
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("Failed to create output directory {}", out_dir.display()))?;

    log::info!("Batch dir:  {}", args.batch_dir.display());
    log::info!("Output dir: {}", out_dir.display());

    let total_start = Instant::now();

    // ── Discover runs ─────────────────────────────────────────────────────────
    let run_paths =
        discover_runs(&args.batch_dir).context("Failed to discover runs in batch directory")?;

    if run_paths.is_empty() {
        anyhow::bail!(
            "No runs found in '{}'. Each run must be a subdirectory containing \
             hills.{{tsv,parquet}} and features.{{tsv,parquet}}.",
            args.batch_dir.display()
        );
    }

    log::info!("Found {} runs:", run_paths.len());
    for rp in &run_paths {
        log::info!(
            "  {} → hills: {}, features: {}",
            rp.name,
            rp.hills_path.display(),
            rp.features_path.display()
        );
    }

    // ── Load runs ─────────────────────────────────────────────────────────────
    let mut runs: Vec<RunInput> = Vec::with_capacity(run_paths.len());

    for rp in &run_paths {
        let t = Instant::now();

        let hills = read_hills(&rp.hills_path)
            .with_context(|| format!("Failed to read hills from {}", rp.hills_path.display()))?;

        let features = read_features(&rp.features_path).with_context(|| {
            format!("Failed to read features from {}", rp.features_path.display())
        })?;

        log::info!(
            "Loaded '{}': {} hills, {} features [{:.2?}]",
            rp.name,
            hills.len(),
            features.len(),
            t.elapsed()
        );

        runs.push(RunInput {
            name: rp.name.clone(),
            features,
            hills,
            scan_times: Vec::new(), // not stored in output files; grid uses rt interpolation
        });
    }

    log_mem("after loading all runs");

    // ── Alignment ─────────────────────────────────────────────────────────────
    log::info!("Running alignment...");
    let t = Instant::now();
    let alignment = align_runs(&runs, &config.alignment);
    log::info!(
        "[timing] alignment: {:.2?} (reference: '{}')",
        t.elapsed(),
        alignment.reference_name
    );

    for (name, al) in &alignment.alignments {
        log::info!(
            "  {}: {} anchors, mass drift {:.3} + {:.4}×RT ppm",
            name,
            al.n_anchors,
            al.mass_drift.intercept,
            al.mass_drift.slope
        );
    }

    // ── LFQ ───────────────────────────────────────────────────────────────────
    log::info!("Running LFQ quantification...");
    let t = Instant::now();
    let matrix = quantify(&runs, &alignment, &config.lfq);
    log::info!(
        "[timing] LFQ: {:.2?} ({} features × {} runs)",
        t.elapsed(),
        matrix.n_features,
        matrix.n_runs
    );

    // ── Write outputs ─────────────────────────────────────────────────────────
    match config.output.format {
        OutputFormat::Tsv => {
            write_consensus_tsv(&matrix, &out_dir, config.output.max_qvalue)?;
            write_matrix_tsv(&matrix, &out_dir)?;
            if config.lfq.run_tdc {
                write_qvalue_matrix_tsv(&matrix, &out_dir)?;
            }
        }
        OutputFormat::Parquet => {
            write_consensus_tsv(&matrix, &out_dir, config.output.max_qvalue)?;
            write_matrix_tsv(&matrix, &out_dir)?;
            if config.lfq.run_tdc {
                write_qvalue_matrix_tsv(&matrix, &out_dir)?;
            }
            log::warn!(
                "Parquet output for the intensity matrix is not yet implemented; \
                 wrote TSV instead."
            );
        }
    }

    // Save the config used
    let cfg_path = out_dir.join("align_config.toml");
    let cfg_str = config
        .to_toml_string()
        .context("Failed to serialize config")?;
    std::fs::write(&cfg_path, cfg_str)
        .with_context(|| format!("Failed to write config to {}", cfg_path.display()))?;

    log::info!("[timing] total: {:.2?}", total_start.elapsed());
    log::info!("Done. Results in {}", out_dir.display());
    Ok(())
}

// ── Output writers ────────────────────────────────────────────────────────────

/// Write `consensus_features.tsv` — one row per consensus feature.
fn write_consensus_tsv(
    matrix: &IntensityMatrix,
    out_dir: &PathBuf,
    max_qvalue: f64,
) -> anyhow::Result<()> {
    let path = out_dir.join("consensus_features.tsv");
    let mut f = std::fs::File::create(&path)?;

    writeln!(
        f,
        "massCalib\tmz\tcharge\trtApex\tim\tscore\tseed_run\tn_contributing_runs\tn_runs_detected"
    )?;

    for feat in 0..matrix.n_features {
        let n_detected = (0..matrix.n_runs)
            .filter(|&run| {
                matrix.intensity(feat, run) > 0.0
                    && matrix.q_value(feat, run) <= max_qvalue
            })
            .count();

        writeln!(
            f,
            "{:.6}\t{:.6}\t{}\t{:.6}\t{}\t{:.6}\t{}\t{}\t{}",
            matrix.feature_mass[feat],
            matrix.feature_mz[feat],
            matrix.feature_charge[feat],
            matrix.feature_rt[feat],
            if matrix.feature_im[feat] != 0.0 {
                format!("{:.6}", matrix.feature_im[feat])
            } else {
                String::new()
            },
            matrix.feature_score[feat],
            matrix.feature_seed_run[feat],
            matrix.feature_n_contributing_runs[feat],
            n_detected,
        )?;
    }

    log::info!("Wrote consensus features to {}", path.display());
    Ok(())
}

/// Write `intensity_matrix.tsv` — features × runs, all intensities unfiltered.
/// Use `qvalue_matrix.tsv` to apply an FDR threshold downstream.
fn write_matrix_tsv(matrix: &IntensityMatrix, out_dir: &PathBuf) -> anyhow::Result<()> {
    let path = out_dir.join("intensity_matrix.tsv");
    let mut f = std::fs::File::create(&path)?;

    write!(f, "massCalib\tmz\tcharge\trtApex\tim\tscore\tseed_run\tn_contributing_runs")?;
    for name in &matrix.run_names {
        write!(f, "\t{}", name)?;
    }
    writeln!(f)?;

    for feat in 0..matrix.n_features {
        write!(
            f,
            "{:.6}\t{:.6}\t{}\t{:.6}\t{}\t{:.6}\t{}\t{}",
            matrix.feature_mass[feat],
            matrix.feature_mz[feat],
            matrix.feature_charge[feat],
            matrix.feature_rt[feat],
            if matrix.feature_im[feat] != 0.0 {
                format!("{:.6}", matrix.feature_im[feat])
            } else {
                String::new()
            },
            matrix.feature_score[feat],
            matrix.feature_seed_run[feat],
            matrix.feature_n_contributing_runs[feat],
        )?;

        for run in 0..matrix.n_runs {
            let v = matrix.intensity(feat, run);
            if v > 0.0 {
                write!(f, "\t{:.5e}", v)?;
            } else {
                write!(f, "\t0")?;
            }
        }
        writeln!(f)?;
    }

    log::info!("Wrote intensity matrix to {}", path.display());
    Ok(())
}

/// Write `qvalue_matrix.tsv` — same layout as intensity_matrix but cells are
/// TDC q-values (1.0 = no signal / TDC not run, 0.0 = perfect score).
/// Use this file to apply an FDR threshold to intensity_matrix.tsv downstream.
fn write_qvalue_matrix_tsv(matrix: &IntensityMatrix, out_dir: &PathBuf) -> anyhow::Result<()> {
    let path = out_dir.join("qvalue_matrix.tsv");
    let mut f = std::fs::File::create(&path)?;

    write!(f, "massCalib\tmz\tcharge\trtApex\tim\tscore\tseed_run\tn_contributing_runs")?;
    for name in &matrix.run_names {
        write!(f, "\t{}", name)?;
    }
    writeln!(f)?;

    for feat in 0..matrix.n_features {
        write!(
            f,
            "{:.6}\t{:.6}\t{}\t{:.6}\t{}\t{:.6}\t{}\t{}",
            matrix.feature_mass[feat],
            matrix.feature_mz[feat],
            matrix.feature_charge[feat],
            matrix.feature_rt[feat],
            if matrix.feature_im[feat] != 0.0 {
                format!("{:.6}", matrix.feature_im[feat])
            } else {
                String::new()
            },
            matrix.feature_score[feat],
            matrix.feature_seed_run[feat],
            matrix.feature_n_contributing_runs[feat],
        )?;

        for run in 0..matrix.n_runs {
            write!(f, "\t{:.4}", matrix.q_value(feat, run))?;
        }
        writeln!(f)?;
    }

    log::info!("Wrote q-value matrix to {}", path.display());
    Ok(())
}
