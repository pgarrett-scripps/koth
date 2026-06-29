use std::io::{BufWriter, Write as IoWrite};
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
    output::{build_align_report, write_align_report, AlignTiming},
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
    let alignment_elapsed = t.elapsed();
    log::info!(
        "[timing] alignment: {:.2?} (reference: '{}')",
        alignment_elapsed,
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
    let lfq_elapsed = t.elapsed();
    log::info!(
        "[timing] LFQ: {:.2?} ({} features × {} runs)",
        lfq_elapsed,
        matrix.n_features,
        matrix.n_runs
    );

    // ── Write outputs ─────────────────────────────────────────────────────────
    if matches!(config.output.format, OutputFormat::Parquet) {
        log::warn!(
            "Parquet output for the intensity matrix is not yet implemented; \
             wrote TSV instead."
        );
    }

    let t = Instant::now();
    write_consensus_tsv(&matrix, &out_dir, config.output.max_qvalue)?;
    log::info!("[timing] write consensus_features.tsv: {:.2?}", t.elapsed());

    let t = Instant::now();
    write_matrix_tsv(&matrix, &out_dir)?;
    log::info!("[timing] write intensity_matrix.tsv: {:.2?}", t.elapsed());

    if config.lfq.run_tdc {
        let t = Instant::now();
        write_qvalue_matrix_tsv(&matrix, &out_dir)?;
        log::info!("[timing] write qvalue_matrix.tsv: {:.2?}", t.elapsed());
        if config.output.export_decoys {
            let t = Instant::now();
            write_decoy_matrix_tsv(&matrix, &out_dir)?;
            log::info!("[timing] write decoy_intensity_matrix.tsv: {:.2?}", t.elapsed());
        }
    }
    if config.output.export_details {
        let t = Instant::now();
        write_lfq_details_tsv(&matrix, &out_dir)?;
        log::info!("[timing] write lfq_details.tsv: {:.2?}", t.elapsed());
    }

    // Save the config used
    let cfg_path = out_dir.join("align_config.toml");
    let cfg_str = config
        .to_toml_string()
        .context("Failed to serialize config")?;
    std::fs::write(&cfg_path, cfg_str)
        .with_context(|| format!("Failed to write config to {}", cfg_path.display()))?;

    // ── Diagnostics report ────────────────────────────────────────────────────
    let total_elapsed = total_start.elapsed();
    let timing = AlignTiming {
        alignment_sec: alignment_elapsed.as_secs_f64(),
        lfq_sec: lfq_elapsed.as_secs_f64(),
        total_sec: total_elapsed.as_secs_f64(),
    };
    let t = Instant::now();
    let report = build_align_report(&runs, &alignment, &matrix, &config, timing);
    log::info!("[timing] build_align_report: {:.2?}", t.elapsed());
    let report_path = out_dir.join("align_report.json");
    let t = Instant::now();
    write_align_report(&report, &report_path).context("Failed to write align report")?;
    log::info!("[timing] write_align_report: {:.2?}", t.elapsed());

    log::info!("[timing] total: {:.2?}", total_elapsed);
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
    let file = std::fs::File::create(&path)?;
    let mut f = BufWriter::with_capacity(1 << 20, file);

    writeln!(
        f,
        "massCalib\tmz\tcharge\trtApex\tim\tcombined_score\tseed_run\tn_contributing_runs\tn_runs_detected"
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
            matrix.feature_combined_score[feat],
            matrix.feature_seed_run[feat],
            matrix.feature_n_contributing_runs[feat],
            n_detected,
        )?;
    }

    f.flush()?;
    log::info!("Wrote consensus features to {}", path.display());
    Ok(())
}

/// Write `intensity_matrix.tsv` — features × runs, all intensities unfiltered.
/// Use `qvalue_matrix.tsv` to apply an FDR threshold downstream.
fn write_matrix_tsv(matrix: &IntensityMatrix, out_dir: &PathBuf) -> anyhow::Result<()> {
    let path = out_dir.join("intensity_matrix.tsv");
    let file = std::fs::File::create(&path)?;
    let mut f = BufWriter::with_capacity(1 << 20, file);

    write!(f, "massCalib\tmz\tcharge\trtApex\tim\tcombined_score\tseed_run\tn_contributing_runs")?;
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
            matrix.feature_combined_score[feat],
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

    f.flush()?;
    log::info!("Wrote intensity matrix to {}", path.display());
    Ok(())
}

/// Write `qvalue_matrix.tsv` — same layout as intensity_matrix but cells are
/// TDC q-values (1.0 = no signal / TDC not run, 0.0 = perfect score).
/// Use this file to apply an FDR threshold to intensity_matrix.tsv downstream.
fn write_qvalue_matrix_tsv(matrix: &IntensityMatrix, out_dir: &PathBuf) -> anyhow::Result<()> {
    let path = out_dir.join("qvalue_matrix.tsv");
    let file = std::fs::File::create(&path)?;
    let mut f = BufWriter::with_capacity(1 << 20, file);

    write!(f, "massCalib\tmz\tcharge\trtApex\tim\tcombined_score\tseed_run\tn_contributing_runs")?;
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
            matrix.feature_combined_score[feat],
            matrix.feature_seed_run[feat],
            matrix.feature_n_contributing_runs[feat],
        )?;

        for run in 0..matrix.n_runs {
            write!(f, "\t{:.4}", matrix.q_value(feat, run))?;
        }
        writeln!(f)?;
    }

    f.flush()?;
    log::info!("Wrote q-value matrix to {}", path.display());
    Ok(())
}

/// Write `decoy_intensity_matrix.tsv` — mirrors `intensity_matrix.tsv` but
/// each cell holds the decoy-pass intensity for that (feature, run).
fn write_decoy_matrix_tsv(matrix: &IntensityMatrix, out_dir: &PathBuf) -> anyhow::Result<()> {
    let path = out_dir.join("decoy_intensity_matrix.tsv");
    let file = std::fs::File::create(&path)?;
    let mut f = BufWriter::with_capacity(1 << 20, file);

    write!(f, "massCalib\tmz\tcharge\trtApex\tim\tcombined_score\tseed_run\tn_contributing_runs")?;
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
            matrix.feature_combined_score[feat],
            matrix.feature_seed_run[feat],
            matrix.feature_n_contributing_runs[feat],
        )?;

        for run in 0..matrix.n_runs {
            let v = matrix.decoy_intensity(feat, run);
            if v > 0.0 {
                write!(f, "\t{:.5e}", v)?;
            } else {
                write!(f, "\t0")?;
            }
        }
        writeln!(f)?;
    }

    f.flush()?;
    log::info!("Wrote decoy intensity matrix to {}", path.display());
    Ok(())
}

/// Write `lfq_details.tsv` — long-format file with one row per
/// (consensus feature, run, target/decoy pass) containing the integration
/// stats, RT diagnostics, and observed monoisotopic m/z + ion mobility.
/// `q_value` is filled for target rows only (decoys get NaN since TDC ranks
/// targets against decoys, not the other way round).
fn write_lfq_details_tsv(matrix: &IntensityMatrix, out_dir: &PathBuf) -> anyhow::Result<()> {
    let path = out_dir.join("lfq_details.tsv");
    let file = std::fs::File::create(&path)?;
    let mut f = BufWriter::with_capacity(1 << 20, file);

    // Column dictionary:
    //   seed_combined_score = seed feature's combined_score (isotope × chromato cosine)
    //   hybrid_score        = cbrt(rt × intensity × spectral_cosine) at apex column
    //   spectral_cosine     = cosine(XIC column, theoretical pattern) at apex
    //                         column (LFQ-stage cosine — NOT the chromatographic
    //                         feature cosine_score)
    writeln!(
        f,
        "feature_idx\tmassCalib\tmz\tcharge\trefRtApex\trefIm\tseed_combined_score\tseed_run\
         \tn_contributing_runs\trun_name\tis_decoy\tis_mbr\tintensity\thybrid_score\tspectral_cosine\
         \tn_isotopes_found\texpected_rt\tapex_rt\trt_diff\tpeak_width_rt\
         \tobserved_mz\tppm_error\tobserved_im\tim_delta\tq_value"
    )?;

    let fmt_f = |v: f64| {
        if v.is_finite() {
            format!("{:.6}", v)
        } else {
            String::from("NaN")
        }
    };

    for entry in &matrix.entries {
        let feat = entry.feature_idx;
        let run = entry.run_idx;
        let intensity_str = if entry.intensity > 0.0 {
            format!("{:.5e}", entry.intensity)
        } else {
            String::from("0")
        };

        // Note: q-values are computed per (feature, run) in target-decoy
        // competition. We surface the target's q-value on the decoy row too
        // for convenience when filtering — same value for both rows in a pair.
        let q_value = if matrix.q_values.is_empty() {
            f64::NAN
        } else {
            matrix.q_value(feat, run)
        };

        // ppm error vs the centre m/z used to build this row's XIC grid.
        // For targets that is the alignment-predicted m/z; for decoys it is
        // the shifted decoy m/z — so target and decoy residuals are on the
        // same "distance from prediction" scale and stay comparable for
        // downstream rescorers.
        let ppm_error = if entry.observed_mz.is_finite() && entry.expected_mz > 0.0 {
            (entry.observed_mz - entry.expected_mz) / entry.expected_mz * 1e6
        } else {
            f64::NAN
        };

        // Per-row IM delta vs the centre IM used to build the grid (same
        // target/decoy symmetry as ppm_error). NaN when either side missing.
        let im_delta = if entry.observed_im.is_finite() && entry.expected_im != 0.0 {
            entry.observed_im - entry.expected_im
        } else {
            f64::NAN
        };

        let rt_diff = if entry.apex_rt.is_finite() && entry.expected_rt.is_finite() {
            entry.apex_rt - entry.expected_rt
        } else {
            f64::NAN
        };

        writeln!(
            f,
            "{feat}\t{mass:.6}\t{mz:.6}\t{charge}\t{rt:.6}\t{im}\t{seed_combined:.6}\t{seed_run}\
             \t{ncont}\t{run_name}\t{decoy}\t{is_mbr}\t{intensity}\t{hybrid:.6}\t{spectral:.6}\
             \t{nslots}\t{exp_rt}\t{apex_rt}\t{rtdiff}\t{pw}\
             \t{obs_mz}\t{ppm}\t{obs_im}\t{im_delta}\t{qv}",
            feat = feat,
            mass = matrix.feature_mass[feat],
            mz = matrix.feature_mz[feat],
            charge = matrix.feature_charge[feat],
            rt = matrix.feature_rt[feat],
            im = if matrix.feature_im[feat] != 0.0 {
                format!("{:.6}", matrix.feature_im[feat])
            } else {
                String::new()
            },
            seed_combined = matrix.feature_combined_score[feat],
            seed_run = matrix.feature_seed_run[feat],
            ncont = matrix.feature_n_contributing_runs[feat],
            run_name = matrix.run_names[run],
            decoy = entry.is_decoy,
            is_mbr = entry.is_mbr,
            intensity = intensity_str,
            hybrid = entry.hybrid_score,
            spectral = entry.spectral_cosine,
            nslots = entry.n_isotopes_found,
            exp_rt = fmt_f(entry.expected_rt),
            apex_rt = fmt_f(entry.apex_rt),
            rtdiff = fmt_f(rt_diff),
            pw = fmt_f(entry.peak_width_rt),
            obs_mz = fmt_f(entry.observed_mz),
            ppm = fmt_f(ppm_error),
            obs_im = fmt_f(entry.observed_im),
            im_delta = fmt_f(im_delta),
            qv = fmt_f(q_value),
        )?;
    }

    f.flush()?;
    log::info!(
        "Wrote LFQ details ({} rows) to {}",
        matrix.entries.len(),
        path.display()
    );
    Ok(())
}
