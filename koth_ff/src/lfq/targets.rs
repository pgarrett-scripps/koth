//! Identification-guided extraction. Search confidence selects donor targets;
//! it is never used as extraction evidence or as a transfer q-value.
use std::{
    collections::{BTreeMap, HashSet},
    path::Path,
};

use anyhow::{bail, ensure, Context, Result};
use serde::Serialize;

use super::{consensus::ConsensusFeature, IntensityMatrix, LfqConfig};
use crate::{
    alignment::{AlignmentResult, RunInput},
    models::PROTON_MASS,
    output::long_matrix::{sha256, source_rows},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TargetFormat {
    Generic,
    Sage,
}

#[derive(Debug, Clone, Serialize)]
pub struct Identification {
    pub modified_peptide: String,
    pub run_idx: usize,
    pub charge: u8,
    pub neutral_mass: f64,
    pub rt_minutes: f64,
    pub im: f64,
    pub id_qvalue: f64,
    pub source_row: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct ImportSummary {
    pub path: String,
    pub sha256: String,
    pub format: TargetFormat,
    pub max_id_qvalue: f64,
    pub ignore_target_im: bool,
    pub rows_read: usize,
    pub rows_filtered: usize,
    pub rows_retained: usize,
}

/// Normalize only known raw-data suffixes, preserving dots inside sample names.
fn run_key(name: &str) -> String {
    let name = name.rsplit(['/', '\\']).next().unwrap_or(name);
    let mut key = name.to_owned();
    for suffix in [".gz", ".mzml", ".raw", ".d"] {
        if key.to_ascii_lowercase().ends_with(suffix) {
            key.truncate(key.len() - suffix.len());
        }
    }
    key
}

pub fn read_identifications(
    path: &Path,
    format: TargetFormat,
    runs: &[RunInput],
    max_q: f64,
    ignore_target_im: bool,
) -> Result<(Vec<Identification>, ImportSummary)> {
    ensure!(
        max_q.is_finite() && (0.0..=1.0).contains(&max_q),
        "max ID q-value must be in [0, 1]"
    );
    let mut run_ids = BTreeMap::new();
    for (i, run) in runs.iter().enumerate() {
        ensure!(
            run_ids.insert(run_key(&run.name), i).is_none(),
            "Ambiguous normalized run name: {}",
            run.name
        );
    }
    let mut summary = ImportSummary {
        path: path.display().to_string(),
        sha256: sha256(path)?,
        format,
        max_id_qvalue: max_q,
        ignore_target_im,
        rows_read: 0,
        rows_filtered: 0,
        rows_retained: 0,
    };
    let mut records = Vec::new();
    source_rows(path, |header, row| {
        summary.rows_read += 1;
        let parse = || -> Result<Option<Identification>> {
            let unique: HashSet<_> = header.iter().collect();
            ensure!(unique.len() == header.len(), "duplicate column names");
            let get = |name: &str| -> Result<&str> {
                let i = header
                    .iter()
                    .position(|h| h == name)
                    .with_context(|| format!("missing column '{name}'"))?;
                row.get(i)
                    .with_context(|| format!("missing value for '{name}'"))
            };
            let number = |name: &str| -> Result<f64> {
                let v: f64 = get(name)?
                    .parse()
                    .with_context(|| format!("invalid {name}"))?;
                ensure!(v.is_finite(), "{name} must be finite");
                Ok(v)
            };
            let qvalue = |name: &str| -> Result<f64> {
                let v = number(name)?;
                ensure!((0.0..=1.0).contains(&v), "{name} must be in [0, 1]");
                Ok(v)
            };
            let (peptide_col, run_col, mass_col, rt_col, q) = match format {
                TargetFormat::Generic => (
                    "modified_peptide",
                    "run",
                    "neutral_mass",
                    "rt_minutes",
                    qvalue("id_qvalue")?,
                ),
                TargetFormat::Sage => {
                    let decoy = if header.iter().any(|h| h == "is_decoy") {
                        match get("is_decoy")?.to_ascii_lowercase().as_str() {
                            "true" | "1" => true,
                            "false" | "0" => false,
                            _ => bail!("is_decoy must be true/false or 1/0"),
                        }
                    } else {
                        match get("label")? {
                            "1" => false,
                            "-1" => true,
                            _ => bail!("label must be 1 or -1"),
                        }
                    };
                    let rank: u32 = get("rank")?.parse().context("invalid rank")?;
                    if decoy || rank != 1 {
                        return Ok(None);
                    }
                    (
                        "peptide",
                        "filename",
                        "calcmass",
                        "rt",
                        qvalue("peptide_q")?.max(qvalue("spectrum_q")?),
                    )
                }
            };
            if q > max_q {
                return Ok(None);
            }
            let modified_peptide = get(peptide_col)?.trim().to_owned();
            ensure!(!modified_peptide.is_empty(), "empty modified peptide");
            let key = run_key(get(run_col)?);
            let run_idx = *run_ids
                .get(&key)
                .with_context(|| format!("unknown run '{}'", get(run_col).unwrap_or("")))?;
            let charge: u8 = get("charge")?.parse().context("invalid charge")?;
            ensure!(charge > 0, "charge must be positive");
            let neutral_mass = number(mass_col)?;
            let rt_minutes = number(rt_col)?;
            ensure!(
                neutral_mass > 0.0 && rt_minutes >= 0.0,
                "mass must be positive and RT nonnegative"
            );
            let im_col = if format == TargetFormat::Sage {
                "ion_mobility"
            } else {
                "im"
            };
            let im = if !ignore_target_im
                && header.iter().any(|h| h == im_col)
                && !get(im_col)?.is_empty()
            {
                number(im_col)?
            } else {
                0.0
            };
            ensure!(im >= 0.0, "IM must be nonnegative native 1/K0, or empty");
            Ok(Some(Identification {
                modified_peptide,
                run_idx,
                charge,
                neutral_mass,
                rt_minutes,
                im,
                id_qvalue: q,
                source_row: summary.rows_read,
            }))
        };
        match parse()
            .with_context(|| format!("{}: data row {}", path.display(), summary.rows_read))?
        {
            Some(id) => records.push(id),
            None => summary.rows_filtered += 1,
        }
        Ok(())
    })?;
    summary.rows_retained = records.len();
    ensure!(
        !records.is_empty(),
        "No identifications passed the input filters"
    );
    Ok((records, summary))
}

#[derive(Debug, Clone, Serialize)]
pub struct SearchTarget {
    pub modified_peptide: String,
    pub charge: u8,
    /// One best accepted PSM per run; source row is a 1-based data-row number.
    pub donors: Vec<Option<Identification>>,
    pub seed_run_idx: usize,
    pub transfer_eligible: bool,
    /// Recipient runs this target must not be transferred into, even when it is
    /// otherwise transfer eligible. Empty unless per-run gating is enabled.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub transfer_blocked_runs: Vec<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rt_rescue: Option<RtRescue>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub inferred_charge: Option<InferredCharge>,
}

#[derive(Debug, Clone, Serialize)]
pub struct InferredCharge {
    /// A real PSM at its original charge, never a fabricated identification.
    pub seed: Identification,
    pub native_anchors: Vec<Option<Identification>>,
    pub mz_eligible: Vec<bool>,
}

#[derive(Debug, Clone, Serialize)]
pub struct RtRescue {
    pub excluded_source_rows: Vec<usize>,
    pub cross_run_rt_ambiguous: bool,
}

/// Bounded, non-chaining RT clusters. A conflicting run is retained only if
/// one cluster has a strict majority of PSMs and at least two observations.
/// Tied/unsupported clusters withhold that run; one target per peptide/charge
/// avoids duplicate quantification of multiple chromatographic components.
fn supported_rt_ids(
    ids: Vec<Identification>,
    runs: &[RunInput],
    fraction: f64,
) -> (Vec<Identification>, Vec<usize>, Vec<bool>) {
    let mut retained = Vec::new();
    let mut excluded = Vec::new();
    // A run that offered identifications but kept none has genuinely ambiguous
    // chromatography for this peptide, so it must not receive a transfer.
    let mut withheld = vec![false; runs.len()];
    for (r, run) in runs.iter().enumerate() {
        let mut rows: Vec<_> = ids.iter().filter(|id| id.run_idx == r).collect();
        rows.sort_by(|a, b| {
            a.rt_minutes
                .total_cmp(&b.rt_minutes)
                .then(a.source_row.cmp(&b.source_row))
        });
        let (lo, hi) = run.rt_range();
        let limit = (hi - lo) * fraction;
        let mut windows = Vec::new();
        for start in 0..rows.len() {
            let end = rows.partition_point(|id| id.rt_minutes <= rows[start].rt_minutes + limit);
            windows.push((end - start, start, end));
        }
        let maximum = windows.iter().map(|w| w.0).max().unwrap_or(0);
        let best: Vec<_> = windows.iter().filter(|w| w.0 == maximum).collect();
        let supported = best.len() == 1
            && (maximum == rows.len() || (maximum >= 2 && maximum * 2 > rows.len()));
        for (i, id) in rows.iter().enumerate() {
            if supported
                && i >= best[0].1
                && i < best[0].2
                && id.rt_minutes >= lo
                && id.rt_minutes <= hi
            {
                retained.push((*id).clone());
            } else {
                excluded.push(id.source_row);
            }
        }
    }
    excluded.sort_unstable();
    for (r, flag) in withheld.iter_mut().enumerate() {
        let offered = ids.iter().any(|id| id.run_idx == r);
        *flag = offered && !retained.iter().any(|id| id.run_idx == r);
    }
    (retained, excluded, withheld)
}

#[derive(Debug, Clone, Serialize)]
pub struct RejectedTarget {
    pub modified_peptide: String,
    pub charge: u8,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct SearchGuidance {
    pub import: ImportSummary,
    pub mbr: bool,
    pub aligned_runs: Vec<bool>,
    pub targets: Vec<SearchTarget>,
    pub rejected_targets: Vec<RejectedTarget>,
}

impl SearchGuidance {
    pub fn is_direct(&self, target: usize, run: usize) -> bool {
        self.targets[target].donors[run].is_some()
    }
    pub fn attempted(&self, target: usize, run: usize) -> bool {
        self.is_direct(target, run)
            || self.native_anchor(target, run).is_some()
            || (self.mbr
                && self.aligned_runs[run]
                && self.targets[target].transfer_eligible
                && !self.targets[target].transfer_blocked_runs.contains(&run)
                && self.targets[target]
                    .inferred_charge
                    .as_ref()
                    .is_none_or(|i| i.mz_eligible[run]))
    }
    pub fn native_anchor(&self, target: usize, run: usize) -> Option<&Identification> {
        self.targets[target].donors[run].as_ref().or_else(|| {
            self.targets[target]
                .inferred_charge
                .as_ref()
                .and_then(|i| i.native_anchors[run].as_ref())
        })
    }
    pub fn is_inferred(&self, target: usize, run: usize) -> bool {
        let t = &self.targets[target];
        !self.is_direct(target, run)
            && t.inferred_charge.as_ref().is_some_and(|i| {
                i.native_anchors[run].is_some() || t.donors.iter().all(Option::is_none)
            })
    }
}

/// Build peptide/charge targets without requiring a detected isotope envelope.
/// Coordinate-only alignment is retained; unreliable alignment blocks transfers.
/// Incompatible RT/IM observations are rejected instead of merging isomers.
pub fn build_targets(
    records: Vec<Identification>,
    import: ImportSummary,
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
    min_anchors: usize,
    mbr: bool,
) -> Result<(Vec<ConsensusFeature>, SearchGuidance)> {
    ensure!(
        !runs.is_empty() && alignment.reference_idx < runs.len(),
        "invalid reference run"
    );
    let aligned_runs: Vec<_> = runs
        .iter()
        .enumerate()
        .map(|(r, run)| {
            r == alignment.reference_idx
                || alignment.alignments.get(&run.name).is_some_and(|a| {
                    a.rt_active.iter().filter(|&&v| v).count() >= min_anchors.max(2)
                })
        })
        .collect();
    let mut grouped: BTreeMap<(String, u8), Vec<Identification>> = BTreeMap::new();
    for id in records {
        grouped
            .entry((id.modified_peptide.clone(), id.charge))
            .or_default()
            .push(id);
    }
    let mut guidance = SearchGuidance {
        import,
        mbr,
        aligned_runs,
        targets: Vec::new(),
        rejected_targets: Vec::new(),
    };
    let mut consensus = Vec::new();
    let range = runs[alignment.reference_idx].rt_range();
    let rt_limit = (range.1 - range.0) * config.consensus.rt_window_pct;
    for ((peptide, charge), mut ids) in grouped {
        // Never use RT rescue to hide incompatible mass or mobility evidence.
        let mass_min = ids
            .iter()
            .map(|id| id.neutral_mass)
            .fold(f64::INFINITY, f64::min);
        let mass_max = ids
            .iter()
            .map(|id| id.neutral_mass)
            .fold(f64::NEG_INFINITY, f64::max);
        let mass_conflict = (mass_max - mass_min) / mass_min * 1e6 > config.consensus.mz_ppm;
        let mut excluded_source_rows = Vec::new();
        let mut withheld_runs = vec![false; runs.len()];
        if config.search_rt_rescue && !mass_conflict && ids.iter().all(|id| id.im == 0.0) {
            (ids, excluded_source_rows, withheld_runs) =
                supported_rt_ids(ids, runs, config.consensus.rt_window_pct);
            if ids.is_empty() {
                guidance.rejected_targets.push(RejectedTarget {
                    modified_peptide: peptide,
                    charge,
                    reason: "no_supported_native_rt_cluster".into(),
                });
                continue;
            }
        }
        ids.sort_by(|a, b| {
            a.id_qvalue
                .total_cmp(&b.id_qvalue)
                .then_with(|| runs[a.run_idx].name.cmp(&runs[b.run_idx].name))
                .then_with(|| a.rt_minutes.total_cmp(&b.rt_minutes))
                .then_with(|| a.im.total_cmp(&b.im))
                .then_with(|| a.neutral_mass.total_cmp(&b.neutral_mass))
                .then(a.source_row.cmp(&b.source_row))
        });
        let transferable = |id: &Identification| {
            let (lo, hi) = runs[id.run_idx].rt_range();
            guidance.aligned_runs[id.run_idx] && id.rt_minutes >= lo && id.rt_minutes <= hi
        };
        let seed = ids.iter().find(|id| transferable(id)).unwrap_or(&ids[0]);
        let project = |id: &Identification| {
            if id.run_idx == alignment.reference_idx {
                (id.rt_minutes, id.im)
            } else {
                let a = &alignment.alignments[&runs[id.run_idx].name];
                (
                    a.warp_rt(id.rt_minutes),
                    if id.im > 0.0 {
                        a.correct_im(id.im, id.rt_minutes)
                    } else {
                        0.0
                    },
                )
            }
        };
        let (ref_rt, ref_im) = project(seed);
        let mut reason = None;
        let mut rt_values = Vec::new();
        let mut im_values = Vec::new();
        for id in &ids {
            if (id.neutral_mass - seed.neutral_mass).abs() / seed.neutral_mass * 1e6
                > config.consensus.mz_ppm
            {
                reason = Some("inconsistent_mass");
            }
            // Within-run ambiguity remains observable even when alignment failed.
            if id.run_idx == seed.run_idx || guidance.aligned_runs[id.run_idx] {
                let (rt, im) = project(id);
                rt_values.push(rt);
                if im > 0.0 {
                    im_values.push(im);
                }
            }
        }
        let span = |v: &[f64]| {
            v.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                - v.iter().copied().fold(f64::INFINITY, f64::min)
        };
        if span(&rt_values) > rt_limit {
            reason = Some("ambiguous_retention_time");
        }
        if span(&im_values) > config.consensus.im_tolerance {
            reason = Some("ambiguous_ion_mobility");
        }
        // Check repeated PSMs in every donor, including unaligned runs.
        for r in 0..runs.len() {
            let native_rt: Vec<_> = ids
                .iter()
                .filter(|id| id.run_idx == r)
                .map(|id| id.rt_minutes)
                .collect();
            let native_im: Vec<_> = ids
                .iter()
                .filter(|id| id.run_idx == r && id.im > 0.0)
                .map(|id| id.im)
                .collect();
            let range = runs[r].rt_range();
            if span(&native_rt) > (range.1 - range.0) * config.consensus.rt_window_pct {
                reason = Some("ambiguous_retention_time");
            }
            if span(&native_im) > config.consensus.im_tolerance {
                reason = Some("ambiguous_ion_mobility");
            }
        }
        let rt_rescued = config.search_rt_rescue
            && !mass_conflict
            && ids.iter().all(|id| id.im == 0.0)
            && reason == Some("ambiguous_retention_time");
        if let Some(reason) = reason.filter(|_| !rt_rescued) {
            guidance.rejected_targets.push(RejectedTarget {
                modified_peptide: peptide,
                charge,
                reason: reason.into(),
            });
            continue;
        }
        let mut donors = vec![None; runs.len()];
        for id in &ids {
            if donors[id.run_idx].is_none() {
                donors[id.run_idx] = Some(id.clone());
            }
        }
        consensus.push(ConsensusFeature {
            ref_mz: seed.neutral_mass / f64::from(charge) + PROTON_MASS,
            ref_rt,
            ref_im,
            charge,
            theoretical_pattern: config
                .isotope_model
                .model()
                .distribution(seed.neutral_mass)
                .to_vec(),
            neutral_mass: seed.neutral_mass,
            seed_combined_score: f64::NAN,
            seed_run_idx: seed.run_idx,
            n_contributing_runs: 0,
            per_run_feature: vec![None; runs.len()],
            members: Vec::new(),
            seed_feature_idx: u32::MAX,
            group_score: f64::NAN,
            group_qvalue: f64::NAN,
        });
        // Cross-run disagreement invalidates the consensus coordinate a transfer
        // would be centred on, so it still withdraws the target everywhere. A
        // pruned row in one run does not, once eligibility is decided per run.
        let per_run = config.search_rt_rescue && config.search_rt_rescue_per_run_transfers;
        let transfer_eligible = transferable(seed)
            && !rt_rescued
            && (per_run || excluded_source_rows.is_empty());
        let transfer_blocked_runs: Vec<usize> = if per_run {
            withheld_runs
                .iter()
                .enumerate()
                .filter(|(_, blocked)| **blocked)
                .map(|(r, _)| r)
                .collect()
        } else {
            Vec::new()
        };
        guidance.targets.push(SearchTarget {
            modified_peptide: peptide,
            charge,
            donors,
            seed_run_idx: seed.run_idx,
            transfer_eligible,
            transfer_blocked_runs,
            rt_rescue: config.search_rt_rescue.then_some(RtRescue {
                excluded_source_rows,
                cross_run_rt_ambiguous: rt_rescued,
            }),
            inferred_charge: None,
        });
    }
    if config.search_expand_charges {
        expand_charges(&mut consensus, &mut guidance, runs, config, rt_limit);
    }
    // Empty candidate sets still produce an auditable rejection report.
    Ok((consensus, guidance))
}

/// Explore 2+ through 4+ only, bounded by each run's observed feature/hill m/z
/// range. Conflicting peptide RT/mass or IM evidence disables expansion.
fn expand_charges(
    consensus: &mut Vec<ConsensusFeature>,
    g: &mut SearchGuidance,
    runs: &[RunInput],
    config: &LfqConfig,
    rt_limit: f64,
) {
    let original = g.targets.clone();
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, t) in original.iter().enumerate() {
        groups
            .entry(t.modified_peptide.clone())
            .or_default()
            .push(i);
    }
    let bounds: Vec<_> = runs
        .iter()
        .map(|r| {
            let mz: Vec<_> = r
                .features
                .iter()
                .map(|f| f.feature.monoisotopic_mz())
                .chain(r.hills.iter().map(|h| h.mz))
                .filter(|v| v.is_finite())
                .collect();
            (
                mz.iter().copied().fold(f64::INFINITY, f64::min),
                mz.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            )
        })
        .collect();
    for (peptide, indices) in groups {
        // Do not silently recreate a rejected charge/isomer under expansion.
        if g.rejected_targets
            .iter()
            .any(|t| t.modified_peptide == peptide)
        {
            continue;
        }
        let seed_idx = indices[0];
        let source = consensus[seed_idx].clone();
        let rt_min = indices
            .iter()
            .map(|&i| consensus[i].ref_rt)
            .fold(f64::INFINITY, f64::min);
        let rt_max = indices
            .iter()
            .map(|&i| consensus[i].ref_rt)
            .fold(f64::NEG_INFINITY, f64::max);
        if rt_max - rt_min > rt_limit
            || indices.iter().any(|&i| {
                !original[i].transfer_eligible
                    || consensus[i].ref_im != 0.0
                    || (consensus[i].neutral_mass - source.neutral_mass).abs() / source.neutral_mass
                        * 1e6
                        > config.consensus.mz_ppm
            })
        {
            continue;
        }
        let Some(seed) = original[seed_idx].donors[source.seed_run_idx].clone() else {
            continue;
        };
        for charge in 2..=4 {
            let mz = source.neutral_mass / f64::from(charge) + PROTON_MASS;
            let anchors: Vec<_> = runs
                .iter()
                .enumerate()
                .map(|(r, run)| {
                    if mz < bounds[r].0 || mz > bounds[r].1 {
                        return None;
                    }
                    let ids: Vec<_> = indices
                        .iter()
                        .filter_map(|&i| original[i].donors[r].as_ref())
                        .collect();
                    let lo = ids
                        .iter()
                        .map(|id| id.rt_minutes)
                        .fold(f64::INFINITY, f64::min);
                    let hi = ids
                        .iter()
                        .map(|id| id.rt_minutes)
                        .fold(f64::NEG_INFINITY, f64::max);
                    let range = run.rt_range();
                    if hi - lo > (range.1 - range.0) * config.consensus.rt_window_pct
                        || ids.iter().any(|id| id.im != 0.0)
                    {
                        return None;
                    }
                    ids.into_iter()
                        .min_by(|a, b| {
                            a.id_qvalue
                                .total_cmp(&b.id_qvalue)
                                .then(a.charge.cmp(&b.charge))
                                .then(a.source_row.cmp(&b.source_row))
                        })
                        .cloned()
                })
                .collect();
            if anchors.iter().all(Option::is_none) {
                continue;
            }
            let provenance = InferredCharge {
                seed: seed.clone(),
                native_anchors: anchors,
                mz_eligible: bounds
                    .iter()
                    .map(|&(lo, hi)| mz >= lo && mz <= hi)
                    .collect(),
            };
            if let Some(&i) = indices.iter().find(|&&i| original[i].charge == charge) {
                g.targets[i].inferred_charge = Some(provenance);
            } else {
                let mut cf = source.clone();
                cf.charge = charge;
                cf.ref_mz = mz;
                consensus.push(cf);
                g.targets.push(SearchTarget {
                    modified_peptide: peptide.clone(),
                    charge,
                    donors: vec![None; runs.len()],
                    seed_run_idx: source.seed_run_idx,
                    transfer_eligible: true,
                    transfer_blocked_runs: Vec::new(),
                    rt_rescue: None,
                    inferred_charge: Some(provenance),
                });
            }
        }
    }
}

fn finite(v: f64) -> String {
    if v.is_finite() {
        v.to_string()
    } else {
        String::new()
    }
}

/// Write peptide identity, donor provenance, and distinct ID/extraction confidence.
/// Raw intensities remain available; the peptide matrix applies an explicit cell gate.
pub fn write_search_outputs(
    matrix: &IntensityMatrix,
    out: &Path,
    max_cell_q: f64,
    run_tdc: bool,
) -> Result<()> {
    let g = matrix
        .search_guidance
        .as_ref()
        .context("missing search guidance")?;
    let mut w = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(out.join("peptide_quant.tsv"))?;
    w.write_record([
        "target_id",
        "modified_peptide",
        "charge",
        "run",
        "evidence",
        "status",
        "intensity",
        "accepted_intensity",
        "extraction_qvalue",
        "direct_id_qvalue",
        "seed_run",
        "seed_id_qvalue",
        "direct_source_row",
        "seed_source_row",
    ])?;
    let mut wide = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(out.join("peptide_intensity_matrix.tsv"))?;
    let mut header = vec![
        "target_id".into(),
        "modified_peptide".into(),
        "charge".into(),
    ];
    header.extend(matrix.run_names.clone());
    wide.write_record(header)?;
    for (i, t) in g.targets.iter().enumerate() {
        let seed = t.donors[t.seed_run_idx]
            .as_ref()
            .or_else(|| t.inferred_charge.as_ref().map(|i| &i.seed))
            .context("target missing seed ID")?;
        let mut row = vec![
            i.to_string(),
            t.modified_peptide.clone(),
            t.charge.to_string(),
        ];
        for r in 0..matrix.n_runs {
            let direct = t.donors[r].as_ref();
            let attempted = g.attempted(i, r);
            let intensity = matrix.intensity(i, r);
            let q = matrix.q_value(i, r);
            let accepted =
                attempted && intensity > 0.0 && run_tdc && q.is_finite() && q <= max_cell_q;
            let value = if accepted { intensity } else { 0.0 };
            row.push(value.to_string());
            let status = if !attempted {
                "not_attempted"
            } else if intensity <= 0.0 {
                "no_signal"
            } else if !run_tdc {
                "unscored"
            } else if accepted {
                "accepted"
            } else {
                "rejected"
            };
            w.write_record(vec![
                i.to_string(),
                t.modified_peptide.clone(),
                t.charge.to_string(),
                matrix.run_names[r].clone(),
                if direct.is_some() {
                    "direct_ms2"
                } else if g.is_inferred(i, r) {
                    if g.native_anchor(i, r).is_some() {
                        "inferred_charge"
                    } else {
                        "mbr_inferred_charge"
                    }
                } else {
                    "mbr"
                }
                .into(),
                status.into(),
                intensity.to_string(),
                value.to_string(),
                if attempted && run_tdc {
                    finite(q)
                } else {
                    String::new()
                },
                direct.map(|d| finite(d.id_qvalue)).unwrap_or_default(),
                matrix.run_names[t.seed_run_idx].clone(),
                finite(seed.id_qvalue),
                direct.map(|d| d.source_row.to_string()).unwrap_or_default(),
                seed.source_row.to_string(),
            ])?;
        }
        wide.write_record(row)?;
    }
    w.flush()?;
    wide.flush()?;
    let mut rejected = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(out.join("search_rejections.tsv"))?;
    rejected.write_record(["modified_peptide", "charge", "reason"])?;
    for t in &g.rejected_targets {
        rejected.write_record([
            t.modified_peptide.clone(),
            t.charge.to_string(),
            t.reason.clone(),
        ])?;
    }
    rejected.flush()?;
    let manifest = serde_json::json!({"schema_version":1,"koth_version":crate::VERSION,"search":g,
        "run_names":matrix.run_names,"max_cell_qvalue":max_cell_q,"run_tdc":run_tdc,
        "confidence":"ID q-values are imported; extraction q-values are exploratory and separately ranked for direct and transferred cells. No calibrated peptide-transfer FDR is claimed.",
        "outputs":{"peptide_quant.tsv":sha256(&out.join("peptide_quant.tsv"))?,
            "peptide_intensity_matrix.tsv":sha256(&out.join("peptide_intensity_matrix.tsv"))?,
            "search_rejections.tsv":sha256(&out.join("search_rejections.tsv"))?}});
    let temp = out.join("search_manifest.json.tmp");
    std::fs::write(&temp, serde_json::to_vec_pretty(&manifest)?)?;
    std::fs::rename(temp, out.join("search_manifest.json"))?;
    Ok(())
}

#[cfg(test)]
#[path = "targets_tests.rs"]
mod tests;
