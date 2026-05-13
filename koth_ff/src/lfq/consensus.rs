use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::alignment::{AlignmentResult, RunInput};
use crate::lfq::LfqConfig;

/// Grouping quality filters for the multi-run consensus feature list.
/// Tolerances (m/z ppm, RT window, IM) are shared with the LFQ config to avoid
/// redundant settings — see `[lfq]` for those values.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsensusConfig {
    /// Minimum feature score for a feature to be eligible as a group member or seed.
    pub min_member_score: f64,
    /// Minimum number of runs that must have detected a feature in a group for it
    /// to be used for LFQ. 1 = include single-run features (default). 2 = require
    /// detection in at least two runs before attempting LFQ across all runs.
    #[serde(default = "default_min_group_size")]
    pub min_group_size: usize,
    /// Minimum isotope-pattern score of the group's seed feature. Groups whose
    /// best-scoring member falls below this threshold are dropped entirely.
    /// 0.0 = keep all groups (default).
    #[serde(default)]
    pub min_seed_score: f64,
}

fn default_min_group_size() -> usize {
    1
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            min_member_score: 0.0,
            min_group_size: 1,
            min_seed_score: 0.0,
        }
    }
}

/// One row in the multi-run consensus feature list.
#[derive(Debug, Clone)]
pub struct ConsensusFeature {
    /// Monoisotopic m/z in reference-run space.
    pub ref_mz: f64,
    /// Retention time apex in reference-run space (minutes).
    pub ref_rt: f64,
    /// Ion mobility apex in reference-run space (0.0 if not available).
    pub ref_im: f64,
    pub charge: u8,
    /// Theoretical averagine isotope pattern (normalised), from the seed's ScoredFeature.
    pub theoretical_pattern: Vec<f64>,
    /// Neutral monoisotopic mass in reference space.
    pub neutral_mass: f64,
    /// Isotope-pattern score of the seed feature.
    pub seed_score: f64,
    /// Index into the `runs` slice that supplied the seed.
    pub seed_run_idx: usize,
    /// How many runs contributed at least one feature to this group.
    pub n_contributing_runs: usize,
}

/// Internal: one feature projected into reference-run coordinate space.
struct ProjectedFeature {
    run_idx: usize,
    ref_mz: f64,
    ref_rt: f64,
    ref_im: f64,
    charge: u8,
    neutral_mass: f64,
    score: f64,
    theoretical_pattern: Vec<f64>,
}

fn emit_group(projected: &[ProjectedFeature], group: &[usize]) -> ConsensusFeature {
    // Keep best-scoring feature per run to avoid double-counting split peaks.
    let mut best_per_run: HashMap<usize, usize> = HashMap::new();
    for &i in group {
        let run_idx = projected[i].run_idx;
        let entry = best_per_run.entry(run_idx).or_insert(i);
        if projected[i].score > projected[*entry].score {
            *entry = i;
        }
    }

    let deduped: Vec<usize> = best_per_run.values().copied().collect();
    let seed_idx = *deduped
        .iter()
        .max_by(|&&a, &&b| projected[a].score.partial_cmp(&projected[b].score).unwrap())
        .unwrap();
    let seed = &projected[seed_idx];

    ConsensusFeature {
        ref_mz: seed.ref_mz,
        ref_rt: seed.ref_rt,
        ref_im: seed.ref_im,
        charge: seed.charge,
        theoretical_pattern: seed.theoretical_pattern.clone(),
        neutral_mass: seed.neutral_mass,
        seed_score: seed.score,
        seed_run_idx: seed.run_idx,
        n_contributing_runs: deduped.len(),
    }
}

/// Apply min_group_size and min_seed_score filters before adding a group to the output.
#[inline]
fn push_if_passes(cf: ConsensusFeature, config: &LfqConfig, out: &mut Vec<ConsensusFeature>) {
    if cf.n_contributing_runs >= config.consensus.min_group_size
        && cf.seed_score >= config.consensus.min_seed_score
    {
        out.push(cf);
    }
}

/// Build a multi-run consensus feature list.
///
/// All features from all runs are projected into reference-run coordinate space
/// using the alignment corrections, then grouped by (charge, neutral mass ± ppm,
/// aligned RT ± window, IM ± tolerance).  Each group emits one `ConsensusFeature`
/// seeded by the highest-scoring member from any run.
pub fn build_consensus(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
) -> Vec<ConsensusFeature> {
    const PROTON: f64 = 1.007_276_466_621;

    let mut projected: Vec<ProjectedFeature> = Vec::new();

    for (run_idx, run) in runs.iter().enumerate() {
        let is_reference = run_idx == alignment.reference_idx;

        for feat in &run.features {
            if feat.feature.charge == 0 || feat.score < config.consensus.min_member_score {
                continue;
            }

            let native_mz = feat.monoisotopic_mz();
            let native_rt = feat.feature.rt_apex();
            let native_im = feat.feature.im_apex();

            let (ref_mz, ref_rt, ref_im) = if is_reference {
                (native_mz, native_rt, native_im)
            } else {
                let al = &alignment.alignments[&run.name];
                (
                    al.correct_mz(native_mz, native_rt),
                    al.warp_rt(native_rt),
                    al.correct_im(native_im, native_rt),
                )
            };

            let neutral_mass =
                ref_mz * feat.feature.charge as f64 - feat.feature.charge as f64 * PROTON;

            projected.push(ProjectedFeature {
                run_idx,
                ref_mz,
                ref_rt,
                ref_im,
                charge: feat.feature.charge,
                neutral_mass,
                score: feat.score,
                theoretical_pattern: feat.theoretical_pattern.clone(),
            });
        }
    }

    // Sort: charge ASC, neutral_mass ASC, ref_im ASC, ref_rt ASC.
    // Including IM before RT ensures features at the same mass but very different
    // ion mobilities (different conformers, or noise) are separated in the list
    // before the sweep runs, preventing an IM-mismatched feature from landing
    // between two same-IM features and breaking the group.
    projected.sort_by(|a, b| {
        a.charge
            .cmp(&b.charge)
            .then(a.neutral_mass.partial_cmp(&b.neutral_mass).unwrap())
            .then(a.ref_im.partial_cmp(&b.ref_im).unwrap())
            .then(a.ref_rt.partial_cmp(&b.ref_rt).unwrap())
    });

    let mut consensus: Vec<ConsensusFeature> = Vec::new();

    if projected.is_empty() {
        return consensus;
    }

    // Convert the fractional RT window to absolute minutes using the reference
    // run's gradient span, computed once before the sweep.
    let ref_rt_range = runs[alignment.reference_idx].rt_range();
    let rt_window_abs = config.rt_window_pct * (ref_rt_range.1 - ref_rt_range.0);


    // Greedy single-linkage sweep: compare each feature against the first member
    // (anchor) of the current group.  Comparing against the anchor rather than the
    // last member prevents drift in group coordinates as members accumulate.
    let mut group: Vec<usize> = vec![0];

    for i in 1..projected.len() {
        let anchor = &projected[group[0]];
        let cur = &projected[i];

        let same_charge = cur.charge == anchor.charge;

        // Use 2× the LFQ extraction ppm for consensus grouping: alignment
        // residual drift can scatter the same peptide's projected mass across
        // >10 ppm, which would fragment one real group into many small ones.
        let mass_ok = if anchor.neutral_mass > 0.0 {
            (cur.neutral_mass - anchor.neutral_mass).abs() / anchor.neutral_mass * 1e6
                <= config.mz_ppm * 2.0
        } else {
            (cur.neutral_mass - anchor.neutral_mass).abs() <= 0.02
        };

        let rt_ok = (cur.ref_rt - anchor.ref_rt).abs() <= rt_window_abs * 2.0;

        let im_ok = anchor.ref_im == 0.0
            || cur.ref_im == 0.0
            || (cur.ref_im - anchor.ref_im).abs() <= config.im_tolerance;

        if same_charge && mass_ok && rt_ok && im_ok {
            group.push(i);
        } else {
            push_if_passes(emit_group(&projected, &group), config, &mut consensus);
            group.clear();
            group.push(i);
        }
    }
    push_if_passes(emit_group(&projected, &group), config, &mut consensus);

    // Second pass: merge consensus groups whose seeds are within tolerance.
    //
    // The first sweep can fragment one real peptide into multiple groups when a
    // run contributes features at the same mass but different retention times
    // (e.g. RT=67 and RT=79 for the same peptide in the same run).  Those
    // "wrong-RT" entries interleave in the sorted projected list and act as
    // group separators, closing the current group prematurely.  The resulting
    // sub-groups have seeds that are trivially close to each other.  A second
    // sweep over the seeds merges them.
    consensus = merge_consensus(consensus, config, rt_window_abs);

    let n_multi_run = consensus.iter().filter(|c| c.n_contributing_runs > 1).count();
    log::info!(
        "Consensus: {} groups from {} projected features ({} runs, {} cross-run)",
        consensus.len(),
        projected.len(),
        runs.len(),
        n_multi_run,
    );

    consensus
}

/// Merge consensus groups whose seeds fall within the same mass/IM/RT window.
///
/// Groups are re-sorted by (charge, neutral_mass, ref_im, ref_rt) and swept
/// with the same 2× tolerances used in the projection sweep.  When two
/// adjacent groups are within tolerance, the higher-scoring seed wins and
/// n_contributing_runs is summed (may overcount runs that appeared in both
/// groups, but this is acceptable).
fn merge_consensus(
    mut groups: Vec<ConsensusFeature>,
    config: &LfqConfig,
    rt_window_abs: f64,
) -> Vec<ConsensusFeature> {
    if groups.len() < 2 {
        return groups;
    }

    groups.sort_by(|a, b| {
        a.charge
            .cmp(&b.charge)
            .then(a.neutral_mass.partial_cmp(&b.neutral_mass).unwrap())
            .then(a.ref_im.partial_cmp(&b.ref_im).unwrap())
            .then(a.ref_rt.partial_cmp(&b.ref_rt).unwrap())
    });

    let mut merged: Vec<ConsensusFeature> = Vec::with_capacity(groups.len());
    let mut cur = groups.remove(0);

    for next in groups {
        let same_charge = next.charge == cur.charge;

        let mass_ok = if cur.neutral_mass > 0.0 {
            (next.neutral_mass - cur.neutral_mass).abs() / cur.neutral_mass * 1e6
                <= config.mz_ppm * 2.0
        } else {
            (next.neutral_mass - cur.neutral_mass).abs() <= 0.02
        };

        let rt_ok = (next.ref_rt - cur.ref_rt).abs() <= rt_window_abs * 2.0;

        let im_ok = cur.ref_im == 0.0
            || next.ref_im == 0.0
            || (next.ref_im - cur.ref_im).abs() <= config.im_tolerance;

        if same_charge && mass_ok && rt_ok && im_ok {
            // Absorb next into cur, keeping the better seed.
            cur.n_contributing_runs += next.n_contributing_runs;
            if next.seed_score > cur.seed_score {
                cur.ref_mz = next.ref_mz;
                cur.ref_rt = next.ref_rt;
                cur.ref_im = next.ref_im;
                cur.neutral_mass = next.neutral_mass;
                cur.theoretical_pattern = next.theoretical_pattern;
                cur.seed_score = next.seed_score;
                cur.seed_run_idx = next.seed_run_idx;
            }
        } else {
            merged.push(cur);
            cur = next;
        }
    }
    merged.push(cur);
    merged
}
