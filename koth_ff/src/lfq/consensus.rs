use std::collections::HashMap;
use std::time::Instant;

use serde::{Deserialize, Serialize};

use crate::alignment::{AlignmentResult, RunInput};
use crate::lfq::LfqConfig;

/// Grouping quality filters for the multi-run consensus feature list.
/// Tolerances (m/z ppm, RT window, IM) are shared with the LFQ config to avoid
/// redundant settings — see `[lfq]` for those values.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConsensusConfig {
    /// Pre-grouping filter: a feature's `combined_score` must clear this to
    /// be eligible as a group member or seed. Features below are excluded
    /// from the projection entirely (don't count toward `n_contributing_runs`).
    pub min_member_combined_score: f64,
    /// Minimum number of distinct runs that must have detected a feature in
    /// a group for it to survive. 1 = MBR on (single-run detections kept).
    #[serde(default = "default_min_group_size")]
    pub min_group_size: usize,
    /// Post-grouping filter: drop groups whose seed (best-`combined_score`
    /// member) is below this. 0.0 = keep all groups.
    #[serde(default)]
    pub min_seed_combined_score: f64,
}

fn default_min_group_size() -> usize {
    1
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            min_member_combined_score: 0.0,
            min_group_size: 1,
            min_seed_combined_score: 0.0,
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
    /// `combined_score` (isotope × chromato cosine) of the seed feature.
    pub seed_combined_score: f64,
    /// Index into the `runs` slice that supplied the seed.
    pub seed_run_idx: usize,
    /// How many runs contributed at least one feature to this group.
    pub n_contributing_runs: usize,
    /// For each run index in `runs`, the per-run feature index that joined this
    /// group (highest combined_score among that run's contributing features) or
    /// `None` if the run contributed nothing. Length = n_runs.
    ///
    /// Used by `lfq::quantify` to take the per-run feature's
    /// `total_intensity()` directly for contributing runs, falling back to XIC
    /// re-integration only for runs that did not contribute (MBR transfer).
    pub per_run_feature: Vec<Option<u32>>,
}

/// Internal: one feature projected into reference-run coordinate space.
/// `combined_score` is the seed-selection / filter axis (isotope × chromato cosine).
struct ProjectedFeature {
    run_idx: usize,
    /// Index into `runs[run_idx].features`. Carried through grouping so each
    /// contributing run's feature intensity can be recovered downstream.
    feature_idx: u32,
    ref_mz: f64,
    ref_rt: f64,
    ref_im: f64,
    charge: u8,
    neutral_mass: f64,
    combined_score: f64,
    theoretical_pattern: Vec<f64>,
}

fn emit_group(projected: &[ProjectedFeature], group: &[usize], n_runs: usize) -> ConsensusFeature {
    // Keep best-combined-scoring feature per run to avoid double-counting split peaks.
    let mut best_per_run: HashMap<usize, usize> = HashMap::new();
    for &i in group {
        let run_idx = projected[i].run_idx;
        let entry = best_per_run.entry(run_idx).or_insert(i);
        if projected[i].combined_score > projected[*entry].combined_score {
            *entry = i;
        }
    }

    // Sort the deduped members before selecting the seed: `best_per_run.values()`
    // iterates a HashMap in nondeterministic order, and `max_by` returns the
    // *last* maximum, so on a combined_score tie the seed (and hence the group's
    // m/z / RT, which drive the downstream merge) would flip run-to-run. Sorting
    // by projected index gives a stable, deterministic tie-break.
    let mut deduped: Vec<usize> = best_per_run.values().copied().collect();
    deduped.sort_unstable();
    let seed_idx = *deduped
        .iter()
        .max_by(|&&a, &&b| {
            projected[a]
                .combined_score
                .partial_cmp(&projected[b].combined_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap();
    let seed = &projected[seed_idx];

    // Densify the per-run-feature map so downstream consumers can do O(1) lookup.
    let mut per_run_feature: Vec<Option<u32>> = vec![None; n_runs];
    for (&run_idx, &proj_i) in &best_per_run {
        per_run_feature[run_idx] = Some(projected[proj_i].feature_idx);
    }

    ConsensusFeature {
        ref_mz: seed.ref_mz,
        ref_rt: seed.ref_rt,
        ref_im: seed.ref_im,
        charge: seed.charge,
        theoretical_pattern: seed.theoretical_pattern.clone(),
        neutral_mass: seed.neutral_mass,
        seed_combined_score: seed.combined_score,
        seed_run_idx: seed.run_idx,
        n_contributing_runs: deduped.len(),
        per_run_feature,
    }
}

/// Apply post-grouping filters before emitting a group.
#[inline]
fn push_if_passes(cf: ConsensusFeature, config: &LfqConfig, out: &mut Vec<ConsensusFeature>) {
    if cf.n_contributing_runs >= config.consensus.min_group_size
        && cf.seed_combined_score >= config.consensus.min_seed_combined_score
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

    let t_project = Instant::now();
    let mut projected: Vec<ProjectedFeature> = Vec::new();

    for (run_idx, run) in runs.iter().enumerate() {
        let is_reference = run_idx == alignment.reference_idx;

        for (feat_idx, feat) in run.features.iter().enumerate() {
            if feat.feature.charge == 0
                || feat.combined_score < config.consensus.min_member_combined_score
            {
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
                feature_idx: feat_idx as u32,
                ref_mz,
                ref_rt,
                ref_im,
                charge: feat.feature.charge,
                neutral_mass,
                combined_score: feat.combined_score,
                theoretical_pattern: feat.theoretical_pattern.clone(),
            });
        }
    }

    log::info!(
        "[timing] consensus project ({} features): {:.2?}",
        projected.len(),
        t_project.elapsed()
    );

    // Sort: charge ASC, neutral_mass ASC, ref_im ASC, ref_rt ASC.
    // Including IM before RT ensures features at the same mass but very different
    // ion mobilities (different conformers, or noise) are separated in the list
    // before the sweep runs, preventing an IM-mismatched feature from landing
    // between two same-IM features and breaking the group.
    let t_sort = Instant::now();
    projected.sort_by(|a, b| {
        a.charge
            .cmp(&b.charge)
            .then(a.neutral_mass.partial_cmp(&b.neutral_mass).unwrap())
            .then(a.ref_im.partial_cmp(&b.ref_im).unwrap())
            .then(a.ref_rt.partial_cmp(&b.ref_rt).unwrap())
    });
    log::info!("[timing] consensus sort: {:.2?}", t_sort.elapsed());

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
    let t_sweep = Instant::now();
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
            push_if_passes(
                emit_group(&projected, &group, runs.len()),
                config,
                &mut consensus,
            );
            group.clear();
            group.push(i);
        }
    }
    push_if_passes(
        emit_group(&projected, &group, runs.len()),
        config,
        &mut consensus,
    );
    log::info!(
        "[timing] consensus sweep ({} pre-merge groups): {:.2?}",
        consensus.len(),
        t_sweep.elapsed()
    );

    // Second pass: merge consensus groups whose seeds are within tolerance.
    //
    // The first sweep can fragment one real peptide into multiple groups when a
    // run contributes features at the same mass but different retention times
    // (e.g. RT=67 and RT=79 for the same peptide in the same run).  Those
    // "wrong-RT" entries interleave in the sorted projected list and act as
    // group separators, closing the current group prematurely.  The resulting
    // sub-groups have seeds that are trivially close to each other.  A second
    // sweep over the seeds merges them.
    let t_merge = Instant::now();
    consensus = merge_consensus(consensus, config, rt_window_abs);
    log::info!(
        "[timing] consensus merge_consensus ({} post-merge groups): {:.2?}",
        consensus.len(),
        t_merge.elapsed()
    );

    let n_multi_run = consensus
        .iter()
        .filter(|c| c.n_contributing_runs > 1)
        .count();
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
/// Groups are sorted by (charge, neutral_mass, ref_im, ref_rt) and a sliding
/// mass window is walked over them. Within each window every pair of groups
/// is tested with the 2× tolerances used in the projection sweep, and matching
/// pairs are unioned in a disjoint-set structure. Each resulting cluster
/// collapses to a single feature whose seed is the highest-scoring member;
/// `n_contributing_runs` is summed across members (may overcount when a run
/// appeared in multiple sub-groups, but this matches the previous behaviour).
///
/// The pairwise sweep (vs. the previous adjacent-only sweep) is what makes
/// the merge transitive: it can rejoin sub-groups even when an unrelated
/// group of similar mass interleaves between them in the sort.
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

    let n = groups.len();
    let ppm_tol = config.mz_ppm * 2.0;
    let rt_tol = rt_window_abs * 2.0;
    let im_tol = config.im_tolerance;

    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], mut a: usize) -> usize {
        while parent[a] != a {
            parent[a] = parent[parent[a]];
            a = parent[a];
        }
        a
    }

    // Sliding mass window: for each i, scan forward j>i while same charge and
    // neutral_mass(j) - neutral_mass(i) within the ppm tolerance.
    for i in 0..n {
        let a_charge = groups[i].charge;
        let a_mass = groups[i].neutral_mass;
        let a_rt = groups[i].ref_rt;
        let a_im = groups[i].ref_im;
        let mass_high = if a_mass > 0.0 {
            a_mass * (1.0 + ppm_tol / 1e6)
        } else {
            a_mass + 0.02
        };

        for j in (i + 1)..n {
            let b = &groups[j];
            // charges are sorted ascending; once b.charge > a.charge no
            // further candidate can match.
            if b.charge != a_charge {
                break;
            }
            if b.neutral_mass > mass_high {
                break;
            }

            let rt_ok = (b.ref_rt - a_rt).abs() <= rt_tol;
            if !rt_ok {
                continue;
            }

            let im_ok = a_im == 0.0 || b.ref_im == 0.0 || (b.ref_im - a_im).abs() <= im_tol;
            if !im_ok {
                continue;
            }

            let pi = find(&mut parent, i);
            let pj = find(&mut parent, j);
            if pi != pj {
                parent[pi] = pj;
            }
        }
    }

    // Bucket members by cluster root.
    let mut cluster_members: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let root = find(&mut parent, i);
        cluster_members.entry(root).or_default().push(i);
    }

    // Collapse each cluster: best-scoring seed wins; per_run_feature is the
    // union across cluster members (per-run, prefer the member whose feature
    // has the highest combined_score — approximated by seed score of the
    // group that owned that run-slot).
    let mut merged: Vec<ConsensusFeature> = Vec::with_capacity(cluster_members.len());
    for (_, members) in cluster_members {
        if members.len() == 1 {
            merged.push(groups[members[0]].clone());
            continue;
        }
        let best = *members
            .iter()
            .max_by(|&&a, &&b| {
                groups[a]
                    .seed_combined_score
                    .partial_cmp(&groups[b].seed_combined_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap();
        let mut out = groups[best].clone();
        // Union per_run_feature: for each run-slot that's None in `out`, take
        // the first non-None entry from any other cluster member. Iterating
        // members in score-descending order ensures we prefer higher-quality
        // owner-groups when more than one member fills the same slot.
        let mut order: Vec<usize> = members.clone();
        order.sort_by(|&a, &b| {
            groups[b]
                .seed_combined_score
                .partial_cmp(&groups[a].seed_combined_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for &k in &order {
            if k == best {
                continue;
            }
            let src = &groups[k].per_run_feature;
            for (slot, dst) in out.per_run_feature.iter_mut().enumerate() {
                if dst.is_none() {
                    if let Some(v) = src.get(slot).and_then(|x| *x) {
                        *dst = Some(v);
                    }
                }
            }
        }
        out.n_contributing_runs = out.per_run_feature.iter().filter(|x| x.is_some()).count();
        merged.push(out);
    }

    // Deterministic ordering for downstream consumers.
    merged.sort_by(|a, b| {
        a.charge
            .cmp(&b.charge)
            .then(a.neutral_mass.partial_cmp(&b.neutral_mass).unwrap())
            .then(a.ref_im.partial_cmp(&b.ref_im).unwrap())
            .then(a.ref_rt.partial_cmp(&b.ref_rt).unwrap())
    });
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    fn proj(
        run_idx: usize,
        feature_idx: u32,
        combined_score: f64,
        ref_mz: f64,
    ) -> ProjectedFeature {
        ProjectedFeature {
            run_idx,
            feature_idx,
            ref_mz,
            ref_rt: 10.0,
            ref_im: 0.0,
            charge: 2,
            neutral_mass: 1000.0,
            combined_score,
            theoretical_pattern: vec![0.6, 0.3, 0.1],
        }
    }

    /// A run contributing two features (a split peak) to one group collapses to
    /// its highest-combined-score feature; `n_contributing_runs` counts distinct
    /// runs, not features.
    #[test]
    fn split_peak_dedup_counts_runs_not_features() {
        let projected = vec![
            proj(0, 10, 0.5, 500.00), // run 0, weaker
            proj(0, 11, 0.9, 500.10), // run 0, stronger -> wins for run 0
            proj(1, 20, 0.7, 500.20), // run 1
        ];
        let cf = emit_group(&projected, &[0, 1, 2], 2);
        assert_eq!(cf.n_contributing_runs, 2, "two runs, not three features");
        assert_eq!(cf.seed_run_idx, 0); // global best (0.9) is run 0
        assert!((cf.ref_mz - 500.10).abs() < 1e-9);
        assert_eq!(cf.per_run_feature[0], Some(11)); // run 0's best-scoring feature
        assert_eq!(cf.per_run_feature[1], Some(20));
    }

    /// Seed selection is deterministic across repeated calls despite
    /// `best_per_run` being a per-instance-randomly-seeded HashMap. Guards the
    /// consensus non-determinism fix (deduped members sorted before `max_by`):
    /// on a combined_score tie the HashMap iteration order would otherwise flip
    /// the seed (and hence the group m/z/RT) run-to-run.
    #[test]
    fn seed_selection_is_deterministic_under_score_tie() {
        let projected = vec![proj(0, 1, 0.8, 500.0), proj(1, 2, 0.8, 600.0)];
        let first = emit_group(&projected, &[0, 1], 2);
        for _ in 0..64 {
            let cf = emit_group(&projected, &[0, 1], 2);
            assert_eq!(
                cf.ref_mz, first.ref_mz,
                "seed m/z must be stable across calls"
            );
            assert_eq!(cf.seed_run_idx, first.seed_run_idx);
        }
    }
}
