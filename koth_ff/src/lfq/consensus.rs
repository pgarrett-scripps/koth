use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

use serde::{Deserialize, Serialize};

use crate::alignment::{AlignmentResult, RunInput};
use crate::lfq::LfqConfig;

#[path = "group_confidence.rs"]
mod group_confidence;
pub use group_confidence::{EvidenceAudit, PermutationAudit, VARIANT_NAMES};

type GroupSpans = ((f64, f64), (f64, f64), Option<(f64, f64)>);

/// Grouping quality filters for the multi-run consensus feature list.
/// Grouping limits are independent of LFQ extraction windows.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ConsensusConfig {
    /// Maximum full group mass span in ppm. Default 20 ppm.
    pub mz_ppm: f64,
    /// Maximum full group RT span as a fraction of the reference gradient.
    pub rt_window_pct: f64,
    /// Maximum full group IM span; missing IM never bridges incompatible values.
    pub im_tolerance: f64,
    /// Experimental group-level target/permuted-RT q-value gate. This is an
    /// empirical grouping statistic, not validated identification or cell FDR.
    /// 1.0 retains all replicated candidates for diagnostics. Default 0.05.
    pub max_group_qvalue: f64,
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            mz_ppm: 20.0,
            rt_window_pct: 0.02,
            im_tolerance: 0.05,
            max_group_qvalue: 0.05,
        }
    }
}

impl ConsensusConfig {
    pub fn validate(&self) -> Result<(), crate::error::KothError> {
        for (name, value) in [
            ("mz_ppm", self.mz_ppm),
            ("rt_window_pct", self.rt_window_pct),
            ("im_tolerance", self.im_tolerance),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(crate::error::KothError::ConfigError(format!(
                    "lfq.consensus.{name} must be finite and nonnegative"
                )));
            }
        }
        if !self.max_group_qvalue.is_finite() || !(0.0..=1.0).contains(&self.max_group_qvalue) {
            return Err(crate::error::KothError::ConfigError(
                "lfq.consensus.max_group_qvalue must be between 0 and 1".into(),
            ));
        }
        Ok(())
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
    /// Every original member, including alternatives to the primary per run.
    pub members: Vec<(usize, u32)>,
    /// Original feature index of the seed in its source run. Search-guided
    /// candidates have no feature seed and use u32::MAX; their ID provenance
    /// is carried separately by SearchGuidance.
    pub seed_feature_idx: u32,
    /// Continuous cross-run evidence; the detector score is not a probability.
    pub group_score: f64,
    /// Experimental permuted-RT group q-value, separate from cell confidence.
    pub group_qvalue: f64,
}

/// Internal: one feature projected into reference-run coordinate space.
/// `combined_score` orders candidate assembly; it is no longer an admission gate.
#[derive(Clone)]
struct ProjectedFeature {
    run_idx: usize,
    /// Rank of the run name; independent of input run order.
    run_rank: usize,
    /// Stable tie-break from original scalar measurements.
    measurement_key: u64,
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
        if quality_order(&projected[i], &projected[*entry]).is_lt() {
            *entry = i;
        }
    }

    let seed_idx = *group
        .iter()
        .min_by(|&&a, &&b| quality_order(&projected[a], &projected[b]))
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
        n_contributing_runs: best_per_run.len(),
        per_run_feature,
        members: group
            .iter()
            .map(|&i| (projected[i].run_idx, projected[i].feature_idx))
            .collect(),
        seed_feature_idx: seed.feature_idx,
        group_score: 0.0,
        group_qvalue: 1.0,
    }
}

/// Build a multi-run consensus feature list.
///
/// All features from all runs are projected into reference-run coordinate space
/// using the alignment corrections, then grouped by (charge, neutral mass ± ppm,
/// aligned RT ± window, IM ± tolerance).  Each group emits one `ConsensusFeature`
/// seeded by the highest-scoring member from any run.
fn project_features(runs: &[RunInput], alignment: &AlignmentResult) -> Vec<ProjectedFeature> {
    const PROTON: f64 = 1.007_276_466_621;

    let t_project = Instant::now();
    let mut projected: Vec<ProjectedFeature> = Vec::new();

    let mut run_order: Vec<usize> = (0..runs.len()).collect();
    run_order.sort_by_key(|&i| &runs[i].name);
    let mut run_ranks = vec![0; runs.len()];
    for (rank, &i) in run_order.iter().enumerate() {
        run_ranks[i] = rank;
    }
    for (run_idx, run) in runs.iter().enumerate() {
        let is_reference = run_idx == alignment.reference_idx;

        for (feat_idx, feat) in run.features.iter().enumerate() {
            if feat.feature.charge == 0 || !feat.combined_score.is_finite() {
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

            if ![ref_mz, ref_rt, ref_im, neutral_mass]
                .iter()
                .all(|x| x.is_finite())
                || neutral_mass <= 0.0
            {
                continue;
            }
            use std::hash::{Hash, Hasher};
            let mut key = std::collections::hash_map::DefaultHasher::new();
            for x in [
                feat.cosine_score,
                feat.isotope_score,
                feat.feature.ppm_error,
                feat.feature.rt_start(),
                feat.feature.rt_end(),
                feat.feature.total_intensity(),
            ] {
                x.to_bits().hash(&mut key);
            }
            feat.neutron_offset.hash(&mut key);
            feat.feature.n_scans_total().hash(&mut key);
            for x in &feat.theoretical_pattern {
                x.to_bits().hash(&mut key);
            }
            projected.push(ProjectedFeature {
                run_idx,
                run_rank: run_ranks[run_idx],
                measurement_key: key.finish(),
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

    projected
}

pub fn build_consensus_candidates(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
) -> Vec<ConsensusFeature> {
    let projected = project_features(runs, alignment);
    let range = runs[alignment.reference_idx].rt_range();
    group_confidence::score_candidates(&projected, &config.consensus, range.1 - range.0, runs.len())
}

/// Research-only scorer/control comparison, independent of production settings.
pub fn audit_consensus(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
    full_controls: bool,
) -> EvidenceAudit {
    let projected = project_features(runs, alignment);
    let range = runs[alignment.reference_idx].rt_range();
    group_confidence::audit_candidates(
        &projected,
        &config.consensus,
        range.1 - range.0,
        runs.len(),
        full_controls,
    )
}

/// Independent RT-permutation diagnostics using the same controls as production.
pub fn audit_permutation_consensus(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
) -> PermutationAudit {
    let projected = project_features(runs, alignment);
    let range = runs[alignment.reference_idx].rt_range();
    group_confidence::audit_permutations(
        &projected,
        &config.consensus,
        range.1 - range.0,
        runs.len(),
    )
}

/// Retain candidates supported by at least two original runs and the group gate.
/// Every accepted group still undergoes independent per-cell LFQ extraction.
pub fn build_consensus(
    runs: &[RunInput],
    alignment: &AlignmentResult,
    config: &LfqConfig,
) -> Vec<ConsensusFeature> {
    config
        .consensus
        .validate()
        .expect("invalid consensus configuration");
    let candidates = build_consensus_candidates(runs, alignment, config);
    let count = candidates.len();
    let retained: Vec<_> = candidates
        .into_iter()
        .filter(|g| {
            g.n_contributing_runs >= 2 && g.group_qvalue <= config.consensus.max_group_qvalue
        })
        .collect();
    log::info!(
        "Consensus confidence: retained {} / {} candidates at group q <= {}",
        retained.len(),
        count,
        config.consensus.max_group_qvalue
    );
    retained
}

fn quality_order(a: &ProjectedFeature, b: &ProjectedFeature) -> std::cmp::Ordering {
    b.combined_score
        .total_cmp(&a.combined_score)
        .then(a.charge.cmp(&b.charge))
        .then(a.neutral_mass.total_cmp(&b.neutral_mass))
        .then(a.ref_rt.total_cmp(&b.ref_rt))
        .then(a.ref_im.total_cmp(&b.ref_im))
        .then(a.run_rank.cmp(&b.run_rank))
        .then(a.measurement_key.cmp(&b.measurement_key))
        .then(a.feature_idx.cmp(&b.feature_idx))
}

struct Group {
    members: Vec<usize>,
    mass: (f64, f64),
    rt: (f64, f64),
    im: Option<(f64, f64)>,
}
impl Group {
    fn new(i: usize, p: &ProjectedFeature) -> Self {
        Self {
            members: vec![i],
            mass: (p.neutral_mass, p.neutral_mass),
            rt: (p.ref_rt, p.ref_rt),
            im: (p.ref_im != 0.0).then_some((p.ref_im, p.ref_im)),
        }
    }
    fn spans(&self, p: &ProjectedFeature) -> GroupSpans {
        let mass = (
            self.mass.0.min(p.neutral_mass),
            self.mass.1.max(p.neutral_mass),
        );
        let rt = (self.rt.0.min(p.ref_rt), self.rt.1.max(p.ref_rt));
        let im = if p.ref_im == 0.0 {
            self.im
        } else {
            Some(self.im.map_or((p.ref_im, p.ref_im), |(lo, hi)| {
                (lo.min(p.ref_im), hi.max(p.ref_im))
            }))
        };
        (mass, rt, im)
    }
}

/// Quality-first clustering: every member must fit the full coordinate span.
/// Fixed seed lookup uses positive mass bits (ordered like positive f64 values).
fn cluster_projected(
    projected: &[ProjectedFeature],
    config: &ConsensusConfig,
    rt_span: f64,
    n_runs: usize,
) -> Vec<ConsensusFeature> {
    assert!(config.mz_ppm.is_finite() && config.mz_ppm >= 0.0);
    assert!(config.rt_window_pct.is_finite() && config.rt_window_pct >= 0.0);
    assert!(config.im_tolerance.is_finite() && config.im_tolerance >= 0.0);
    let rt_tol = config.rt_window_pct * rt_span;
    let ratio = 1.0 + config.mz_ppm / 1e6;
    let mut order: Vec<usize> = (0..projected.len()).collect();
    order.sort_by(|&a, &b| quality_order(&projected[a], &projected[b]));
    let mut index: BTreeMap<(u8, u64), Vec<usize>> = BTreeMap::new();
    let mut groups: Vec<Group> = Vec::new();
    for i in order {
        let p = &projected[i];
        let lo = (p.neutral_mass / ratio).to_bits();
        let hi = (p.neutral_mass * ratio).to_bits();
        let mut best: Option<(f64, usize)> = None;
        for ids in index
            .range((p.charge, lo)..=(p.charge, hi))
            .map(|(_, ids)| ids)
        {
            for &g in ids {
                let (mass, rt, im) = groups[g].spans(p);
                if (mass.1 - mass.0) / mass.0 * 1e6 > config.mz_ppm
                    || rt.1 - rt.0 > rt_tol
                    || im.is_some_and(|(lo, hi)| hi - lo > config.im_tolerance)
                {
                    continue;
                }
                let seed = &projected[groups[g].members[0]];
                let residual = ((p.neutral_mass - seed.neutral_mass) / seed.neutral_mass * 1e6
                    / config.mz_ppm.max(f64::EPSILON))
                .powi(2)
                    + ((p.ref_rt - seed.ref_rt) / rt_tol.max(f64::EPSILON)).powi(2)
                    + if p.ref_im != 0.0 && seed.ref_im != 0.0 {
                        ((p.ref_im - seed.ref_im) / config.im_tolerance.max(f64::EPSILON)).powi(2)
                    } else {
                        0.0
                    };
                if best.is_none_or(|(r, j)| residual < r || (residual == r && g < j)) {
                    best = Some((residual, g));
                }
            }
        }
        if let Some((_, g)) = best {
            let (mass, rt, im) = groups[g].spans(p);
            groups[g].mass = mass;
            groups[g].rt = rt;
            groups[g].im = im;
            groups[g].members.push(i);
        } else {
            index
                .entry((p.charge, p.neutral_mass.to_bits()))
                .or_default()
                .push(groups.len());
            groups.push(Group::new(i, p));
        }
    }
    let mut out: Vec<_> = groups
        .iter()
        .map(|g| emit_group(projected, &g.members, n_runs))
        .collect();
    out.sort_by(|a, b| {
        a.charge
            .cmp(&b.charge)
            .then(a.neutral_mass.total_cmp(&b.neutral_mass))
            .then(a.ref_im.total_cmp(&b.ref_im))
            .then(a.ref_rt.total_cmp(&b.ref_rt))
    });
    out
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
            run_rank: run_idx,
            measurement_key: 0,
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
    pub(super) fn at(run: usize, id: u32, mass: f64, rt: f64, score: f64) -> ProjectedFeature {
        let mut p = proj(run, id, score, mass / 2.0 + 1.007276466621);
        p.neutral_mass = mass;
        p.ref_rt = rt;
        p
    }
    fn cfg() -> ConsensusConfig {
        ConsensusConfig::default()
    }
    fn replicated(
        p: &[ProjectedFeature],
        c: &ConsensusConfig,
        span: f64,
        runs: usize,
    ) -> Vec<ConsensusFeature> {
        cluster_projected(p, c, span, runs)
            .into_iter()
            .filter(|g| g.n_contributing_runs >= 2)
            .collect()
    }
    #[test]
    fn interleaved_low_seed_keeps_its_valid_run_support() {
        let p = vec![
            at(0, 0, 1000.0, 10.0, 0.7),
            at(0, 1, 1000.001, 40.0, 0.9),
            at(1, 0, 1000.002, 10.0, 0.9),
        ];
        let g = replicated(&p, &cfg(), 100.0, 2);
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].n_contributing_runs, 2);
        assert_eq!(g[0].members.len(), 2);
    }
    #[test]
    fn weak_members_are_candidates_without_a_seed_override() {
        let p = vec![
            at(0, 0, 1000.0, 10.0, 0.3),
            at(1, 0, 1000.002, 10.1, 0.2),
            at(0, 1, 1100.0, 20.0, 0.4),
        ];
        let g = replicated(&p, &cfg(), 100.0, 2);
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].seed_combined_score, 0.3);
    }
    #[test]
    fn replicated_support_is_counted_after_interleaved_members_join() {
        let p = vec![
            at(0, 0, 1000.0, 10.0, 0.9),
            at(0, 1, 1000.001, 40.0, 0.9),
            at(1, 0, 1000.002, 10.0, 0.9),
        ];
        let g = replicated(&p, &cfg(), 100.0, 2);
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].n_contributing_runs, 2);
    }
    #[test]
    fn chain_cannot_exceed_mass_or_rt_span() {
        let p = vec![
            at(0, 0, 1000.0, 10.0, 0.99),
            at(1, 0, 1000.015, 11.5, 0.9),
            at(2, 0, 1000.030, 13.0, 0.9),
        ];
        let g = cluster_projected(&p, &cfg(), 100.0, 3);
        assert_eq!(g.len(), 2);
        assert_eq!(g.iter().map(|g| g.members.len()).sum::<usize>(), 3);
        for g in g {
            assert!(g.members.len() <= 2);
        }
    }
    // Diagnostic witnesses: group-level confidence does not reconcile aliases.
    #[test]
    fn rt_outliers_can_form_two_replicated_groups_at_the_same_mass() {
        let p = vec![
            at(0, 0, 1000.0, 10.0, 0.99),
            at(1, 0, 1000.0, 10.1, 0.95),
            at(2, 0, 1000.0, 10.4, 0.90),
            at(3, 0, 1000.0, 10.5, 0.85),
        ];
        let g = cluster_projected(&p, &cfg(), 15.0, 4);
        assert_eq!(g.len(), 2);
        assert!(g.iter().all(|g| g.n_contributing_runs == 2));
        let ownership: std::collections::HashSet<_> =
            g.iter().flat_map(|g| g.members.iter().copied()).collect();
        assert_eq!(ownership.len(), p.len());
        // Unique observations do not imply one group per underlying analyte.
        assert_eq!(g[0].neutral_mass, g[1].neutral_mass);
    }

    #[test]
    fn residual_monoisotope_error_forms_a_separate_group() {
        let c13 = crate::lfq::grid::C13_NEUTRON;
        let p = vec![
            at(0, 0, 1000.0, 10.0, 0.99),
            at(1, 0, 1000.0, 10.0, 0.95),
            at(2, 0, 1000.0 + c13, 10.0, 0.90),
            at(3, 0, 1000.0 + c13, 10.0, 0.85),
        ];
        let g = cluster_projected(&p, &cfg(), 15.0, 4);
        assert_eq!(g.len(), 2);
        assert!(g.iter().all(|g| g.n_contributing_runs == 2));
        assert!((g[1].neutral_mass - g[0].neutral_mass - c13).abs() < 1e-9);
    }

    #[test]
    fn per_run_winner_uses_its_own_score_and_keeps_alternatives() {
        let p = vec![
            at(0, 0, 1000.0, 10.0, 0.99),
            at(1, 10, 1000.001, 10.0, 0.51),
            at(1, 20, 1000.002, 10.0, 0.90),
        ];
        let g = replicated(&p, &cfg(), 100.0, 2);
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].per_run_feature[1], Some(20));
        assert_eq!(g[0].members.len(), 3);
        assert_eq!(g[0].n_contributing_runs, 2);
    }
    #[test]
    fn missing_im_does_not_bridge_incompatible_conformers() {
        let mut p = vec![
            at(0, 0, 1000.0, 10.0, 0.99),
            at(1, 0, 1000.0, 10.0, 0.95),
            at(2, 0, 1000.0, 10.0, 0.9),
        ];
        p[0].ref_im = 1.0;
        p[1].ref_im = 0.0;
        p[2].ref_im = 1.2;
        assert_eq!(cluster_projected(&p, &cfg(), 100.0, 3).len(), 2);
    }
    #[test]
    fn reordering_observations_and_run_indices_preserves_group_coordinates() {
        let mut p = vec![
            at(0, 0, 1000.0, 10.0, 0.9),
            at(1, 0, 1000.015, 11.5, 0.9),
            at(2, 0, 1000.03, 13.0, 0.9),
        ];
        let signature = |g: Vec<ConsensusFeature>| {
            g.iter()
                .map(|f| (f.neutral_mass, f.ref_rt, f.n_contributing_runs))
                .collect::<Vec<_>>()
        };
        let expected = signature(cluster_projected(&p, &cfg(), 100.0, 3));
        p.reverse();
        for x in &mut p {
            x.run_idx = 2 - x.run_idx;
        }
        assert_eq!(expected, signature(cluster_projected(&p, &cfg(), 100.0, 3)));
    }
}
