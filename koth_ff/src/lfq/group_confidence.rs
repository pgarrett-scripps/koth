//! Experimental cross-run evidence and matched-pipeline permutation controls.
//!
//! Scores rank groups; they are not probabilities. Independent RT permutations
//! within each run/charge/quality stratum preserve exact RT density. Ten cohorts
//! undergo the identical grouping and evidence calculation. This stabilizes the
//! tested control estimates but does not establish biological or member-link FDR.
use super::{
    cluster_projected, quality_order, ConsensusConfig, ConsensusFeature, ProjectedFeature,
};

const NULL_COHORTS: usize = 5;
/// Production choice; experimental alternatives are confined to the audit API.
pub const DEFAULT_VARIANT: usize = 0;
pub const VARIANT_NAMES: [&str; 4] = ["average", "average_no_penalty", "tree", "tree_no_penalty"];

/// Research-only comparison; no corresponding user configuration switches.
pub struct EvidenceAudit {
    pub groups: Vec<ConsensusFeature>,
    pub scores: Vec<[f64; 4]>,
    pub q_global: Vec<[f64; 4]>,
    pub q_repeat: Vec<[f64; 4]>,
    pub q_density: Vec<[f64; 4]>,
    pub q_pooled: Vec<[f64; 4]>,
    pub q_conditioned: Vec<[f64; 4]>,
    pub null_group_counts: Vec<usize>,
}

/// Research-only independent-permutation control audit.
pub struct PermutationAudit {
    pub groups: Vec<ConsensusFeature>,
    pub q_first: Vec<f64>,
    pub q_repeat: Vec<f64>,
    pub q_pooled: Vec<f64>,
    pub null_group_counts: Vec<usize>,
}

fn distance(a: &ProjectedFeature, b: &ProjectedFeature, cfg: &ConsensusConfig, rt: f64) -> f64 {
    let scaled = |delta: f64, width: f64| {
        if width > 0.0 {
            (delta / width).powi(2)
        } else if delta == 0.0 {
            0.0
        } else {
            f64::INFINITY
        }
    };
    scaled(
        (a.neutral_mass - b.neutral_mass) / a.neutral_mass.min(b.neutral_mass) * 1e6,
        cfg.mz_ppm,
    ) + scaled(a.ref_rt - b.ref_rt, cfg.rt_window_pct * rt)
        + if a.ref_im != 0.0 && b.ref_im != 0.0 {
            scaled(a.ref_im - b.ref_im, cfg.im_tolerance)
        } else {
            0.0
        }
}

/// Sum pair evidence, divided by n-1 so support grows linearly with run count.
/// A lone excellent feature supplies no cross-run evidence. Each run supplies
/// one primary; with that primary fixed, alternatives only dilute its evidence.
fn evidence_variants<const AUDIT: bool>(
    p: &[ProjectedFeature],
    members: &[usize],
    n_runs: usize,
    cfg: &ConsensusConfig,
    rt: f64,
) -> [f64; 4] {
    let mut primary: Vec<Option<usize>> = vec![None; n_runs];
    let mut quality_sum = vec![0.0; primary.len()];
    for &i in members {
        let run = p[i].run_idx;
        let q = p[i].combined_score.clamp(0.0, 1.0);
        quality_sum[run] += q;
        if primary[run].is_none_or(|j| quality_order(&p[i], &p[j]).is_lt()) {
            primary[run] = Some(i);
        }
    }
    let mut rows: Vec<usize> = primary.into_iter().flatten().collect();
    rows.sort_by_key(|&i| p[i].run_rank);
    if rows.len() < 2 {
        return [0.0; 4];
    }
    let quality = |i: usize| {
        let q = p[i].combined_score.clamp(0.0, 1.0);
        if quality_sum[p[i].run_idx] > 0.0 {
            q * q / quality_sum[p[i].run_idx]
        } else {
            0.0
        }
    };
    let n = rows.len();
    let mut edges = if AUDIT {
        [vec![0.0; n * n], vec![0.0; n * n]]
    } else {
        [Vec::new(), Vec::new()]
    };
    let mut sum = [0.0; 2];
    for (j, &a) in rows.iter().enumerate() {
        for (k, &b) in rows.iter().enumerate().skip(j + 1) {
            let affinity = (-2.0 * distance(&p[a], &p[b], cfg, rt)).exp();
            let weights = [
                (quality(a) * quality(b)).sqrt() * affinity,
                (p[a].combined_score.clamp(0.0, 1.0) * p[b].combined_score.clamp(0.0, 1.0)).sqrt()
                    * affinity,
            ];
            for v in 0..if AUDIT { 2 } else { 1 } {
                sum[v] += weights[v];
                if AUDIT {
                    edges[v][j * n + k] = weights[v];
                    edges[v][k * n + j] = weights[v];
                }
            }
        }
    }
    [
        sum[0] / (n - 1) as f64,
        sum[1] / (n - 1) as f64,
        if AUDIT {
            maximum_tree(&edges[0], n)
        } else {
            0.0
        },
        if AUDIT {
            maximum_tree(&edges[1], n)
        } else {
            0.0
        },
    ]
}

fn evidence_all(
    p: &[ProjectedFeature],
    members: &[usize],
    n_runs: usize,
    cfg: &ConsensusConfig,
    rt: f64,
) -> [f64; 4] {
    evidence_variants::<true>(p, members, n_runs, cfg, rt)
}

/// Maximum spanning tree: n-1 nonnegative edges, not all n(n-1)/2 pairs.
/// Extending the same primary set cannot reduce its score. This is a ranking
/// statistic; correlated runs and edges are not assumed statistically independent.
fn maximum_tree(edges: &[f64], n: usize) -> f64 {
    let mut used = vec![false; n];
    let mut best = vec![0.0f64; n];
    let mut total = 0.0;
    for _ in 0..n {
        let i = (0..n)
            .filter(|&i| !used[i])
            .max_by(|&a, &b| best[a].total_cmp(&best[b]).then(b.cmp(&a)))
            .unwrap();
        used[i] = true;
        total += best[i];
        for j in 0..n {
            if !used[j] {
                best[j] = best[j].max(edges[i * n + j]);
            }
        }
    }
    total
}

fn evidence(
    p: &[ProjectedFeature],
    members: &[usize],
    n_runs: usize,
    cfg: &ConsensusConfig,
    rt: f64,
) -> f64 {
    evidence_variants::<false>(p, members, n_runs, cfg, rt)[DEFAULT_VARIANT]
}

fn evidence_for_groups(
    p: &[ProjectedFeature],
    groups: &mut [ConsensusFeature],
    cfg: &ConsensusConfig,
    rt: f64,
) {
    let index: std::collections::HashMap<_, _> = p
        .iter()
        .enumerate()
        .map(|(i, f)| ((f.run_idx, f.feature_idx), i as u32))
        .collect();
    for g in groups {
        let members: Vec<_> = g.members.iter().map(|m| index[m] as usize).collect();
        g.group_score = evidence(p, &members, g.per_run_feature.len(), cfg, rt);
    }
}

/// Stable integer mixing avoids RNG state and input-order dependence.
fn mix(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

fn shifted_control(
    p: &[ProjectedFeature],
    rt_span: f64,
    replicate: usize,
    local: bool,
) -> Vec<ProjectedFeature> {
    if local {
        return shifted_rank(p, replicate);
    }
    let origin = p.iter().map(|p| p.ref_rt).fold(f64::INFINITY, f64::min);
    let end = p.iter().map(|p| p.ref_rt).fold(f64::NEG_INFINITY, f64::max);
    // Cover projected endpoints too; warps can extend beyond the reference.
    let period = rt_span.max(end - origin);
    p.iter()
        .map(|f| {
            let mut x = f.clone();
            let bits = mix((f.run_rank as u64 + 1).wrapping_mul(0x9e3779b97f4a7c15)
                ^ (replicate as u64 + 1).wrapping_mul(0xd1b54a32d192ed03));
            let fraction = (bits >> 11) as f64 / ((1u64 << 53) as f64);
            x.ref_rt = origin + (f.ref_rt - origin + fraction * period).rem_euclid(period);
            x
        })
        .collect()
}

/// Translate ranks on each run's empirical RT distribution. This preserves its
/// exact RT marginal, including dense regions and ties, without a bin-width knob.
fn shifted_rank(p: &[ProjectedFeature], replicate: usize) -> Vec<ProjectedFeature> {
    shifted_rank_stratified(p, replicate, false)
}

// Research sensitivity check: fixed score deciles preserve quality/RT association.
fn shifted_rank_stratified(
    p: &[ProjectedFeature],
    replicate: usize,
    conditioned: bool,
) -> Vec<ProjectedFeature> {
    let mut out = p.to_vec();
    let mut by_run: std::collections::BTreeMap<(usize, usize), Vec<usize>> =
        std::collections::BTreeMap::new();
    for (i, f) in p.iter().enumerate() {
        let stratum = if conditioned {
            (f.combined_score.clamp(0.0, 1.0) * 10.0).floor().min(9.0) as usize
        } else {
            0
        };
        by_run.entry((f.run_idx, stratum)).or_default().push(i);
    }
    for indices in by_run.values_mut() {
        indices.sort_by(|&a, &b| {
            p[a].ref_rt
                .total_cmp(&p[b].ref_rt)
                .then_with(|| quality_order(&p[a], &p[b]))
        });
        let bits = mix(
            (p[indices[0]].run_rank as u64 + 1).wrapping_mul(0x9e3779b97f4a7c15)
                ^ (replicate as u64 + 1).wrapping_mul(0xd1b54a32d192ed03),
        );
        let shift = ((bits >> 11) as f64 / (1u64 << 53) as f64 * indices.len() as f64) as usize;
        for (j, &i) in indices.iter().enumerate() {
            out[i].ref_rt = p[indices[(j + shift) % indices.len()]].ref_rt;
        }
    }
    out
}

/// Independently assign existing RT positions within run/charge/quality strata.
/// This preserves exact RT density and avoids coherent whole-run displacements.
/// Fixed points are retained: excluding near-target RTs would bias controls.
/// Mass–RT and IM–RT dependence within a stratum is not preserved.
fn permuted_rt(p: &[ProjectedFeature], replicate: usize) -> Vec<ProjectedFeature> {
    let mut out = p.to_vec();
    let mut strata = std::collections::BTreeMap::new();
    for (i, f) in p.iter().enumerate() {
        let bin = (f.combined_score.clamp(0.0, 1.0) * 10.0).floor().min(9.0) as usize;
        strata
            .entry((f.run_idx, f.charge, bin))
            .or_insert_with(Vec::new)
            .push(i);
    }
    for ((_, charge, bin), indices) in &mut strata {
        indices.sort_by(|&a, &b| {
            p[a].ref_rt
                .total_cmp(&p[b].ref_rt)
                .then_with(|| quality_order(&p[a], &p[b]))
        });
        let mut assigned = indices.clone();
        let mut state = mix(
            (p[indices[0]].run_rank as u64 + 1).wrapping_mul(0x9e3779b97f4a7c15)
                ^ (replicate as u64 + 1).wrapping_mul(0xd1b54a32d192ed03)
                ^ (*charge as u64).wrapping_mul(0x94d049bb133111eb)
                ^ (*bin as u64).wrapping_mul(0xbf58476d1ce4e5b9),
        );
        for j in (1..assigned.len()).rev() {
            state = state.wrapping_add(0x9e3779b97f4a7c15);
            // Rejection sampling avoids modulo bias in the bounded draw.
            let bound = (j + 1) as u64;
            let threshold = bound.wrapping_neg() % bound;
            let mut draw = mix(state);
            while draw < threshold {
                state = state.wrapping_add(0x9e3779b97f4a7c15);
                draw = mix(state);
            }
            assigned.swap(j, (draw % bound) as usize);
        }
        for (&i, &source) in indices.iter().zip(&assigned) {
            out[i].ref_rt = p[source].ref_rt;
        }
    }
    out
}

pub(super) fn audit_permutations(
    p: &[ProjectedFeature],
    cfg: &ConsensusConfig,
    rt: f64,
    runs: usize,
) -> PermutationAudit {
    let mut groups = cluster_projected(p, cfg, rt, runs);
    evidence_for_groups(p, &mut groups, cfg, rt);
    let replicated: Vec<_> = groups
        .iter()
        .enumerate()
        .filter(|(_, g)| g.n_contributing_runs >= 2)
        .map(|(i, _)| i)
        .collect();
    let scores: Vec<_> = replicated.iter().map(|&i| groups[i].group_score).collect();
    let valid = !p.is_empty() && rt.is_finite() && rt > 0.0 && runs >= 2;
    let mut null = [Vec::new(), Vec::new()];
    let mut counts = Vec::new();
    if valid {
        for rep in 0..2 * NULL_COHORTS {
            let shuffled = permuted_rt(p, rep);
            let mut gs = cluster_projected(&shuffled, cfg, rt, runs);
            evidence_for_groups(&shuffled, &mut gs, cfg, rt);
            let values: Vec<_> = gs
                .iter()
                .filter(|g| g.n_contributing_runs >= 2)
                .map(|g| g.group_score)
                .collect();
            counts.push(values.len());
            null[rep / NULL_COHORTS].extend(values);
            log::info!(
                "Independent RT permutation {}/{}: {} replicated groups",
                rep + 1,
                2 * NULL_COHORTS,
                counts.last().unwrap()
            );
        }
    }
    let calibrate = |ns: &mut [f64], cohorts| {
        let mut q = vec![1.0; groups.len()];
        if valid {
            for (&i, value) in replicated.iter().zip(group_qvalues(&scores, ns, cohorts)) {
                q[i] = value;
            }
        }
        q
    };
    let q_first = calibrate(&mut null[0], NULL_COHORTS);
    let q_repeat = calibrate(&mut null[1], NULL_COHORTS);
    let mut pooled = null[0].clone();
    pooled.extend_from_slice(&null[1]);
    let q_pooled = calibrate(&mut pooled, 2 * NULL_COHORTS);
    for (g, &q) in groups.iter_mut().zip(&q_pooled) {
        g.group_qvalue = q;
    }
    PermutationAudit {
        groups,
        q_first,
        q_repeat,
        q_pooled,
        null_group_counts: counts,
    }
}

/// Conservative +1 expected-null correction, equal-score blocks and reverse
/// cumulative minimum. No total-count rescaling: each null cohort has the same
/// original feature/run opportunities, but may produce more separate groups.
fn group_qvalues(target: &[f64], null: &mut [f64], cohorts: usize) -> Vec<f64> {
    null.sort_by(|a, b| b.total_cmp(a));
    let mut order: Vec<_> = (0..target.len()).collect();
    order.sort_by(|&a, &b| target[b].total_cmp(&target[a]));
    let mut q = vec![1.0; target.len()];
    let (mut start, mut nd) = (0, 0);
    while start < order.len() {
        let score = target[order[start]];
        // Zero evidence is never a discovery opportunity. Including such rows
        // in the denominator could lower q-values of real scored candidates.
        if score <= 0.0 {
            break;
        }
        let mut end = start + 1;
        while end < order.len() && target[order[end]] == score {
            end += 1;
        }
        while nd < null.len() && null[nd] >= score {
            nd += 1;
        }
        let fdp = ((1.0 + nd as f64 / cohorts as f64) / end as f64).min(1.0);
        for &i in &order[start..end] {
            q[i] = fdp;
        }
        start = end;
    }
    let mut best: f64 = 1.0;
    for &i in order.iter().rev() {
        best = best.min(q[i]);
        q[i] = if target[i] > 0.0 { best } else { 1.0 };
    }
    q
}

pub(super) fn score_candidates(
    p: &[ProjectedFeature],
    cfg: &ConsensusConfig,
    rt: f64,
    runs: usize,
) -> Vec<ConsensusFeature> {
    let audit = audit_permutations(p, cfg, rt, runs);
    let mut retained = [0usize; 3];
    for (i, g) in audit.groups.iter().enumerate() {
        if g.n_contributing_runs >= 2 {
            for (j, q) in [audit.q_first[i], audit.q_repeat[i], audit.q_pooled[i]]
                .into_iter()
                .enumerate()
            {
                retained[j] += usize::from(q <= cfg.max_group_qvalue);
            }
        }
    }
    log::info!("Consensus permutation sensitivity at q <= {}: first five={}, second five={}, pooled ten={}", cfg.max_group_qvalue, retained[0], retained[1], retained[2]);
    audit.groups
}

/// Full controls use disjoint global shifts, marginal RT-rank shifts, and
/// quality-conditioned RT-rank shifts. No null is presumed calibrated.
pub(super) fn audit_candidates(
    p: &[ProjectedFeature],
    cfg: &ConsensusConfig,
    rt: f64,
    runs: usize,
    full_controls: bool,
) -> EvidenceAudit {
    let mut groups = cluster_projected(p, cfg, rt, runs);
    let score = |p: &[ProjectedFeature], gs: &[ConsensusFeature]| {
        let index: std::collections::HashMap<_, _> = p
            .iter()
            .enumerate()
            .map(|(i, f)| ((f.run_idx, f.feature_idx), i))
            .collect();
        gs.iter()
            .map(|g| {
                let ids: Vec<_> = g.members.iter().map(|m| index[m]).collect();
                evidence_all(p, &ids, runs, cfg, rt)
            })
            .collect::<Vec<_>>()
    };
    let scores = score(p, &groups);
    let batches = if full_controls { 5 } else { 1 };
    let mut null: Vec<[Vec<f64>; 4]> = (0..batches)
        .map(|_| std::array::from_fn(|_| Vec::new()))
        .collect();
    let mut counts = Vec::new();
    if !p.is_empty() && rt.is_finite() && rt > 0.0 && runs >= 2 {
        for (batch, values) in null.iter_mut().enumerate() {
            for rep in 0..NULL_COHORTS {
                let shifted = if batch >= 3 {
                    shifted_rank_stratified(p, rep + (batch - 3) * NULL_COHORTS, true)
                } else {
                    shifted_control(
                        p,
                        rt,
                        if batch == 1 { rep + NULL_COHORTS } else { rep },
                        batch == 2,
                    )
                };
                let gs = cluster_projected(&shifted, cfg, rt, runs);
                let ss = score(&shifted, &gs);
                let mut count = 0;
                for (g, s) in gs.iter().zip(ss) {
                    if g.n_contributing_runs >= 2 {
                        count += 1;
                        for v in 0..4 {
                            values[v].push(s[v]);
                        }
                    }
                }
                counts.push(count);
                log::info!(
                    "audit control batch {} replicate {}: {} groups",
                    batch,
                    rep,
                    count
                );
            }
        }
    }
    let replicated: Vec<_> = groups
        .iter()
        .enumerate()
        .filter(|(_, g)| g.n_contributing_runs >= 2)
        .map(|(i, _)| i)
        .collect();
    let calibrate = |ns: &mut [Vec<f64>; 4], cohorts: usize| {
        let mut out = vec![[1.0; 4]; groups.len()];
        if p.is_empty() || !rt.is_finite() || rt <= 0.0 || runs < 2 {
            return out;
        }
        for v in 0..4 {
            let target: Vec<_> = replicated.iter().map(|&i| scores[i][v]).collect();
            let q = group_qvalues(&target, &mut ns[v], cohorts);
            for (&i, q) in replicated.iter().zip(q) {
                out[i][v] = q;
            }
        }
        out
    };
    let q_global = calibrate(&mut null[0], NULL_COHORTS);
    let (q_repeat, q_density, q_pooled) = if full_controls {
        let repeat = calibrate(&mut null[1], NULL_COHORTS);
        let density = calibrate(&mut null[2], NULL_COHORTS);
        let mut pooled: [Vec<f64>; 4] = std::array::from_fn(|v| {
            let mut x = null[0][v].clone();
            x.extend_from_slice(&null[1][v]);
            x
        });
        (repeat, density, calibrate(&mut pooled, 2 * NULL_COHORTS))
    } else {
        (
            vec![[f64::NAN; 4]; groups.len()],
            vec![[f64::NAN; 4]; groups.len()],
            vec![[f64::NAN; 4]; groups.len()],
        )
    };
    let q_conditioned = if full_controls {
        let mut pooled = std::array::from_fn(|v| {
            let mut x = null[3][v].clone();
            x.extend_from_slice(&null[4][v]);
            x
        });
        calibrate(&mut pooled, 2 * NULL_COHORTS)
    } else {
        vec![[f64::NAN; 4]; groups.len()]
    };
    for (i, g) in groups.iter_mut().enumerate() {
        g.group_score = scores[i][DEFAULT_VARIANT];
        g.group_qvalue = q_global[i][DEFAULT_VARIANT];
    }
    EvidenceAudit {
        groups,
        scores,
        q_global,
        q_repeat,
        q_density,
        q_pooled,
        q_conditioned,
        null_group_counts: counts,
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::at;
    use super::*;

    #[test]
    fn independent_permutation_preserves_strata_and_is_order_invariant() {
        let p: Vec<_> = (0..200)
            .flat_map(|id| {
                (0..3).map(move |r| {
                    let mut f = at(r, id, 800.0 + id as f64, id as f64, (id % 10) as f64 / 10.0);
                    f.charge = 2 + (id % 2) as u8;
                    f
                })
            })
            .collect();
        let null = permuted_rt(&p, 0);
        let stratum_rts = |features: &[ProjectedFeature]| {
            let mut out = std::collections::BTreeMap::new();
            for f in features {
                out.entry((f.run_rank, f.charge, (f.combined_score * 10.0) as usize))
                    .or_insert_with(Vec::new)
                    .push(f.ref_rt);
            }
            for values in out.values_mut() {
                values.sort_by(f64::total_cmp);
            }
            out
        };
        assert_eq!(stratum_rts(&p), stratum_rts(&null));
        for (a, b) in p.iter().zip(&null) {
            assert_eq!(
                (a.neutral_mass, a.combined_score, a.ref_im, a.feature_idx),
                (b.neutral_mass, b.combined_score, b.ref_im, b.feature_idx)
            );
        }
        let signature = |features: &[ProjectedFeature]| {
            let mut rows: Vec<_> = features
                .iter()
                .map(|f| (f.run_rank, f.feature_idx, f.ref_rt))
                .collect();
            rows.sort_by_key(|f| (f.0, f.1));
            rows
        };
        assert_ne!(signature(&p), signature(&null));
        assert_ne!(signature(&null), signature(&permuted_rt(&p, 1)));
        let mut reordered = p.clone();
        reordered.reverse();
        for f in &mut reordered {
            f.run_idx = 2 - f.run_idx;
        }
        assert_eq!(signature(&null), signature(&permuted_rt(&reordered, 0)));
    }

    #[test]
    fn permutation_cannot_claim_confidence_when_rt_has_no_discriminating_power() {
        let p: Vec<_> = (0..100)
            .flat_map(|id| (0..4).map(move |run| at(run, id, 800.0 + id as f64, 30.0, 0.3)))
            .collect();
        let audit = audit_permutations(&p, &ConsensusConfig::default(), 100.0, 4);
        assert_eq!(audit.groups.len(), 100);
        assert!(audit.groups.iter().all(|g| g.group_qvalue == 1.0));
    }

    #[test]
    fn permutation_independent_noise_is_not_rescued() {
        let p: Vec<_> = (0..2000)
            .flat_map(|id| {
                (0..4).map(move |run| {
                    let bits = mix((id as u64 * 4 + run as u64 + 1) ^ 987654321);
                    at(
                        run,
                        id,
                        800.0 + id as f64 * 0.1,
                        100.0 * (bits >> 11) as f64 / (1u64 << 53) as f64,
                        0.3,
                    )
                })
            })
            .collect();
        let audit = audit_permutations(&p, &ConsensusConfig::default(), 100.0, 4);
        assert!(audit.groups.iter().all(|g| g.group_qvalue > 0.05));
    }

    #[test]
    fn tree_keeps_strong_core_when_a_weak_run_joins() {
        let cfg = ConsensusConfig::default();
        let mut p = vec![at(0, 0, 1000.0, 10.0, 0.9), at(1, 0, 1000.0, 10.0, 0.9)];
        let pair = evidence_all(&p, &[0, 1], 3, &cfg, 100.0);
        p.push(at(2, 0, 1000.0, 10.0, 0.01));
        let triple = evidence_all(&p, &[0, 1, 2], 3, &cfg, 100.0);
        assert_eq!(evidence(&p, &[0, 1, 2], 3, &cfg, 100.0), triple[0]);
        assert!(triple[0] < pair[0]);
        assert!(triple[2] >= pair[2]);
        assert!((pair[2] - 0.9).abs() < 1e-12);
    }

    #[test]
    fn rank_controls_preserve_each_runs_rt_marginal() {
        let p: Vec<_> = (0..20)
            .flat_map(|id| (0..3).map(move |r| at(r, id, 800.0 + id as f64, (id / 4) as f64, 0.3)))
            .collect();
        let null = shifted_rank(&p, 0);
        for run in 0..3 {
            let values = |p: &[ProjectedFeature]| {
                let mut v: Vec<_> = p
                    .iter()
                    .filter(|f| f.run_idx == run)
                    .map(|f| f.ref_rt)
                    .collect();
                v.sort_by(f64::total_cmp);
                v
            };
            assert_eq!(values(&p), values(&null));
        }
        for (a, b) in p.iter().zip(&null) {
            assert_eq!(a.neutral_mass, b.neutral_mass);
            assert_eq!(a.combined_score, b.combined_score);
        }
    }

    #[test]
    fn conditioned_controls_preserve_rt_within_quality_strata() {
        let p: Vec<_> = (0..200)
            .flat_map(|id| {
                (0..3)
                    .map(move |r| at(r, id, 800.0 + id as f64, id as f64, (id % 10) as f64 / 10.0))
            })
            .collect();
        let shifted = shifted_rank_stratified(&p, 2, true);
        for (a, b) in p.iter().zip(&shifted) {
            assert_eq!(a.combined_score, b.combined_score);
            assert_eq!(a.neutral_mass, b.neutral_mass);
        }
        for run in 0..3 {
            for bin in 0..10 {
                let marginal = |features: &[ProjectedFeature]| {
                    let mut rts: Vec<_> = features
                        .iter()
                        .filter(|f| f.run_idx == run && f.combined_score == bin as f64 / 10.0)
                        .map(|f| f.ref_rt)
                        .collect();
                    rts.sort_by(f64::total_cmp);
                    rts
                };
                assert_eq!(marginal(&p), marginal(&shifted));
            }
        }
    }

    #[test]
    #[ignore = "bounded null-model audit; run explicitly in release mode"]
    fn refinement_null_stress_audit() {
        for scenario in [
            "uniform",
            "clustered",
            "recurring_artifact",
            "absent_members",
        ] {
            let mut totals = vec![vec![[0usize; 4]; 4]; 3];
            let mut permutation_totals = [0usize; 4];
            for trial in 0..5u64 {
                let mut p = Vec::new();
                for id in 0..100 {
                    for run in 0..6 {
                        let bits = mix(id as u64 * 6 + run as u64 + trial * 100000);
                        let u = (bits >> 11) as f64 / (1u64 << 53) as f64;
                        let absent = scenario == "absent_members" && run >= 3;
                        p.push(at(
                            run,
                            id,
                            800.0 + id as f64,
                            if absent {
                                100.0 * u
                            } else {
                                10.0 + id as f64 * 0.7 + run as f64 * 0.005
                            },
                            if absent {
                                0.3
                            } else if scenario == "absent_members" {
                                0.9
                            } else {
                                0.35
                            },
                        ));
                    }
                }
                for id in 100..600 {
                    for run in 0..6 {
                        let bits = mix(id as u64 * 6 + run as u64 + trial * 100000);
                        let u = (bits >> 11) as f64 / (1u64 << 53) as f64;
                        let rt = match scenario {
                            "clustered" => 35.0 + 2.0 * u,
                            "recurring_artifact" => 10.0 + (id % 100) as f64 * 0.7 + u * 0.01,
                            _ => 100.0 * u,
                        };
                        p.push(at(run, id, 1800.0 + id as f64, rt, 0.1 + 0.8 * u));
                    }
                }
                let perm = audit_permutations(&p, &ConsensusConfig::default(), 100.0, 6);
                for g in perm
                    .groups
                    .iter()
                    .filter(|g| g.n_contributing_runs >= 2 && g.group_qvalue <= 0.05)
                {
                    let false_group = g.neutral_mass >= 1000.0;
                    permutation_totals[usize::from(false_group)] += 1;
                    if scenario == "absent_members" && !false_group {
                        permutation_totals[2] += g.members.iter().filter(|&&(r, _)| r >= 3).count();
                        permutation_totals[3] += g.members.len();
                    }
                }
                let audit = audit_candidates(&p, &ConsensusConfig::default(), 100.0, 6, true);
                for (control, qs) in [&audit.q_global, &audit.q_density, &audit.q_conditioned]
                    .iter()
                    .enumerate()
                {
                    for (g, q) in audit.groups.iter().zip(qs.iter()) {
                        for v in 0..4 {
                            if q[v] <= 0.05 && g.n_contributing_runs >= 2 {
                                let false_group = g.neutral_mass >= 1000.0;
                                totals[control][v][usize::from(false_group)] += 1;
                                if scenario == "absent_members" && !false_group {
                                    totals[control][v][2] +=
                                        g.members.iter().filter(|&&(r, _)| r >= 3).count();
                                    totals[control][v][3] += g.members.len();
                                }
                            }
                        }
                    }
                }
            }
            eprintln!(
                "STRESS {}",
                serde_json::json!({"scenario":scenario,"trials":5,"independent_permutation":permutation_totals,
                "columns":["true_groups","false_groups_or_artifacts","absent_member_links","all_links_in_true_groups"],
                "variants":VARIANT_NAMES,"global":totals[0],"density":totals[1],"conditioned":totals[2]})
            );
        }
    }

    #[test]
    fn ties_have_equal_q_and_finite_sample_floor() {
        let q = group_qvalues(&[3.0, 3.0, 2.0, 0.0], &mut [3.0; 5], 5);
        assert_eq!(q[0], q[1]);
        assert!(q[0] > 0.0);
        assert!(q[0] <= q[2]);
        assert_eq!(q[3], 1.0);
        assert_eq!(group_qvalues(&[1.0], &mut [], 5), vec![1.0]);
        assert_eq!(
            group_qvalues(&[1.0, 0.0, 0.0, 0.0], &mut [], 5),
            vec![1.0; 4]
        );
    }

    #[test]
    fn consistent_weak_population_survives_and_singletons_do_not() {
        let cfg = ConsensusConfig::default();
        let mut p = Vec::new();
        for id in 0..100 {
            for run in 0..4 {
                p.push(at(
                    run,
                    id,
                    800.0 + id as f64,
                    10.0 + id as f64 * 0.7 + run as f64 * 0.01,
                    0.3,
                ));
            }
        }
        p.push(at(0, 100, 1500.0, 30.0, 0.99));
        let groups = score_candidates(&p, &cfg, 100.0, 4);
        assert_eq!(
            groups.iter().filter(|g| g.group_qvalue <= 0.05).count(),
            100
        );
        let singleton = groups.iter().find(|g| g.n_contributing_runs == 1).unwrap();
        assert_eq!(singleton.group_qvalue, 1.0);
        let first: Vec<_> = groups
            .iter()
            .map(|g| (g.neutral_mass, g.group_score, g.group_qvalue))
            .collect();
        p.reverse();
        for f in &mut p {
            f.run_idx = 3 - f.run_idx;
        }
        let second: Vec<_> = score_candidates(&p, &cfg, 100.0, 4)
            .iter()
            .map(|g| (g.neutral_mass, g.group_score, g.group_qvalue))
            .collect();
        assert_eq!(first, second);
    }

    #[test]
    fn duplicate_observations_do_not_create_support() {
        let p = vec![at(0, 0, 1000.0, 10.0, 0.9), at(0, 1, 1000.0, 10.0, 0.9)];
        let g = score_candidates(&p, &ConsensusConfig::default(), 100.0, 2);
        assert_eq!(g[0].group_score, 0.0);
        assert_eq!(g[0].group_qvalue, 1.0);
    }

    #[test]
    fn tight_agreement_beats_near_boundary_match_and_ambiguity() {
        let cfg = ConsensusConfig::default();
        let p = vec![at(0, 0, 1000.0, 10.0, 0.3), at(1, 0, 1000.001, 10.01, 0.3)];
        let score = |p: &[ProjectedFeature]| {
            let mut g = cluster_projected(p, &cfg, 100.0, 2);
            evidence_for_groups(p, &mut g, &cfg, 100.0);
            g[0].group_score
        };
        let tight = score(&p);
        let mut loose = p.clone();
        loose[1].ref_rt = 11.9;
        assert!(tight > score(&loose));
        let mut ambiguous = p.clone();
        let mut alt = p[1].clone();
        alt.feature_idx = 1;
        ambiguous.push(alt);
        assert!(tight > score(&ambiguous));
    }

    #[test]
    fn independent_uniform_noise_is_not_rescued() {
        let mut p = Vec::new();
        for id in 0..2000 {
            for run in 0..4 {
                let bits = mix((id as u64 * 4 + run as u64 + 1) ^ 987654321);
                let rt = 100.0 * (bits >> 11) as f64 / (1u64 << 53) as f64;
                p.push(at(run, id, 800.0 + id as f64 * 0.1, rt, 0.3));
            }
        }
        let groups = score_candidates(&p, &ConsensusConfig::default(), 100.0, 4);
        assert_eq!(groups.iter().filter(|g| g.group_qvalue <= 0.05).count(), 0);
    }

    #[test]
    fn mixed_population_audit_recovers_weak_signals_without_excess_false_groups() {
        let (mut recovered, mut false_groups) = (0, 0);
        for trial in 0..20u64 {
            let mut p = Vec::new();
            for id in 0..100 {
                for run in 0..6 {
                    p.push(at(
                        run,
                        id,
                        800.0 + id as f64,
                        40.0 + 0.01 * run as f64,
                        0.35,
                    ));
                }
            }
            for id in 100..1100 {
                for run in 0..6 {
                    let bits = mix((id as u64 * 6 + run as u64) ^ mix(trial + 12345));
                    let rt = 100.0 * (bits >> 11) as f64 / (1u64 << 53) as f64;
                    let quality = 0.1 + 0.8 * (mix(bits) >> 11) as f64 / (1u64 << 53) as f64;
                    p.push(at(run, id, 1800.0 + id as f64, rt, quality));
                }
            }
            for g in score_candidates(&p, &ConsensusConfig::default(), 100.0, 6) {
                if g.group_qvalue <= 0.05 {
                    if g.neutral_mass < 1000.0 {
                        recovered += 1;
                    } else {
                        false_groups += 1;
                    }
                }
            }
        }
        eprintln!(
            "mixed uniform-noise audit: {recovered}/2000 weak signals, {false_groups} false groups"
        );
        assert!(recovered >= 1900);
        assert!((false_groups as f64) / (recovered + false_groups) as f64 <= 0.075);
    }
}
