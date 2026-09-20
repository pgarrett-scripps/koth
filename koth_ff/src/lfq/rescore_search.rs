//! Experimental search-guided scoring. Legacy scoring remains unchanged.
use super::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SearchScoring {
    #[default]
    Legacy,
    PeptideGrouped,
    ChargeStratified,
    ChargeSigned,
    ChargeSignedQuality,
}

/// Stable peptide-only fold: all runs, charges, targets and paired decoys stay
/// together. FNV-1a is explicitly specified, independent of HashMap's seed.
pub fn peptide_fold(peptide: &str) -> usize {
    let hash = peptide.bytes().fold(0xcbf29ce484222325u64, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x100000001b3)
    });
    (hash % N_FOLDS as u64) as usize
}

fn charge_bin(charge: u8) -> u8 {
    match charge {
        2 => 2,
        3 => 3,
        _ => 0,
    }
}

/// Count whole score ties together; a +1 decoy correction prevents unopposed
/// targets (including a stratum with no decoys) from receiving zero q-values.
fn conservative_qvalues(scores: &[f64], decoy: &[bool]) -> Vec<f64> {
    if !decoy.iter().any(|&d| d) {
        return vec![1.0; scores.len()];
    }
    let mut order: Vec<_> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
    let mut q = vec![1.0; scores.len()];
    let (mut nt, mut nd, mut start) = (0usize, 0usize, 0usize);
    let mut blocks = Vec::new();
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() && scores[order[end]] == scores[order[start]] {
            end += 1;
        }
        for &i in &order[start..end] {
            if decoy[i] {
                nd += 1;
            } else {
                nt += 1;
            }
        }
        blocks.push((start, end, ((nd + 1) as f64 / nt.max(1) as f64).min(1.0)));
        start = end;
    }
    let mut minimum: f64 = 1.0;
    for (start, end, fdr) in blocks.into_iter().rev() {
        minimum = minimum.min(fdr);
        for &i in &order[start..end] {
            if !decoy[i] {
                q[i] = minimum;
            }
        }
    }
    q
}

struct Transform<const N: usize = NF> {
    median: [f64; N],
    mean: [f64; N],
    sd: [f64; N],
}
impl<const N: usize> Transform<N> {
    /// All preprocessing is learned from the training fold only.
    fn fit(raw: &[[f64; N]], train: &[usize]) -> Self {
        let mut t = Self {
            median: [0.0; N],
            mean: [0.0; N],
            sd: [1.0; N],
        };
        for k in 0..N {
            let mut values: Vec<_> = train
                .iter()
                .map(|&i| raw[i][k])
                .filter(|v| v.is_finite())
                .collect();
            values.sort_by(f64::total_cmp);
            t.median[k] = values.get(values.len() / 2).copied().unwrap_or(0.0);
            let value = |i: usize| {
                if raw[i][k].is_finite() {
                    raw[i][k]
                } else {
                    t.median[k]
                }
            };
            t.mean[k] = train.iter().map(|&i| value(i)).sum::<f64>() / train.len().max(1) as f64;
            t.sd[k] = (train
                .iter()
                .map(|&i| (value(i) - t.mean[k]).powi(2))
                .sum::<f64>()
                / train.len().max(1) as f64)
                .sqrt();
            if t.sd[k] < 1e-12 {
                t.sd[k] = 1.0;
            }
        }
        t
    }
    fn apply(&self, raw: &[f64; N]) -> [f64; N] {
        std::array::from_fn(|k| {
            let value = if raw[k].is_finite() {
                raw[k]
            } else {
                self.median[k]
            };
            (value - self.mean[k]) / self.sd[k]
        })
    }
}

fn supported(entries: &[&LfqEntry], peptides: &[String]) -> bool {
    (0..N_FOLDS).all(|fold| {
        let training: Vec<_> = entries
            .iter()
            .filter(|e| peptide_fold(&peptides[e.feature_idx]) != fold)
            .collect();
        let targets = training.iter().filter(|e| !e.is_decoy).count();
        let decoys = training.len() - targets;
        let groups: BTreeSet<_> = training.iter().map(|e| &peptides[e.feature_idx]).collect();
        targets >= 100 && decoys >= 100 && groups.len() >= 20
    })
}

fn search_features<const N: usize>(e: &LfqEntry, mode: SearchScoring) -> [f64; N] {
    let mut base = feat_of(e);
    if matches!(
        mode,
        SearchScoring::ChargeSigned | SearchScoring::ChargeSignedQuality
    ) {
        base[0] = if e.observed_mz.is_finite() && e.expected_mz > 0.0 {
            (e.observed_mz - e.expected_mz) / e.expected_mz * 1e6
        } else {
            f64::NAN
        };
        base[1] = if e.apex_rt.is_finite() && e.expected_rt.is_finite() {
            e.apex_rt - e.expected_rt
        } else {
            f64::NAN
        };
    }
    std::array::from_fn(|k| match k {
        0..=4 => base[k],
        5 => f64::from(e.n_isotopes_found).ln_1p(),
        6 => {
            if e.peak_width_rt.is_finite() {
                e.peak_width_rt.max(0.0).ln_1p()
            } else {
                f64::NAN
            }
        }
        7 => f64::from(e.preceding_signal_fraction),
        _ => unreachable!("unsupported search feature dimension"),
    })
}

fn score_stratum<const N: usize>(
    entries: &[&LfqEntry],
    peptides: &[String],
    mode: SearchScoring,
) -> HashMap<(usize, usize), f64> {
    let raw: Vec<_> = entries
        .iter()
        .map(|e| search_features::<N>(e, mode))
        .collect();
    let dec: Vec<_> = entries.iter().map(|e| e.is_decoy).collect();
    let fold: Vec<_> = entries
        .iter()
        .map(|e| peptide_fold(&peptides[e.feature_idx]))
        .collect();
    let hybrid: Vec<_> = entries
        .iter()
        .map(|e| {
            if e.hybrid_score.is_finite() {
                e.hybrid_score as f64
            } else {
                0.0
            }
        })
        .collect();
    let mut heldout = hybrid.clone();
    // Use hybrid for the ENTIRE stratum if any fold cannot be fitted: mixing
    // hybrid and likelihood-ratio score scales would invalidate pooled ranking.
    let mut complete = true;
    for f in 0..N_FOLDS {
        let train: Vec<_> = (0..entries.len()).filter(|&i| fold[i] != f).collect();
        let test: Vec<_> = (0..entries.len()).filter(|&i| fold[i] == f).collect();
        if test.is_empty() {
            continue;
        }
        let transform = Transform::fit(&raw, &train);
        let feat: Vec<_> = train.iter().map(|&i| transform.apply(&raw[i])).collect();
        let train_dec: Vec<_> = train.iter().map(|&i| dec[i]).collect();
        let neg: Vec<_> = (0..train.len()).filter(|&i| train_dec[i]).collect();
        let mut scores: Vec<_> = train.iter().map(|&i| hybrid[i]).collect();
        let mut model = None;
        for _ in 0..N_ITER {
            let q = conservative_qvalues(&scores, &train_dec);
            let mut pos: Vec<_> = (0..train.len())
                .filter(|&i| !train_dec[i] && q[i] <= TRAIN_FDR)
                .collect();
            if pos.len() < MIN_POS {
                pos = (0..train.len()).filter(|&i| !train_dec[i]).collect();
                pos.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]).then(a.cmp(&b)));
                pos.truncate(MIN_POS);
            }
            let Some(m) = fit_qda(&feat, &pos, &neg) else {
                break;
            };
            scores = feat.iter().map(|x| m.score(x)).collect();
            model = Some(m);
        }
        if let Some(m) = model {
            for i in test {
                heldout[i] = m.score(&transform.apply(&raw[i]));
            }
        } else {
            complete = false;
        }
    }
    if !complete || heldout.iter().any(|s| !s.is_finite()) {
        heldout = hybrid;
    }
    let q = conservative_qvalues(&heldout, &dec);
    entries
        .iter()
        .enumerate()
        .filter(|(_, e)| !e.is_decoy)
        .map(|(i, e)| ((e.feature_idx, e.run_idx), q[i]))
        .collect()
}

/// Call separately for direct and transferred cells. Metadata is indexed by
/// feature_idx and shared by each target/decoy pair. Sparse charge bins cause
/// the whole evidence class to use peptide-grouped pooled scoring.
pub fn compute(
    entries: &[LfqEntry],
    peptides: &[String],
    charges: &[u8],
    mode: SearchScoring,
) -> HashMap<(usize, usize), f64> {
    if mode == SearchScoring::Legacy {
        return compute_qvalues_qda(entries);
    }
    assert_eq!(peptides.len(), charges.len());
    let positive: Vec<_> = entries
        .iter()
        .filter(|e| e.intensity.is_finite() && e.intensity > 0.0)
        .collect();
    let mut strata: BTreeMap<u8, Vec<&LfqEntry>> = BTreeMap::new();
    for e in &positive {
        strata
            .entry(charge_bin(charges[e.feature_idx]))
            .or_default()
            .push(e);
    }
    let mut out: HashMap<_, _> = entries
        .iter()
        .filter(|e| !e.is_decoy)
        .map(|e| ((e.feature_idx, e.run_idx), 1.0))
        .collect();
    let score = |subset: &[&LfqEntry]| {
        if mode == SearchScoring::ChargeSignedQuality {
            score_stratum::<8>(subset, peptides, mode)
        } else {
            score_stratum::<NF>(subset, peptides, mode)
        }
    };
    if matches!(
        mode,
        SearchScoring::ChargeStratified
            | SearchScoring::ChargeSigned
            | SearchScoring::ChargeSignedQuality
    ) && strata.values().all(|s| supported(s, peptides))
    {
        for (charge, subset) in strata {
            log::info!(
                "search QDA charge bin {charge}: {} positive cells",
                subset.len()
            );
            out.extend(score(&subset));
        }
    } else {
        log::info!(
            "search QDA peptide-grouped pooled scoring: {} positive cells",
            positive.len()
        );
        out.extend(score(&positive));
    }
    out
}

#[cfg(test)]
#[path = "rescore_search_tests.rs"]
mod tests;
