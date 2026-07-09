//! Semi-supervised QDA rescorer for match-between-runs target-decoy FDR.
//!
//! Replaces the plain "rank all cells by `hybrid_score`" TDC ([`super::tdc`])
//! with a learned discriminant over five symmetric per-cell features:
//!   0. |ppm error|   — mass distance from the cell's own predicted m/z
//!   1. |RT diff|     — RT distance from the cell's own predicted RT
//!   2. Bhattacharyya — observed vs theoretical isotope PATTERN match
//!   3. coelution     — cosine between the matched isotopes' XIC traces
//!   4. |IM delta|    — ion-mobility distance from the cell's own predicted IM
//!                      (inert on Orbitrap, where IM is absent → imputed)
//!
//! Protocol (Percolator-style, but with a QDA core):
//!   * Competition = ALL target cells (detected + match-between-runs) + decoy
//!     cells (all with signal). Every target is scored and gets a q, so a
//!     background-contaminated cell in a depleted well is gatable regardless of
//!     whether it was detected or transferred — required for fold-change rescue
//!     on large-dynamic-range designs (e.g. timsTOF HYE).
//!   * 3-fold cross-validation by `feature_idx` (a cell is never scored by a
//!     model trained on its own feature).
//!   * Within each training fold, iterate: q-values by the current score →
//!     confident targets (train-FDR ≤ 1%, or the top-N by score as a fallback)
//!     as positives, decoys as negatives → fit QDA → rescore.
//!   * Final q-values are computed by TDC over the held-out QDA scores.
//!
//! QDA (vs LDA) fits a SEPARATE covariance per class, which matches this
//! problem: real transfers cluster tightly near the origin (low ppm/RT, high
//! pattern/coelution) while decoys are a diffuse cloud — a shared-covariance
//! assumption is wrong here. Deterministic: no RNG, all ties broken by index.

use std::collections::HashMap;

use super::LfqEntry;

/// Number of discriminant features.
const NF: usize = 5;
/// Semi-supervised refinement iterations per fold.
const N_ITER: usize = 10;
/// Cross-validation folds (partitioned by `feature_idx`).
const N_FOLDS: usize = 3;
/// Train-FDR for selecting confident positive examples.
const TRAIN_FDR: f64 = 0.01;
/// Fallback: if fewer than this many targets clear the train-FDR, take the
/// top-N targets by the current score instead (keeps early iterations stable).
const MIN_POS: usize = 1000;
/// Covariance ridge (features are standardised, so unit-scale). Guarantees a
/// positive-definite matrix for the Cholesky factorisation.
const RIDGE: f64 = 1e-2;

fn sub(a: &[f64; NF], b: &[f64; NF]) -> [f64; NF] {
    let mut o = [0.0; NF];
    for k in 0..NF {
        o[k] = a[k] - b[k];
    }
    o
}

/// Lower-triangular Cholesky factor of a symmetric matrix; `None` if not
/// positive-definite.
fn cholesky(a: &[[f64; NF]; NF]) -> Option<[[f64; NF]; NF]> {
    let mut l = [[0.0f64; NF]; NF];
    for i in 0..NF {
        for j in 0..=i {
            let mut s = a[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i][j] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    Some(l)
}

/// Squared Mahalanobis distance `(x-μ)ᵀ Σ⁻¹ (x-μ)` via forward substitution on
/// the Cholesky factor: solve `L y = diff`, return `yᵀy`.
fn mahalanobis(l: &[[f64; NF]; NF], diff: [f64; NF]) -> f64 {
    let mut y = [0.0f64; NF];
    for i in 0..NF {
        let mut s = diff[i];
        for k in 0..i {
            s -= l[i][k] * y[k];
        }
        y[i] = s / l[i][i];
    }
    y.iter().map(|v| v * v).sum()
}

fn logdet_from_l(l: &[[f64; NF]; NF]) -> f64 {
    2.0 * (0..NF).map(|i| l[i][i].ln()).sum::<f64>()
}

/// Mean, Cholesky factor of the (ridge-regularised) covariance, and log-det for
/// one class. `None` if too few rows to estimate a covariance.
fn fit_gaussian(feat: &[[f64; NF]], idx: &[usize]) -> Option<([f64; NF], [[f64; NF]; NF], f64)> {
    let n = idx.len();
    if n < NF + 2 {
        return None;
    }
    let mut mean = [0.0f64; NF];
    for &i in idx {
        for k in 0..NF {
            mean[k] += feat[i][k];
        }
    }
    for k in 0..NF {
        mean[k] /= n as f64;
    }
    let mut cov = [[0.0f64; NF]; NF];
    for &i in idx {
        let d = sub(&feat[i], &mean);
        for a in 0..NF {
            for b in 0..NF {
                cov[a][b] += d[a] * d[b];
            }
        }
    }
    let denom = (n - 1) as f64;
    for a in 0..NF {
        for b in 0..NF {
            cov[a][b] /= denom;
        }
        cov[a][a] += RIDGE;
    }
    let l = cholesky(&cov)?;
    let logdet = logdet_from_l(&l);
    Some((mean, l, logdet))
}

/// QDA target-vs-decoy log-likelihood-ratio model over standardised features.
struct Qda {
    mean_t: [f64; NF],
    l_t: [[f64; NF]; NF],
    logdet_t: f64,
    mean_d: [f64; NF],
    l_d: [[f64; NF]; NF],
    logdet_d: f64,
    /// Orientation multiplier so higher score = more target-like.
    sign: f64,
}

impl Qda {
    fn score(&self, x: &[f64; NF]) -> f64 {
        let mt = mahalanobis(&self.l_t, sub(x, &self.mean_t));
        let md = mahalanobis(&self.l_d, sub(x, &self.mean_d));
        self.sign * (-0.5 * mt - 0.5 * self.logdet_t + 0.5 * md + 0.5 * self.logdet_d)
    }
}

fn fit_qda(feat: &[[f64; NF]], pos: &[usize], neg: &[usize]) -> Option<Qda> {
    let (mean_t, l_t, logdet_t) = fit_gaussian(feat, pos)?;
    let (mean_d, l_d, logdet_d) = fit_gaussian(feat, neg)?;
    let mut m = Qda {
        mean_t,
        l_t,
        logdet_t,
        mean_d,
        l_d,
        logdet_d,
        sign: 1.0,
    };
    let pos_mean = pos.iter().map(|&i| m.score(&feat[i])).sum::<f64>() / pos.len().max(1) as f64;
    let neg_mean = neg.iter().map(|&i| m.score(&feat[i])).sum::<f64>() / neg.len().max(1) as f64;
    if pos_mean < neg_mean {
        m.sign = -1.0;
    }
    Some(m)
}

/// Running target-decoy q-values from scores (higher = target-like). Returns a
/// q per row; decoy rows keep 1.0 (unused). Ties broken by index → deterministic.
fn qvalues(scores: &[f64], is_decoy: &[bool]) -> Vec<f64> {
    let n = scores.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        scores[b]
            .partial_cmp(&scores[a])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(&b))
    });
    let mut nt = 0usize;
    let mut nd = 0usize;
    let mut fdr = vec![1.0f64; n];
    let mut targets_in_rank: Vec<usize> = Vec::new();
    for &i in &order {
        if is_decoy[i] {
            nd += 1;
        } else {
            nt += 1;
            fdr[i] = nd as f64 / nt as f64;
            targets_in_rank.push(i);
        }
    }
    // Monotonise: each target's q = min FDR among itself and all worse-scoring targets.
    let mut min_q = 1.0f64;
    for &i in targets_in_rank.iter().rev() {
        min_q = min_q.min(fdr[i]);
        fdr[i] = min_q;
    }
    fdr
}

fn feat_of(e: &LfqEntry) -> [f64; NF] {
    let ppm = if e.observed_mz.is_finite() && e.expected_mz > 0.0 {
        ((e.observed_mz - e.expected_mz) / e.expected_mz * 1e6).abs()
    } else {
        f64::NAN
    };
    let rt = if e.apex_rt.is_finite() && e.expected_rt.is_finite() {
        (e.apex_rt - e.expected_rt).abs()
    } else {
        f64::NAN
    };
    // |IM delta| — NaN on Orbitrap (no IM), where it is imputed to the column
    // median and standardises to a constant, so it contributes nothing there.
    let im = if e.observed_im.is_finite() && e.expected_im != 0.0 {
        (e.observed_im - e.expected_im).abs()
    } else {
        f64::NAN
    };
    [
        ppm,
        rt,
        e.spectral_bhattacharyya as f64,
        e.coelution as f64,
        im,
    ]
}

/// Compute per-(feature, run) q-values for target cells via the semi-supervised
/// QDA rescorer. Every target cell with signal (detected or MBR) competes against
/// the decoys and gets a learned q; a target cell with no signal → q = 1.
pub fn compute_qvalues_qda(entries: &[LfqEntry]) -> HashMap<(usize, usize), f64> {
    struct Raw {
        feat: [f64; NF],
        is_decoy: bool,
        feature_idx: usize,
        run_idx: usize,
        hybrid: f64,
    }

    let mut q_out: HashMap<(usize, usize), f64> = HashMap::new();
    let mut raw: Vec<Raw> = Vec::new();
    for e in entries {
        if e.intensity > 0.0 {
            // Every cell with signal — detected target, MBR target, or decoy —
            // competes; detected cells are NOT given a free pass, so a
            // background-contaminated detected cell in a depleted well is gatable.
            raw.push(Raw {
                feat: feat_of(e),
                is_decoy: e.is_decoy,
                feature_idx: e.feature_idx,
                run_idx: e.run_idx,
                hybrid: e.hybrid_score as f64,
            });
        } else if !e.is_decoy {
            // Target cell with no signal — nothing to keep.
            q_out.insert((e.feature_idx, e.run_idx), 1.0);
        }
    }

    let n = raw.len();
    let n_dec = raw.iter().filter(|r| r.is_decoy).count();
    if n == 0 || n_dec == 0 || n_dec == n {
        // No competition possible — fall back to hybrid ranking for MBR cells.
        let scores: Vec<f64> = raw.iter().map(|r| r.hybrid).collect();
        let dec: Vec<bool> = raw.iter().map(|r| r.is_decoy).collect();
        let q = qvalues(&scores, &dec);
        for (i, r) in raw.iter().enumerate() {
            if !r.is_decoy {
                q_out.insert((r.feature_idx, r.run_idx), q[i]);
            }
        }
        return q_out;
    }

    // Impute non-finite features with column medians, then standardise.
    let mut medians = [0.0f64; NF];
    for k in 0..NF {
        let mut vals: Vec<f64> = raw
            .iter()
            .map(|r| r.feat[k])
            .filter(|v| v.is_finite())
            .collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        medians[k] = if vals.is_empty() { 0.0 } else { vals[vals.len() / 2] };
    }
    let mut feat: Vec<[f64; NF]> = Vec::with_capacity(n);
    for r in &raw {
        let mut f = r.feat;
        for k in 0..NF {
            if !f[k].is_finite() {
                f[k] = medians[k];
            }
        }
        feat.push(f);
    }
    let mut mean = [0.0f64; NF];
    for f in &feat {
        for k in 0..NF {
            mean[k] += f[k];
        }
    }
    for k in 0..NF {
        mean[k] /= n as f64;
    }
    let mut sdv = [0.0f64; NF];
    for f in &feat {
        for k in 0..NF {
            let d = f[k] - mean[k];
            sdv[k] += d * d;
        }
    }
    for k in 0..NF {
        sdv[k] = (sdv[k] / n as f64).sqrt();
        if sdv[k] < 1e-12 {
            sdv[k] = 1.0;
        }
    }
    for f in feat.iter_mut() {
        for k in 0..NF {
            f[k] = (f[k] - mean[k]) / sdv[k];
        }
    }

    let dec: Vec<bool> = raw.iter().map(|r| r.is_decoy).collect();
    let hybrid: Vec<f64> = raw.iter().map(|r| r.hybrid).collect();
    let fold: Vec<usize> = raw.iter().map(|r| r.feature_idx % N_FOLDS).collect();

    // Cross-validated semi-supervised QDA. Held-out scores only.
    let mut perc = vec![f64::NAN; n];
    for f in 0..N_FOLDS {
        let tr: Vec<usize> = (0..n).filter(|&i| fold[i] != f).collect();
        let te: Vec<usize> = (0..n).filter(|&i| fold[i] == f).collect();
        if tr.is_empty() || te.is_empty() {
            continue;
        }
        let tr_feat: Vec<[f64; NF]> = tr.iter().map(|&i| feat[i]).collect();
        let tr_dec: Vec<bool> = tr.iter().map(|&i| dec[i]).collect();
        let neg: Vec<usize> = (0..tr.len()).filter(|&j| tr_dec[j]).collect();
        let mut score: Vec<f64> = tr.iter().map(|&i| hybrid[i]).collect();
        let mut model: Option<Qda> = None;

        for _ in 0..N_ITER {
            let q = qvalues(&score, &tr_dec);
            // Confident targets (any cell with a real feature identity or a
            // high-scoring transfer) as positives; decoys as negatives.
            let mut pos: Vec<usize> = (0..tr.len())
                .filter(|&j| !tr_dec[j] && q[j] <= TRAIN_FDR)
                .collect();
            if pos.len() < MIN_POS {
                let mut tgt: Vec<usize> = (0..tr.len()).filter(|&j| !tr_dec[j]).collect();
                tgt.sort_by(|&a, &b| {
                    score[b]
                        .partial_cmp(&score[a])
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then(a.cmp(&b))
                });
                pos = tgt.into_iter().take(MIN_POS).collect();
            }
            let m = match fit_qda(&tr_feat, &pos, &neg) {
                Some(m) => m,
                None => break,
            };
            for j in 0..tr.len() {
                score[j] = m.score(&tr_feat[j]);
            }
            model = Some(m);
        }

        match &model {
            Some(m) => {
                for &i in &te {
                    perc[i] = m.score(&feat[i]);
                }
            }
            None => {
                for &i in &te {
                    perc[i] = hybrid[i];
                }
            }
        }
    }
    for i in 0..n {
        if !perc[i].is_finite() {
            perc[i] = hybrid[i];
        }
    }

    let q = qvalues(&perc, &dec);
    for (i, r) in raw.iter().enumerate() {
        if !r.is_decoy {
            q_out.insert((r.feature_idx, r.run_idx), q[i]);
        }
    }
    q_out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(feature_idx: usize, run_idx: usize, is_decoy: bool, ppm: f64, rt: f64, bc: f32) -> LfqEntry {
        LfqEntry {
            feature_idx,
            run_idx,
            intensity: 1.0,
            hybrid_score: bc,
            spectral_bhattacharyya: bc,
            n_isotopes_found: 3,
            rt_score: 1.0,
            int_score: 1.0,
            coelution: bc,
            is_decoy,
            is_mbr: !is_decoy, // targets are MBR here
            expected_rt: 10.0,
            apex_rt: 10.0 + rt,
            peak_width_rt: 0.1,
            expected_mz: 500.0,
            observed_mz: 500.0 * (1.0 + ppm / 1e6),
            expected_im: 0.0,
            observed_im: f64::NAN,
        }
    }

    #[test]
    fn qda_separates_tight_targets_from_diffuse_decoys() {
        // Targets: tight cluster (near-0 ppm/RT, high BC). Decoys: diffuse.
        let mut entries = Vec::new();
        for i in 0..400 {
            let jitter = ((i * 7) % 11) as f64 / 11.0 - 0.5;
            entries.push(entry(i, i % 20, false, 0.2 * jitter, 0.01 * jitter, 0.9));
            // decoys spread across the ppm/RT window with low BC
            let spread = ((i * 13) % 20) as f64 / 20.0;
            let dbc = 0.3f32 + 0.4 * ((i % 3) as f32) / 3.0;
            entries.push(entry(10_000 + i, i % 20, true, 8.0 * spread, 0.3 * spread, dbc));
        }
        let q = compute_qvalues_qda(&entries);
        // Every target cell got a q-value; most should be small (well separated).
        let mut n_small = 0;
        for i in 0..400 {
            let key = (i, i % 20);
            let qi = *q.get(&key).expect("target should have a q-value");
            if qi <= 0.05 {
                n_small += 1;
            }
        }
        assert!(
            n_small > 300,
            "QDA should recover most tight targets at q<=0.05, got {n_small}/400"
        );
    }

    #[test]
    fn qda_is_deterministic() {
        let mut entries = Vec::new();
        for i in 0..300 {
            let jitter = ((i * 7) % 11) as f64 / 11.0 - 0.5;
            entries.push(entry(i, i % 10, false, 0.2 * jitter, 0.01 * jitter, 0.9));
            let spread = ((i * 13) % 20) as f64 / 20.0;
            entries.push(entry(10_000 + i, i % 10, true, 8.0 * spread, 0.3 * spread, 0.4));
        }
        let a = compute_qvalues_qda(&entries);
        let b = compute_qvalues_qda(&entries);
        assert_eq!(a.len(), b.len());
        for (k, v) in &a {
            assert_eq!(b.get(k), Some(v), "QDA rescorer must be deterministic");
        }
    }
}
