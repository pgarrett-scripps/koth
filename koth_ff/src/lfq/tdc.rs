use std::collections::HashMap;

use super::LfqEntry;

/// Compute q-values for all target entries via target-decoy competition.
///
/// All entries (target + decoy) are ranked by hybrid_score descending.
/// The running FDR = n_decoy / n_target is computed at each target entry.
/// Q-values are then monotonised by replacing each with the minimum FDR at
/// or below its rank (backward pass).
///
/// Returns a map of `(feature_idx, run_idx) → q_value` for target entries only.
pub fn compute_qvalues(entries: &[LfqEntry]) -> HashMap<(usize, usize), f64> {
    // Sort all entries by score descending
    let mut sorted: Vec<&LfqEntry> = entries.iter().collect();
    sorted.sort_by(|a, b| {
        b.hybrid_score
            .partial_cmp(&a.hybrid_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Forward pass: assign running FDR to each target entry
    let mut n_target = 0usize;
    let mut n_decoy = 0usize;
    // (feature_idx, run_idx, raw_fdr)
    let mut raw: Vec<(usize, usize, f64)> = Vec::new();

    for entry in &sorted {
        if entry.is_decoy {
            n_decoy += 1;
        } else {
            n_target += 1;
            let fdr = if n_target > 0 {
                n_decoy as f64 / n_target as f64
            } else {
                1.0
            };
            raw.push((entry.feature_idx, entry.run_idx, fdr));
        }
    }

    // Backward pass: monotonise (each q-value = min FDR at or below rank)
    let mut min_q = 1.0f64;
    let mut q_values: HashMap<(usize, usize), f64> = HashMap::new();
    for (feat_idx, run_idx, fdr) in raw.iter().rev() {
        min_q = min_q.min(*fdr);
        q_values.insert((*feat_idx, *run_idx), min_q);
    }

    q_values
}
