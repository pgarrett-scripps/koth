//! Exclusive LFQ extraction. A sample is identified in the immutable input hill,
//! not by an isotope row or an RT bin. Target and decoy worlds never deplete one
//! another. Group membership and its permutation statistic are not pooled here.
use std::collections::{HashMap, HashSet};

use super::{consensus::ConsensusFeature, grid::SampleId, LfqEntry};

pub(super) struct CellCandidate {
    pub entry: LfqEntry,
    pub support: Vec<SampleId>,
    pub context: Vec<SampleId>,
    pub peak_start: f64,
    pub peak_end: f64,
    pub peak_rows: usize,
}

/// Rank native envelope/co-elution evidence, penalized for unexplained
/// preceding-isotope signal, consistently across runs. Neither
/// target labels, group q-values nor cell q-values participate in this decision.
/// Sorting before summation makes the preference independent of run iteration.
pub(super) fn global_priority(
    entries: &[LfqEntry],
    groups: &[ConsensusFeature],
) -> [Vec<usize>; 2] {
    std::array::from_fn(|side| {
        let mut evidence = vec![Vec::new(); groups.len()];
        for e in entries
            .iter()
            .filter(|e| usize::from(e.is_decoy) == side && e.intensity > 0.0)
        {
            let fit = f64::from(e.spectral_bhattacharyya)
                * f64::from(e.coelution)
                * (1.0 - f64::from(e.preceding_signal_fraction));
            if fit.is_finite() {
                evidence[e.feature_idx].push(fit);
            }
        }
        let scores: Vec<f64> = evidence
            .iter_mut()
            .map(|values| {
                values.sort_by(f64::total_cmp);
                values.iter().sum::<f64>() / values.len().max(1) as f64
            })
            .collect();
        let mut order: Vec<_> = (0..groups.len()).collect();
        order.sort_by(|&a, &b| {
            scores[b]
                .total_cmp(&scores[a])
                .then_with(|| evidence[b].len().cmp(&evidence[a].len()))
                .then_with(|| groups[a].charge.cmp(&groups[b].charge))
                .then_with(|| groups[a].ref_mz.total_cmp(&groups[b].ref_mz))
                .then_with(|| groups[a].ref_rt.total_cmp(&groups[b].ref_rt))
                .then_with(|| groups[a].ref_im.total_cmp(&groups[b].ref_im))
                .then_with(|| a.cmp(&b))
        });
        let mut rank = vec![0; groups.len()];
        for (r, i) in order.into_iter().enumerate() {
            rank[i] = r;
        }
        rank
    })
}

/// Partition a hill only at a valley below half of both neighboring peak
/// heights. A short triangular smoother prevents individual noisy scans from
/// creating aliases. This affects ownership boundaries only, never intensity.
/// A unimodal hill is one segment; well-resolved peaks on a long hill remain
/// separately available. The half-height convention matches peak expansion's
/// existing internal half-score rule; it is not a new user setting.
fn peak_segments(profile: &[f32]) -> Vec<(usize, usize)> {
    let n = profile.len();
    if n == 0 {
        return Vec::new();
    }
    let value = |i: usize| {
        let v = profile[i];
        if v.is_finite() {
            v.max(0.0)
        } else {
            0.0
        }
    };
    let smooth: Vec<_> = (0..n)
        .map(|i| (value(i.saturating_sub(1)) + 2.0 * value(i) + value((i + 1).min(n - 1))) / 4.0)
        .collect();
    let peaks: Vec<_> = (0..n)
        .filter(|&i| {
            smooth[i] > 0.0
                && (i == 0 || smooth[i] >= smooth[i - 1])
                && (i + 1 == n || smooth[i] > smooth[i + 1])
        })
        .collect();
    let mut bounds = Vec::new();
    let mut start = 0;
    if let Some(&first) = peaks.first() {
        let mut previous = first;
        for &next in peaks.iter().skip(1) {
            let valley = (previous..=next)
                .min_by(|&a, &b| smooth[a].total_cmp(&smooth[b]))
                .unwrap();
            if smooth[valley] * 2.0 <= smooth[previous].min(smooth[next]) {
                let split = valley + 1;
                bounds.push((start, split));
                start = split;
                previous = next;
            } else if smooth[next] > smooth[previous] {
                previous = next;
            }
        }
    }
    bounds.push((start, n));
    bounds
}

fn clear(c: &mut CellCandidate, status: &'static str) {
    c.support.clear();
    c.entry.intensity = 0.0;
    c.entry.hybrid_score = 0.0;
    c.entry.spectral_bhattacharyya = 0.0;
    c.entry.n_isotopes_found = 0;
    c.entry.rt_score = 0.0;
    c.entry.int_score = 0.0;
    c.entry.coelution = 0.0;
    c.entry.apex_rt = f64::NAN;
    c.entry.peak_width_rt = 0.0;
    c.entry.observed_mz = f64::NAN;
    c.entry.observed_im = f64::NAN;
    c.entry.ownership_status = status;
}

/// Competing candidates are evaluated before any claims. After exclusions the
/// grid, isotope scores and integration are rebuilt from the remaining samples.
/// A majority-overlapped original peak is not rediscovered as a truncated tail;
/// a resolved second peak may survive. The majority rule is an internal
/// ambiguity safeguard, not an admission threshold on detector quality.
pub(super) fn resolve_run(
    mut candidates: Vec<CellCandidate>,
    priority: &[Vec<usize>; 2],
    hills: &[crate::models::Hill],
    extract: impl Fn(usize, bool, &HashSet<SampleId>) -> CellCandidate,
    contest_on_peak_overlap: bool,
) -> Vec<LfqEntry> {
    candidates.sort_by_key(|c| {
        (
            c.entry.is_decoy,
            priority[usize::from(c.entry.is_decoy)][c.entry.feature_idx],
        )
    });
    let mut segments = HashMap::<u32, Vec<(usize, usize)>>::new();
    let mut claimed = HashSet::new();
    let mut owner = HashMap::new();
    let mut current_side = false;
    let mut out = Vec::with_capacity(candidates.len());
    for mut c in candidates {
        if c.entry.is_decoy != current_side {
            claimed.clear();
            owner.clear();
            current_side = c.entry.is_decoy;
        }
        let overlap = c.support.iter().filter(|s| claimed.contains(*s)).count();
        let context_overlap = c.context.iter().filter(|s| claimed.contains(*s)).count();
        // A competitor elsewhere in the window does not make this cell's own
        // peak ambiguous, and re-extracting against it manufactures a residual
        // where a clean measurement already existed.
        let contested = if contest_on_peak_overlap {
            overlap > 0
        } else {
            context_overlap > 0
        };
        if contested {
            let competitor = c
                .context
                .iter()
                .filter_map(|s| owner.get(s).copied())
                .min_by_key(|&i| priority[usize::from(current_side)][i]);
            let old_start = c.peak_start;
            let old_end = c.peak_end;
            let old_apex_error = (c.entry.apex_rt - c.entry.expected_rt).abs();
            let majority = overlap * 2 >= c.support.len();
            let old_excluded = context_overlap;
            c = extract(c.entry.feature_idx, c.entry.is_decoy, &claimed);
            c.entry.excluded_samples = old_excluded;
            c.entry.competing_feature = competitor;
            c.entry.ownership_status = "residual";
            if c.entry.intensity <= 0.0 {
                clear(&mut c, "shared_signal_only");
            } else if c.peak_rows < 2
                || (majority
                    && ((c.peak_start < old_end && c.peak_end > old_start)
                        || (c.entry.apex_rt - c.entry.expected_rt).abs() >= old_apex_error))
            {
                clear(&mut c, "ambiguous_residual");
            }
        }
        c.entry.owned_samples = c.support.len();
        if c.entry.intensity > 0.0 {
            for &s in &c.support {
                assert!(
                    !claimed.contains(&s),
                    "LFQ final cells reused an original signal sample"
                );
            }
            // Reserve the native chromatographic basin containing each used
            // sample. This prevents a second group from taking the next scan of
            // the same peak merely because its grid bins differ.
            let mut reservations = Vec::new();
            for &(hill, sample) in &c.support {
                let profile = &hills[hill as usize].intensity_profile;
                let bounds = segments
                    .entry(hill)
                    .or_insert_with(|| peak_segments(profile));
                let &(lo, hi) = bounds
                    .iter()
                    .find(|&&(lo, hi)| sample >= lo && sample < hi)
                    .expect("sample belongs to a native peak segment");
                reservations.push((hill, lo, hi));
            }
            reservations.sort_unstable();
            reservations.dedup();
            for (hill, lo, hi) in reservations {
                let profile = &hills[hill as usize].intensity_profile;
                for (i, &value) in profile.iter().enumerate().take(hi).skip(lo) {
                    if value > 0.0 && value.is_finite() {
                        claimed.insert((hill, i));
                        owner.insert((hill, i), c.entry.feature_idx);
                    }
                }
            }
        }
        out.push(c.entry);
    }
    out.sort_by_key(|e| (e.feature_idx, e.is_decoy));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lfq::{
        grid::{SortedHills, XicGrid, C13_NEUTRON},
        quantify_cell,
        score::ColumnScores,
        LfqConfig, PhaseTimers,
    };
    use crate::models::Hill;
    use std::sync::Arc;
    fn hill(
        mz: f64,
        rt_start: f64,
        rt_end: f64,
        im: f64,
        scan_start: usize,
        profile: Vec<f32>,
    ) -> Hill {
        let n = profile.len();
        let int_sum: f64 = profile.iter().map(|&x| x as f64).sum();
        let int_max: f64 = profile.iter().copied().fold(0.0f32, f32::max) as f64;
        Hill {
            hill_id: 0,
            mz,
            mz_std: 0.0,
            mz_se: 0.0,
            rt: (rt_start + rt_end) / 2.0,
            rt_start,
            rt_end,
            rt_width: rt_end - rt_start,
            im,
            im_std: 0.0,
            scan_start,
            scan_apex: scan_start + n / 2,
            scan_end: scan_start + n.saturating_sub(1),
            n_scans: n,
            skipped_scans: 0,
            intensity_sum: int_sum,
            intensity_max: int_max,
            hill_score: 1.0,
            intensity_profile: Arc::from(profile.as_slice()),
            isolation_window: None,
            faims_cv: None,
        }
    }

    fn group(mz: f64, rt: f64, charge: u8) -> ConsensusFeature {
        ConsensusFeature {
            ref_mz: mz,
            ref_rt: rt,
            ref_im: 0.0,
            charge,
            theoretical_pattern: vec![0.6, 0.3, 0.1],
            neutral_mass: (mz - 1.007276466621) * f64::from(charge),
            seed_combined_score: 0.2,
            seed_run_idx: 0,
            n_contributing_runs: 2,
            per_run_feature: vec![Some(0); 2],
            members: vec![(0, 0), (1, 0)],
            seed_feature_idx: 0,
            group_score: 1.0,
            group_qvalue: 0.01,
        }
    }

    fn envelope(mz: f64, charge: u8, peaks: &[(f64, f32)]) -> Vec<Hill> {
        [0.6, 0.3, 0.1]
            .iter()
            .enumerate()
            .map(|(iso, scale)| {
                let profile = (0..100)
                    .map(|i| {
                        let rt = (i as f64 + 0.5) / 100.0;
                        peaks
                            .iter()
                            .map(|&(centre, height)| {
                                if (rt - centre).abs() > 0.06 {
                                    0.0
                                } else {
                                    height * (-((rt - centre) / 0.02).powi(2) / 2.0).exp() as f32
                                }
                            })
                            .sum::<f32>()
                            * scale
                    })
                    .collect();
                hill(
                    mz + iso as f64 * C13_NEUTRON / f64::from(charge),
                    0.005,
                    0.995,
                    0.0,
                    0,
                    profile,
                )
            })
            .collect()
    }

    fn extract(
        hills: &[Hill],
        g: &ConsensusFeature,
        i: usize,
        decoy: bool,
        excluded: &HashSet<SampleId>,
    ) -> CellCandidate {
        let config = LfqConfig::default();
        let mut grid = XicGrid::empty(config.n_isotopes, config.grid_cols, 0.0, 1.0);
        let mut scores = ColumnScores::new(config.grid_cols);
        let mut totals = vec![0.0; config.grid_cols];
        let mut obs = vec![0.0; config.n_isotopes];
        let template = config.isotope_model.model().distribution(g.neutral_mass);
        quantify_cell(
            &mut grid,
            &mut scores,
            &mut totals,
            &mut obs,
            hills,
            &SortedHills::from_hills(hills),
            &[],
            &g.theoretical_pattern,
            &template,
            &config,
            &PhaseTimers::new(),
            i,
            0,
            g.charge,
            g.ref_mz,
            g.ref_rt,
            0.0,
            0.5,
            None,
            decoy,
            true,
            excluded,
        )
    }

    fn resolve(hills: &[Hill], groups: &[ConsensusFeature], reverse: bool) -> Vec<LfqEntry> {
        let empty = HashSet::new();
        let mut candidates: Vec<_> = groups
            .iter()
            .enumerate()
            .map(|(i, g)| extract(hills, g, i, false, &empty))
            .collect();
        let evidence: Vec<_> = candidates.iter().map(|c| c.entry.clone()).collect();
        let priority = global_priority(&evidence, groups);
        if reverse {
            candidates.reverse();
        }
        resolve_run(
            candidates,
            &priority,
            hills,
            |i, d, x| extract(hills, &groups[i], i, d, x),
            false,
        )
    }

    #[test]
    fn duplicate_grids_report_one_peak_and_preserve_raw_input() {
        let hills = envelope(500.0, 2, &[(0.5, 1000.0)]);
        let before: Vec<_> = hills.iter().map(|h| h.intensity_profile.to_vec()).collect();
        let groups = [group(500.0, 0.5, 2), group(500.00001, 0.5, 2)];
        let out = resolve(&hills, &groups, false);
        assert_eq!(out.iter().filter(|e| e.intensity > 0.0).count(), 1);
        assert!(out
            .iter()
            .any(|e| e.excluded_samples > 0 && e.competing_feature.is_some()));
        assert_eq!(
            before,
            hills
                .iter()
                .map(|h| h.intensity_profile.to_vec())
                .collect::<Vec<_>>()
        );
        let reversed = resolve(&hills, &groups, true);
        assert_eq!(
            out.iter().map(|e| e.intensity).collect::<Vec<_>>(),
            reversed.iter().map(|e| e.intensity).collect::<Vec<_>>()
        );
    }

    #[test]
    fn isotope_alias_cannot_report_reused_envelope_as_a_second_peak() {
        let hills = envelope(500.0, 2, &[(0.5, 1000.0)]);
        let groups = [
            group(500.0, 0.5, 2),
            group(500.0 + C13_NEUTRON / 2.0, 0.5, 2),
        ];
        let out = resolve(&hills, &groups, false);
        assert!(out[0].intensity > 0.0);
        assert_eq!(out[1].intensity, 0.0);
    }

    #[test]
    fn resolved_second_peak_survives_after_stronger_shared_peak_is_excluded() {
        let hills = envelope(500.0, 2, &[(0.3, 1000.0), (0.7, 700.0)]);
        let groups = [group(500.0, 0.3, 2), group(500.0, 0.7, 2)];
        // Both broad windows initially select the stronger peak at 0.3.
        // Fix canonical preference to the group whose expectation matches it.
        let empty = HashSet::new();
        let candidates = groups
            .iter()
            .enumerate()
            .map(|(i, g)| extract(&hills, g, i, false, &empty))
            .collect();
        let out = resolve_run(
            candidates,
            &[vec![0, 1], vec![0, 1]],
            &hills,
            |i, d, x| extract(&hills, &groups[i], i, d, x),
            false,
        );
        assert!(out[0].intensity > 0.0 && out[1].intensity > 0.0);
        assert!((out[0].apex_rt - 0.3).abs() < 0.03);
        assert!((out[1].apex_rt - 0.7).abs() < 0.03);
        assert_eq!(out[1].ownership_status, "residual");
    }

    #[test]
    fn peak_overlap_rule_still_detects_a_genuine_collision() {
        // Both groups select the same stronger peak first, so the second cell's
        // own peak really is claimed. The narrower rule must still contest it,
        // otherwise it would be licensing double counting rather than avoiding
        // a needless re-extraction.
        let hills = envelope(500.0, 2, &[(0.3, 1000.0), (0.7, 700.0)]);
        let groups = [group(500.0, 0.3, 2), group(500.0, 0.7, 2)];
        let empty = HashSet::new();
        let build = || {
            groups
                .iter()
                .enumerate()
                .map(|(i, g)| extract(&hills, g, i, false, &empty))
                .collect::<Vec<_>>()
        };
        for peak_rule in [false, true] {
            let out = resolve_run(
                build(),
                &[vec![0, 1], vec![0, 1]],
                &hills,
                |i, d, x| extract(&hills, &groups[i], i, d, x),
                peak_rule,
            );
            assert_eq!(out[1].ownership_status, "residual");
            assert!((out[1].apex_rt - 0.7).abs() < 0.03);
            assert!(out[0].intensity > 0.0 && out[1].intensity > 0.0);
        }
    }

    #[test]
    fn distinct_charge_ions_and_nearby_resolved_masses_survive() {
        let mut hills = envelope(500.0, 2, &[(0.5, 1000.0)]);
        hills.extend(envelope(333.67151774, 3, &[(0.5, 800.0)]));
        hills.extend(envelope(500.05, 2, &[(0.5, 600.0)]));
        let groups = [
            group(500.0, 0.5, 2),
            group(333.67151774, 0.5, 3),
            group(500.05, 0.5, 2),
        ];
        assert!(resolve(&hills, &groups, false)
            .iter()
            .all(|e| e.intensity > 0.0 && e.excluded_samples == 0));
    }

    #[test]
    fn target_and_decoy_worlds_are_symmetric_and_independent() {
        let hills = envelope(500.0, 2, &[(0.5, 1000.0)]);
        let groups = [group(500.0, 0.5, 2), group(500.00001, 0.5, 2)];
        let mut candidates = Vec::new();
        for i in 0..2 {
            for decoy in [false, true] {
                candidates.push(extract(&hills, &groups[i], i, decoy, &HashSet::new()));
            }
        }
        let evidence: Vec<_> = candidates.iter().map(|c| c.entry.clone()).collect();
        let priority = global_priority(&evidence, &groups);
        let out = resolve_run(
            candidates,
            &priority,
            &hills,
            |i, d, x| extract(&hills, &groups[i], i, d, x),
            false,
        );
        let targets: Vec<_> = out
            .iter()
            .filter(|e| !e.is_decoy)
            .map(|e| (e.intensity, e.excluded_samples))
            .collect();
        let decoys: Vec<_> = out
            .iter()
            .filter(|e| e.is_decoy)
            .map(|e| (e.intensity, e.excluded_samples))
            .collect();
        assert_eq!(targets, decoys);
        assert!(targets[0].0 > 0.0);
    }

    #[test]
    fn no_signal_stays_empty() {
        let groups = [group(500.0, 0.5, 2), group(500.00001, 0.5, 2)];
        assert!(resolve(&[], &groups, false)
            .iter()
            .all(|e| e.intensity == 0.0 && e.owned_samples == 0));
    }
    #[test]
    fn nearby_slices_of_one_sparse_native_peak_cannot_be_reused() {
        let mut hills = envelope(500.0, 2, &[(0.5, 1000.0)]);
        // Sparse profiles reproduce the real Bruker case: integration can stop
        // at empty grid columns and initially select just one scan per hill.
        for h in &mut hills {
            let scale = h.intensity_max as f32 / 1000.0;
            h.intensity_profile = Arc::from(vec![
                20.0 * scale,
                80.0 * scale,
                100.0 * scale,
                90.0 * scale,
                40.0 * scale,
            ]);
            h.rt_start = 0.42;
            h.rt_end = 0.58;
        }
        let groups = [group(500.0, 0.48, 2), group(500.0, 0.52, 2)];
        let out = resolve(&hills, &groups, false);
        assert_eq!(out.iter().filter(|e| e.intensity > 0.0).count(), 1);
    }

    #[test]
    fn peak_segments_ignore_scan_noise_but_preserve_deep_valleys() {
        assert_eq!(peak_segments(&[1.0, 0.0, 1.0]), vec![(0, 3)]);
        assert_eq!(
            peak_segments(&[0.0, 5.0, 10.0, 8.0, 9.0, 3.0, 0.0]),
            vec![(0, 7)]
        );
        let segments = peak_segments(&[0.0, 5.0, 10.0, 5.0, 0.0, 0.0, 0.0, 3.0, 7.0, 3.0, 0.0]);
        assert_eq!(segments.len(), 2);
        assert_eq!(segments[0].0, 0);
        assert_eq!(segments[0].1, segments[1].0);
        assert_eq!(segments[1].1, 11);
    }

    #[test]
    fn cross_run_preference_is_stable_when_runs_and_group_input_are_reordered() {
        let hills = envelope(500.0, 2, &[(0.5, 1000.0)]);
        let groups = vec![
            group(500.0, 0.5, 2),
            group(500.0 + C13_NEUTRON / 2.0, 0.5, 2),
        ];
        let mut entries = Vec::new();
        for run in 0..3 {
            for (i, g) in groups.iter().enumerate() {
                let mut e = extract(&hills, g, i, false, &HashSet::new()).entry;
                e.run_idx = run;
                entries.push(e);
            }
        }
        let expected = global_priority(&entries, &groups);
        entries.reverse();
        assert_eq!(expected, global_priority(&entries, &groups));
        for e in &mut entries {
            e.feature_idx = 1 - e.feature_idx;
        }
        let reversed = global_priority(&entries, &[groups[1].clone(), groups[0].clone()]);
        assert_eq!(expected[0][0], reversed[0][1]);
        assert_eq!(expected[0][1], reversed[0][0]);
    }
    #[test]
    fn higher_mass_wins_when_lower_mono_hypothesis_has_no_m_peak() {
        let hills = envelope(500.0, 2, &[(0.5, 1000.0)]);
        let groups = [
            group(500.0 - C13_NEUTRON / 2.0, 0.5, 2),
            group(500.0, 0.5, 2),
        ];
        let out = resolve(&hills, &groups, false);
        assert_eq!(out[0].intensity, 0.0);
        assert!(out[1].intensity > 0.0);
    }

    #[test]
    fn preceding_peak_penalty_requires_coelution() {
        let base = group(500.0 + C13_NEUTRON / 2.0, 0.5, 2);
        let aligned = envelope(500.0, 2, &[(0.5, 1000.0)]);
        let mut separated = aligned.clone();
        separated[0] = envelope(500.0, 2, &[(0.2, 1000.0)]).remove(0);
        let a = extract(&aligned, &base, 0, false, &HashSet::new());
        let b = extract(&separated, &base, 0, false, &HashSet::new());
        assert!(a.entry.preceding_signal_fraction > 0.3);
        assert_eq!(b.entry.preceding_signal_fraction, 0.0);
    }
}
