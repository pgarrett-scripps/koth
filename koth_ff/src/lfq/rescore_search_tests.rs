use super::*;

fn entry(i: usize, decoy: bool) -> LfqEntry {
    LfqEntry {
        feature_idx: i,
        run_idx: 0,
        is_decoy: decoy,
        is_mbr: true,
        intensity: 1.0,
        hybrid_score: if decoy { 0.2 } else { 0.9 },
        spectral_bhattacharyya: if decoy { 0.2 } else { 0.9 },
        coelution: if decoy { 0.2 } else { 0.9 },
        n_isotopes_found: 3,
        rt_score: 1.0,
        int_score: 1.0,
        expected_rt: 10.0,
        apex_rt: if decoy { 10.3 } else { 10.01 },
        peak_width_rt: 0.1,
        expected_mz: 500.0,
        observed_mz: if decoy { 500.003 } else { 500.0001 },
        expected_im: 0.0,
        observed_im: f64::NAN,
        owned_samples: 0,
        excluded_samples: 0,
        competing_feature: None,
        ownership_status: "exclusive",
        preceding_signal_fraction: 0.0,
    }
}

#[test]
fn tied_null_never_accepts_by_input_order() {
    let q = conservative_qvalues(&[1.0; 400], &(0..400).map(|i| i >= 200).collect::<Vec<_>>());
    assert!(q.iter().all(|&v| v == 1.0));
    assert_eq!(conservative_qvalues(&[1.0], &[false]), vec![1.0]);
}

#[test]
fn preprocessing_cannot_see_heldout_values() {
    let mut raw = vec![[1.0; NF], [3.0; NF], [f64::NAN; NF]];
    let first = Transform::fit(&raw, &[0, 1]);
    raw[2] = [1e10; NF];
    let second = Transform::fit(&raw, &[0, 1]);
    assert_eq!(first.mean, second.mean);
    assert_eq!(first.median, second.median);
    assert_eq!(first.sd, second.sd);
}

#[test]
fn sparse_charge_falls_back_and_legacy_is_exact() {
    let entries: Vec<_> = (0..150)
        .flat_map(|i| [entry(i, false), entry(i, true)])
        .collect();
    let peptides: Vec<_> = (0..150).map(|i| format!("PEPTIDE{i}")).collect();
    let mut charges = vec![2; 150];
    charges[0] = 4;
    assert_eq!(
        compute(&entries, &peptides, &charges, SearchScoring::Legacy),
        compute_qvalues_qda(&entries)
    );
    assert_eq!(
        compute(&entries, &peptides, &charges, SearchScoring::PeptideGrouped),
        compute(
            &entries,
            &peptides,
            &charges,
            SearchScoring::ChargeStratified
        )
    );
}

#[test]
fn separated_charges_recover_signal_and_zero_signal_stays_rejected() {
    let mut entries: Vec<_> = (0..900)
        .flat_map(|i| [entry(i, false), entry(i, true)])
        .collect();
    let peptides: Vec<_> = (0..900).map(|i| format!("PEPTIDE{}", i / 3)).collect();
    let charges: Vec<_> = (0..900).map(|i| 2 + (i % 3) as u8).collect();
    entries[0].intensity = 0.0;
    let q = compute(
        &entries,
        &peptides,
        &charges,
        SearchScoring::ChargeStratified,
    );
    assert_eq!(q[&(0, 0)], 1.0);
    assert!(q.values().filter(|&&v| v <= 0.01).count() > 800);
    assert_eq!(
        q,
        compute(
            &entries,
            &peptides,
            &charges,
            SearchScoring::ChargeStratified
        )
    );
    // The fold API deliberately has no charge/run/decoy argument.
    assert_eq!(peptide_fold(&peptides[0]), peptide_fold(&peptides[2]));
}

#[test]
fn signed_features_keep_residual_direction_and_quality_inputs() {
    let mut e = entry(0, false);
    e.observed_mz = e.expected_mz - 0.001;
    e.apex_rt = e.expected_rt - 0.2;
    e.n_isotopes_found = 3;
    e.peak_width_rt = 0.5;
    e.preceding_signal_fraction = 0.25;
    let old = search_features::<5>(&e, SearchScoring::ChargeStratified);
    let signed = search_features::<8>(&e, SearchScoring::ChargeSignedQuality);
    assert!(old[0] > 0.0 && old[1] > 0.0);
    assert_eq!(signed[0], -old[0]);
    assert_eq!(signed[1], -old[1]);
    assert_eq!(signed[5], 4.0_f64.ln());
    assert_eq!(signed[6], 1.5_f64.ln());
    assert_eq!(signed[7], 0.25);
    e.observed_mz = f64::NAN;
    e.peak_width_rt = f64::NAN;
    let missing = search_features::<8>(&e, SearchScoring::ChargeSignedQuality);
    assert!(missing[0].is_nan() && missing[6].is_nan());
}

#[test]
fn signed_modes_preserve_grouping_and_reject_zero_signal() {
    let mut entries: Vec<_> = (0..900)
        .flat_map(|i| [entry(i, false), entry(i, true)])
        .collect();
    let peptides: Vec<_> = (0..900).map(|i| format!("PEPTIDE{}", i / 3)).collect();
    let charges: Vec<_> = (0..900).map(|i| 2 + (i % 3) as u8).collect();
    entries[0].intensity = 0.0;
    for mode in [
        SearchScoring::ChargeSigned,
        SearchScoring::ChargeSignedQuality,
    ] {
        let q = compute(&entries, &peptides, &charges, mode);
        assert_eq!(q[&(0, 0)], 1.0);
        assert!(q.values().filter(|&&v| v <= 0.01).count() > 800);
        assert_eq!(q, compute(&entries, &peptides, &charges, mode));
    }
}
