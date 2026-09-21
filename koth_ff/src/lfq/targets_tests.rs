use super::*;
use crate::{
    alignment::{align_runs, AlignmentConfig},
    models::Hill,
};
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};

static NEXT: AtomicU64 = AtomicU64::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "koth-targets-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&p).unwrap();
        Self(p)
    }
    fn input(&self, text: &str) -> PathBuf {
        let p = self.0.join("input.tsv");
        std::fs::write(&p, text).unwrap();
        p
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn hill(mz: f64, rt: f64, scale: f32) -> Hill {
    let profile: Vec<f32> = [0., 10., 40., 80., 100., 80., 40., 10., 0.]
        .iter()
        .map(|v| v * scale)
        .collect();
    Hill {
        hill_id: 0,
        mz,
        mz_std: 0.,
        mz_se: 0.,
        rt,
        rt_start: rt - 0.2,
        rt_end: rt + 0.2,
        rt_width: 0.4,
        im: 0.,
        im_std: 0.,
        scan_start: 0,
        scan_apex: 4,
        scan_end: 8,
        n_scans: 9,
        skipped_scans: 0,
        intensity_sum: profile.iter().map(|&v| v as f64).sum(),
        intensity_max: 100. * scale as f64,
        hill_score: 1.,
        intensity_profile: profile.into(),
        isolation_window: None,
        faims_cv: None,
    }
}
fn runs() -> Vec<RunInput> {
    ["A", "B"]
        .into_iter()
        .map(|name| {
            let mut h = hill(100., 5., 0.);
            h.rt_start = 0.;
            h.rt_end = 10.;
            RunInput {
                name: name.into(),
                features: vec![],
                hills: vec![h],
                scan_times: vec![],
                rt_bounds: None,
            }
        })
        .collect()
}
fn alignment(runs: &[RunInput], reliable: bool) -> AlignmentResult {
    let mut c = AlignmentConfig::default();
    c.reference_run = Some("A".into());
    let mut a = align_runs(runs, &c);
    if reliable {
        a.alignments.get_mut("B").unwrap().rt_active = vec![true; 10];
    }
    a
}
fn id(run: usize) -> Identification {
    Identification {
        modified_peptide: "PEPTIDE".into(),
        run_idx: run,
        charge: 2,
        neutral_mass: 1000.,
        rt_minutes: 5.,
        im: 0.,
        id_qvalue: 0.001,
        source_row: 1,
    }
}
fn summary() -> ImportSummary {
    ImportSummary {
        path: "test".into(),
        sha256: String::new(),
        format: TargetFormat::Generic,
        max_id_qvalue: 0.01,
        ignore_target_im: false,
        rows_read: 1,
        rows_filtered: 0,
        rows_retained: 1,
    }
}
fn config() -> LfqConfig {
    let mut c = LfqConfig::default();
    c.rt_window_pct = 0.05;
    c.grid_cols = 40;
    c
}
const HEADER: &str = "modified_peptide\trun\tcharge\tneutral_mass\trt_minutes\tid_qvalue\n";

#[test]
fn inferred_charge_preserves_real_psm_charge_and_no_mbr_semantics() {
    let mut r = runs();
    for run in &mut r {
        run.hills.push(hill(800.0, 5.0, 1.0));
    }
    let a = alignment(&r, true);
    let mut c = config();
    c.search_expand_charges = true;
    let (cf, g) = build_targets(vec![id(0)], summary(), &r, &a, &c, 10, false).unwrap();
    assert_eq!(g.targets.len(), 3);
    let i = g.targets.iter().position(|t| t.charge == 3).unwrap();
    assert!(!g.is_direct(i, 0));
    assert!(g.is_inferred(i, 0));
    assert!(g.attempted(i, 0));
    assert!(!g.attempted(i, 1));
    assert_eq!(g.native_anchor(i, 0).unwrap().charge, 2);
    assert_eq!(
        g.targets[i].inferred_charge.as_ref().unwrap().seed.charge,
        2
    );
    assert!((cf[i].ref_mz - (1000.0 / 3.0 + PROTON_MASS)).abs() < 1e-10);
    let matrix = super::super::quantify_guided(&r, &a, &c, cf, g, |_| vec![]);
    let tmp = Scratch::new();
    write_search_outputs(&matrix, &tmp.0, 0.01, true, false).unwrap();
    let text = std::fs::read_to_string(tmp.0.join("peptide_quant.tsv")).unwrap();
    assert!(text.contains("inferred_charge"));
}

#[test]
fn charge_expansion_never_recreates_a_rejected_charge() {
    let mut r = runs();
    for run in &mut r {
        run.hills.push(hill(800.0, 5.0, 1.0));
    }
    let a = alignment(&r, true);
    let mut c = config();
    c.search_expand_charges = true;
    let mut bad = id(0);
    bad.charge = 3;
    let mut conflict = bad.clone();
    conflict.rt_minutes = 8.0;
    let (_, g) =
        build_targets(vec![id(0), bad, conflict], summary(), &r, &a, &c, 10, true).unwrap();
    assert_eq!(g.targets.len(), 1);
    assert_eq!(g.targets[0].charge, 2);
    assert!(g.targets[0].inferred_charge.is_none());
}

#[test]
fn inferred_charge_cannot_quantify_a_lone_isotope() {
    let mut r = runs();
    for run in &mut r {
        run.hills.push(hill(800.0, 5.0, 1.0));
    }
    let a = alignment(&r, true);
    let mut c = config();
    c.search_expand_charges = true;
    let (cf, g) = build_targets(vec![id(0)], summary(), &r, &a, &c, 10, false).unwrap();
    let i = g.targets.iter().position(|t| t.charge == 3).unwrap();
    let matrix = super::super::quantify_guided(&r, &a, &c, cf, g, |_| {
        vec![hill(1000.0 / 3.0 + PROTON_MASS, 5.0, 10.0)]
    });
    assert!(matrix
        .entries
        .iter()
        .any(|e| e.feature_idx == i && !e.is_decoy));
    assert_eq!(matrix.intensity(i, 0), 0.0);
}

#[test]
fn rt_rescue_retains_direct_cells_but_blocks_ambiguous_transfers() {
    let r = runs();
    let a = alignment(&r, true);
    let mut b = id(1);
    b.rt_minutes = 8.0;
    let mut c = config();
    let (_, legacy) =
        build_targets(vec![id(0), b.clone()], summary(), &r, &a, &c, 10, true).unwrap();
    assert!(legacy.targets.is_empty());
    c.search_rt_rescue = true;
    let (_, rescued) = build_targets(vec![id(0), b], summary(), &r, &a, &c, 10, true).unwrap();
    assert_eq!(rescued.targets.len(), 1);
    assert!(rescued.is_direct(0, 0) && rescued.is_direct(0, 1));
    assert!(!rescued.targets[0].transfer_eligible);
    assert!(
        rescued.targets[0]
            .rt_rescue
            .as_ref()
            .unwrap()
            .cross_run_rt_ambiguous
    );
}

#[test]
fn bounded_cluster_rescue_records_outliers_and_rejects_ties() {
    let r = runs();
    let mut close = id(0);
    close.rt_minutes = 5.05;
    close.source_row = 2;
    let mut far = id(0);
    far.rt_minutes = 8.0;
    far.source_row = 3;
    let ranges: Vec<(f64, f64)> = r.iter().map(|x| x.rt_range()).collect();
    let (kept, excluded, _) =
        supported_rt_ids(vec![id(0), close, far.clone()], &r, &ranges, 0.02);
    assert_eq!(kept.len(), 2);
    assert_eq!(excluded, vec![3]);
    let (kept, _, _) = supported_rt_ids(vec![id(0), far], &r, &ranges, 0.02);
    assert!(kept.is_empty());
}

#[test]
fn rt_rescue_never_rescues_mass_or_mobility_conflicts() {
    let r = runs();
    let a = alignment(&r, true);
    let mut c = config();
    c.search_rt_rescue = true;
    for mass in [true, false] {
        let mut first = id(0);
        let mut second = id(1);
        second.rt_minutes = 8.0;
        if mass {
            second.neutral_mass += 10.0;
        } else {
            first.im = 1.0;
            second.im = 1.5;
        }
        let (_, g) = build_targets(vec![first, second], summary(), &r, &a, &c, 10, true).unwrap();
        assert!(g.targets.is_empty());
    }
}

#[test]
fn generic_import_resolves_paths_and_preserves_modifications() {
    let tmp = Scratch::new();
    let p=tmp.input(&format!("{HEADER}PEP[+15.9949]TIDE\tC:\\data\\A.mzML.gz\t2\t1000\t5\t0.001\nBAD\tA\t2\t900\t5\t0.1\n"));
    let (ids, s) = read_identifications(&p, TargetFormat::Generic, &runs(), 0.01, false).unwrap();
    assert_eq!(ids.len(), 1);
    assert_eq!(ids[0].run_idx, 0);
    assert_eq!(ids[0].modified_peptide, "PEP[+15.9949]TIDE");
    assert_eq!(s.rows_filtered, 1);
    assert_eq!(s.sha256.len(), 64);
}
#[test]
fn sage_filters_decoys_rank_and_both_confidences_and_uses_calculated_mass() {
    let tmp = Scratch::new();
    let p=tmp.input("peptide\tfilename\tcharge\tcalcmass\texpmass\trt\taligned_rt\tspectrum_q\tpeptide_q\trank\tlabel\nPEPTIDE\tA.raw\t2\t1000\t1001\t5\t0.5\t0.001\t0.002\t1\t1\nDECOY\tA\t2\t1000\t1000\t5\t0.5\t0.001\t0.001\t1\t-1\nRANK\tA\t2\t1000\t1000\t5\t0.5\t0.001\t0.001\t2\t1\nBADPSM\tA\t2\t1000\t1000\t5\t0.5\t0.1\t0.001\t1\t1\nBADPEP\tA\t2\t1000\t1000\t5\t0.5\t0.001\t0.1\t1\t1\n");
    let (ids, s) = read_identifications(&p, TargetFormat::Sage, &runs(), 0.01, false).unwrap();
    assert_eq!(ids.len(), 1);
    assert_eq!(ids[0].neutral_mass, 1000.);
    assert_eq!(ids[0].rt_minutes, 5.);
    assert_eq!(ids[0].id_qvalue, 0.002);
    assert_eq!(s.rows_filtered, 4);
}
#[test]
fn rejects_invalid_inputs_and_ambiguous_run_names() {
    let tmp = Scratch::new();
    for row in [
        "P\tUNKNOWN\t2\t1000\t5\t0.001",
        "P\tA\t0\t1000\t5\t0.001",
        "P\tA\t2\tNaN\t5\t0.001",
        "P\tA\t2\t1000\t-1\t0.001",
        "P\tA\t2\t1000\t5\t-0.1",
    ] {
        let p = tmp.input(&format!("{HEADER}{row}\n"));
        let e = read_identifications(&p, TargetFormat::Generic, &runs(), 0.01, false).unwrap_err();
        assert!(format!("{e:#}").contains("data row 1"));
    }
    let mut r = runs();
    r[1].name = "A.mzML".into();
    let p = tmp.input(&format!("{HEADER}P\tA\t2\t1000\t5\t0.001\n"));
    assert!(read_identifications(&p, TargetFormat::Generic, &r, 0.01, false).is_err());
    assert!(read_identifications(&p, TargetFormat::Generic, &runs(), f64::NAN, false).is_err());
}
#[test]
fn parquet_sage_boolean_decoys_are_supported() {
    use arrow::{
        array::{ArrayRef, BooleanArray, Float64Array, Int32Array, StringArray},
        datatypes::{DataType, Field, Schema},
        record_batch::RecordBatch,
    };
    use std::sync::Arc;
    let tmp = Scratch::new();
    let path = tmp.0.join("results.sage.parquet");
    let schema = Arc::new(Schema::new(vec![
        Field::new("peptide", DataType::Utf8, false),
        Field::new("filename", DataType::Utf8, false),
        Field::new("charge", DataType::Int32, false),
        Field::new("calcmass", DataType::Float64, false),
        Field::new("rt", DataType::Float64, false),
        Field::new("spectrum_q", DataType::Float64, false),
        Field::new("peptide_q", DataType::Float64, false),
        Field::new("rank", DataType::Int32, false),
        Field::new("is_decoy", DataType::Boolean, false),
    ]));
    let arrays: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(vec!["P", "D"])),
        Arc::new(StringArray::from(vec!["A", "A"])),
        Arc::new(Int32Array::from(vec![2, 2])),
        Arc::new(Float64Array::from(vec![1000., 1000.])),
        Arc::new(Float64Array::from(vec![5., 5.])),
        Arc::new(Float64Array::from(vec![0.001, 0.001])),
        Arc::new(Float64Array::from(vec![0.001, 0.001])),
        Arc::new(Int32Array::from(vec![1, 1])),
        Arc::new(BooleanArray::from(vec![false, true])),
    ];
    let batch = RecordBatch::try_new(schema.clone(), arrays).unwrap();
    let mut w =
        parquet::arrow::ArrowWriter::try_new(std::fs::File::create(&path).unwrap(), schema, None)
            .unwrap();
    w.write(&batch).unwrap();
    w.close().unwrap();
    let (ids, s) = read_identifications(&path, TargetFormat::Sage, &runs(), 0.01, false).unwrap();
    assert_eq!(ids.len(), 1);
    assert_eq!(s.rows_filtered, 1);
}
#[test]
fn deduplicates_psms_but_keeps_charge_and_modified_sequence_separate() {
    let r = runs();
    let a = alignment(&r, true);
    let c = config();
    let mut duplicate = id(0);
    duplicate.id_qvalue = 0.005;
    duplicate.source_row = 2;
    let mut modified = id(0);
    modified.modified_peptide = "PEP[+16]TIDE".into();
    modified.neutral_mass += 16.;
    let mut charged = id(0);
    charged.charge = 3;
    let (cf, g) = build_targets(
        vec![duplicate, id(0), id(1), modified, charged],
        summary(),
        &r,
        &a,
        &c,
        10,
        true,
    )
    .unwrap();
    assert_eq!(cf.len(), 3);
    assert_eq!(g.targets[1].donors[0].as_ref().unwrap().id_qvalue, 0.001);
    assert!(cf
        .iter()
        .all(|c| c.members.is_empty() && c.group_qvalue.is_nan()));
    assert_eq!(
        g.targets.iter().filter(|t| t.donors[1].is_some()).count(),
        1
    );
}
#[test]
fn rejects_incompatible_mass_rt_and_mobility_instead_of_merging() {
    let r = runs();
    let a = alignment(&r, true);
    let c = config();
    for which in 0..3 {
        let mut other = id(1);
        let mut first = id(0);
        match which {
            0 => other.neutral_mass += 1.,
            1 => other.rt_minutes += 2.,
            _ => {
                first.im = 1.;
                other.im = 1.2;
            }
        }
        let (cf, g) = build_targets(vec![first, other], summary(), &r, &a, &c, 10, true).unwrap();
        assert!(cf.is_empty());
        assert_eq!(g.rejected_targets.len(), 1);
    }
}
#[test]
fn missing_alignment_blocks_transfers_but_preserves_direct_extraction() {
    let r = runs();
    let a = alignment(&r, false);
    let c = config();
    let (_, g) = build_targets(vec![id(0)], summary(), &r, &a, &c, 10, true).unwrap();
    assert!(g.attempted(0, 0));
    assert!(!g.attempted(0, 1));
    let (_, g) = build_targets(vec![id(1)], summary(), &r, &a, &c, 10, true).unwrap();
    assert!(!g.attempted(0, 0));
    assert!(g.attempted(0, 1));
}
fn signal(rt: f64) -> Vec<Hill> {
    [1., 0.55, 0.17]
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            hill(
                500. + PROTON_MASS + i as f64 * super::super::grid::C13_NEUTRON / 2.,
                rt,
                v,
            )
        })
        .collect()
}
#[test]
fn targeted_and_transferred_signals_work_without_complete_features_and_follow_rt_warp() {
    let r = runs();
    let mut a = alignment(&r, true);
    let c = config();
    a.alignments.get_mut("B").unwrap().rt_warp.knot_y = vec![-0.1, -0.1];
    for mbr in [false, true] {
        let (cf, g) = build_targets(vec![id(0)], summary(), &r, &a, &c, 10, mbr).unwrap();
        let matrix =
            super::super::quantify_guided(&r, &a, &c, cf, g, |run| signal(5. + run as f64));
        assert!(matrix.intensity(0, 0) > 0.);
        assert_eq!(matrix.intensity(0, 1) > 0., mbr);
        let direct = matrix
            .entries
            .iter()
            .find(|e| !e.is_decoy && e.run_idx == 0)
            .unwrap();
        assert!(!direct.is_mbr);
        if mbr {
            let transfer = matrix
                .entries
                .iter()
                .find(|e| !e.is_decoy && e.run_idx == 1)
                .unwrap();
            assert!(transfer.is_mbr);
            assert!((transfer.expected_rt - 6.).abs() < 1e-9);
            assert_eq!(matrix.entries.iter().filter(|e| e.is_decoy).count(), 2);
        } else {
            assert!(matrix.entries.iter().all(|e| e.run_idx == 0));
        }
    }
}
#[test]
fn donor_native_rt_and_missing_mobility_are_preserved_for_paired_decoys() {
    let r = runs();
    let mut a = alignment(&r, true);
    let c = config();
    a.alignments.get_mut("B").unwrap().im_drift.intercept = 0.1;
    let mut second = id(1);
    second.rt_minutes = 5.1;
    let (cf, g) = build_targets(vec![id(0), second], summary(), &r, &a, &c, 10, true).unwrap();
    let matrix = super::super::quantify_guided(&r, &a, &c, cf, g, |_| signal(5.1));
    let direct = matrix
        .entries
        .iter()
        .find(|e| e.run_idx == 1 && !e.is_decoy)
        .unwrap();
    let decoy = matrix
        .entries
        .iter()
        .find(|e| e.run_idx == 1 && e.is_decoy)
        .unwrap();
    assert_eq!(direct.expected_rt, 5.1);
    assert_eq!(direct.expected_im, 0.);
    assert!((decoy.expected_rt - (5.1 - c.decoy_rt_shift_pct * 10.)).abs() < 1e-9);
    assert!(!direct.is_mbr);
}
#[test]
fn separate_confidence_and_missing_signal_are_exported() {
    let r = runs();
    let a = alignment(&r, true);
    let c = config();
    let (cf, g) = build_targets(vec![id(0)], summary(), &r, &a, &c, 10, true).unwrap();
    let mut matrix =
        super::super::quantify_guided(
            &r,
            &a,
            &c,
            cf,
            g,
            |i| if i == 0 { signal(5.) } else { vec![] },
        );
    matrix.q_values[0] = 0.2;
    let tmp = Scratch::new();
    write_search_outputs(&matrix, &tmp.0, 0.01, true, false).unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(tmp.0.join("peptide_quant.tsv"))
        .unwrap();
    let rows: Vec<_> = reader.records().map(Result::unwrap).collect();
    assert_eq!(&rows[0][4], "direct_ms2");
    assert_eq!(&rows[0][5], "rejected");
    assert_eq!(&rows[0][7], "0");
    assert_eq!(&rows[0][8], "0.2");
    assert_eq!(&rows[0][9], "0.001");
    assert_eq!(&rows[1][4], "mbr");
    assert_eq!(&rows[1][5], "no_signal");
    assert_eq!(&rows[1][9], "");
    let manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(tmp.0.join("search_manifest.json")).unwrap())
            .unwrap();
    assert_eq!(
        manifest["search"]["targets"][0]["modified_peptide"],
        "PEPTIDE"
    );
}
#[test]
fn guided_long_bundle_has_no_fabricated_feature_links() {
    let r = runs();
    let a = alignment(&r, true);
    let c = config();
    let tmp = Scratch::new();
    let (cf, g) = build_targets(vec![id(0)], summary(), &r, &a, &c, 10, true).unwrap();
    let matrix = super::super::quantify_guided(&r, &a, &c, cf, g, |_| signal(5.));
    let sources: Vec<_> = (0..2)
        .map(|i| {
            let p = tmp.0.join(format!("{i}.tsv"));
            std::fs::write(&p, "charge\tmz\trtApex\n").unwrap();
            p
        })
        .collect();
    crate::output::long_matrix::write_bundle(
        &matrix,
        &r,
        &a,
        &sources,
        &crate::config::AlignConfig::default(),
        &tmp.0,
    )
    .unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(tmp.0.join("lfq_matrix.long.tsv"))
        .unwrap();
    let header = reader.headers().unwrap().clone();
    let index = |name| header.iter().position(|h| h == name).unwrap();
    for row in reader.records() {
        let row = row.unwrap();
        assert!(row[index("observation_id")].is_empty());
        assert!(row[index("seed_observation_id")].is_empty());
    }
}

#[test]
fn disabling_rescoring_never_labels_unscored_signal_accepted() {
    let r = runs();
    let a = alignment(&r, true);
    let mut c = config();
    c.run_tdc = false;
    let (cf, g) = build_targets(vec![id(0)], summary(), &r, &a, &c, 10, true).unwrap();
    let matrix = super::super::quantify_guided(&r, &a, &c, cf, g, |_| signal(5.));
    let tmp = Scratch::new();
    write_search_outputs(&matrix, &tmp.0, 1.0, false, false).unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(tmp.0.join("peptide_quant.tsv"))
        .unwrap();
    for row in reader.records() {
        let row = row.unwrap();
        assert_eq!(&row[5], "unscored");
        assert!(row[6].parse::<f64>().unwrap() > 0.0);
        assert_eq!(&row[7], "0");
        assert_eq!(&row[8], "");
    }
}

#[test]
fn ignoring_target_im_skips_non_native_or_invalid_coordinates() {
    let tmp = Scratch::new();
    let p = tmp.input("modified_peptide\trun\tcharge\tneutral_mass\trt_minutes\tid_qvalue\tim\nP\tA\t2\t1000\t5\t0.001\tNaN\n");
    assert!(read_identifications(&p, TargetFormat::Generic, &runs(), 0.01, false).is_err());
    let (ids, summary) =
        read_identifications(&p, TargetFormat::Generic, &runs(), 0.01, true).unwrap();
    assert_eq!(ids[0].im, 0.0);
    assert!(summary.ignore_target_im);
}

#[test]
fn out_of_range_donor_is_direct_only_instead_of_transferring_clipped_rt() {
    let r = runs();
    let a = alignment(&r, true);
    let mut donor = id(0);
    donor.rt_minutes = 20.0;
    let (_, g) = build_targets(vec![donor], summary(), &r, &a, &config(), 10, true).unwrap();
    assert!(g.attempted(0, 0));
    assert!(!g.attempted(0, 1));
}

#[test]
fn entirely_ambiguous_input_emits_one_rejection_without_fabricated_targets() {
    let r = runs();
    let a = alignment(&r, true);
    let c = config();
    let mut conflict = id(1);
    conflict.neutral_mass += 1.0;
    let (cf, g) = build_targets(vec![id(0), conflict], summary(), &r, &a, &c, 10, true).unwrap();
    let matrix = super::super::quantify_guided(&r, &a, &c, cf, g, |_| vec![]);
    let tmp = Scratch::new();
    write_search_outputs(&matrix, &tmp.0, 0.01, true, false).unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(tmp.0.join("search_rejections.tsv"))
        .unwrap();
    let records: Vec<_> = reader.records().map(Result::unwrap).collect();
    assert_eq!(records.len(), 1);
    assert_eq!(&records[0][0], "PEPTIDE");
    assert_eq!(&records[0][2], "inconsistent_mass");
    assert_eq!(matrix.n_features, 0);
}
