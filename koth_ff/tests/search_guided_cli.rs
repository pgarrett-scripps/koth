//! End-to-end search-guided CLI tests, including legacy-mode output parity.
use koth_ff::{
    models::{Feature, Hill, Polarity, ScoredFeature, PROTON_MASS},
    output::{write_features_tsv, write_hills_tsv},
};
use std::{
    path::{Path, PathBuf},
    process::Command,
    sync::atomic::{AtomicU64, Ordering},
};
static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture(PathBuf);
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn hill(mz: f64, rt: f64, scale: f32) -> Hill {
    let profile: Vec<f32> = (0..41)
        .map(|i| {
            let x = (i as f32 - 20.) / 5.;
            scale * 100. * (-0.5 * x * x).exp()
        })
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
        scan_apex: 20,
        scan_end: 40,
        n_scans: 41,
        skipped_scans: 0,
        intensity_sum: profile.iter().map(|&x| x as f64).sum(),
        intensity_max: 100. * scale as f64,
        hill_score: 1.,
        intensity_profile: profile.into(),
        isolation_window: None,
        faims_cv: None,
    }
}
impl Fixture {
    fn new(features: bool) -> Self {
        let root = std::env::temp_dir().join(format!(
            "koth-guided-cli-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        for (run, shift) in [("A", 0.), ("B", 1.)] {
            let p = root.join("batch").join(run);
            std::fs::create_dir_all(&p).unwrap();
            let mut hills: Vec<_> = [1., 0.55, 0.17]
                .iter()
                .enumerate()
                .map(|(i, &scale)| {
                    hill(
                        500. + PROTON_MASS + i as f64 * 1.003354835 / 2.,
                        5. + shift,
                        scale,
                    )
                })
                .collect();
            let anchors: Vec<_> = (0..12)
                .map(|i| {
                    let h = hill(600. + i as f64 * 20., 0.5 + i as f64 * 0.7 + shift, 1.);
                    hills.push(h.clone());
                    ScoredFeature {
                        feature: Feature {
                            hills: vec![h],
                            charge: 2,
                            cosine_score: 0.99,
                            ppm_error: 0.,
                        },
                        neutron_offset: 0,
                        isotope_score: 0.99,
                        cosine_score: 0.99,
                        combined_score: 0.98,
                        theoretical_pattern: vec![0.7, 0.2, 0.1],
                    }
                })
                .collect();
            write_hills_tsv(&hills, &p.join("hills.tsv")).unwrap();
            write_features_tsv(
                if features { &anchors } else { &[] },
                &p.join("features.tsv"),
                Polarity::Positive,
            )
            .unwrap();
        }
        std::fs::write(root.join("targets.tsv"),"modified_peptide\trun\tcharge\tneutral_mass\trt_minutes\tid_qvalue\nPEPTIDE\tA.mzML\t2\t1000\t5\t0.001\n").unwrap();
        std::fs::write(root.join("sage.tsv"),"peptide\tfilename\tcharge\tcalcmass\trt\tspectrum_q\tpeptide_q\trank\tlabel\nPEPTIDE\tA.mzML\t2\t1000\t5\t0.001\t0.001\t1\t1\n").unwrap();
        std::fs::write(root.join("config.toml"),"[alignment]\nreference_run = \"A\"\n[lfq]\nrt_window_pct = 0.05\n[lfq.consensus]\nmax_group_qvalue = 1.0\n[output]\nexport_long = true\n").unwrap();
        Self(root)
    }
    fn run(&self, out: &str, flags: &[&str]) -> PathBuf {
        let result = Command::new(env!("CARGO_BIN_EXE_koth_align"))
            .current_dir(&self.0)
            .args([
                "batch",
                "--output",
                out,
                "--config",
                "config.toml",
                "--log-level",
                "error",
            ])
            .args(flags)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        self.0.join(out)
    }
}
fn rows(path: &Path) -> Vec<csv::StringRecord> {
    csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap()
        .records()
        .map(Result::unwrap)
        .collect()
}
#[test]
fn cli_sage_and_generic_targets_agree_and_no_mbr_disables_transfers() {
    let f = Fixture::new(true);
    let generic = f.run("generic", &["--targets", "targets.tsv"]);
    let sage = f.run("sage", &["--sage-psms", "sage.tsv"]);
    assert_eq!(
        std::fs::read(generic.join("peptide_quant.tsv")).unwrap(),
        std::fs::read(sage.join("peptide_quant.tsv")).unwrap()
    );
    let r = rows(&generic.join("peptide_quant.tsv"));
    assert_eq!(r.len(), 2);
    assert_eq!(&r[0][4], "direct_ms2");
    assert_eq!(&r[1][4], "mbr");
    assert!(r[0][6].parse::<f64>().unwrap() > 0.);
    assert!(r[1][6].parse::<f64>().unwrap() > 0.);
    assert_eq!(&r[1][9], "");
    let direct = f.run("direct", &["--targets", "targets.tsv", "--no-mbr"]);
    let r = rows(&direct.join("peptide_quant.tsv"));
    assert_eq!(&r[1][5], "not_attempted");
    assert_eq!(&r[1][6], "0");
    let manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(generic.join("matrix_manifest.json")).unwrap())
            .unwrap();
    assert_eq!(manifest["schema_version"], 4);
}
#[test]
fn cli_extracts_direct_targets_with_zero_detected_features() {
    let f = Fixture::new(false);
    let out = f.run("out", &["--targets", "targets.tsv"]);
    let r = rows(&out.join("peptide_quant.tsv"));
    assert!(r[0][6].parse::<f64>().unwrap() > 0.);
    assert_eq!(&r[1][5], "not_attempted");
}
#[test]
fn search_defaults_optionally_match_frozen_release_byte_for_byte() {
    let Some(binary) = std::env::var_os("KOTH_LEGACY_ALIGN") else {
        return;
    };
    let f = Fixture::new(true);
    for (name, flags) in [
        ("mbr", vec!["--targets", "targets.tsv"]),
        ("direct", vec!["--targets", "targets.tsv", "--no-mbr"]),
    ] {
        let out = f.run(name, &flags);
        let old = format!("legacy_{name}");
        let result = Command::new(&binary)
            .current_dir(&f.0)
            .args([
                "batch",
                "--output",
                &old,
                "--config",
                "config.toml",
                "--log-level",
                "error",
            ])
            .args(&flags)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        for file in [
            "peptide_quant.tsv",
            "peptide_intensity_matrix.tsv",
            "search_rejections.tsv",
            "lfq_details.tsv",
            "qvalue_matrix.tsv",
        ] {
            assert_eq!(
                std::fs::read(out.join(file)).unwrap(),
                std::fs::read(f.0.join(&old).join(file)).unwrap(),
                "default parity {name}/{file}"
            );
        }
        // Build provenance must differ for an uncommitted development binary;
        // all scientific fields, input hashes and output hashes must agree.
        let read_manifest = |path: PathBuf| {
            let mut value: serde_json::Value =
                serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
            value.as_object_mut().unwrap().remove("koth_version");
            value
        };
        assert_eq!(
            read_manifest(out.join("search_manifest.json")),
            read_manifest(f.0.join(&old).join("search_manifest.json"))
        );
    }
}
#[test]
fn id_free_mode_keeps_existing_schema_and_optionally_matches_legacy_binary() {
    let f = Fixture::new(true);
    let out = f.run("unguided", &[]);
    assert!(!out.join("search_manifest.json").exists());
    let manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(out.join("matrix_manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["schema_version"], 3);
    if let Some(binary) = std::env::var_os("KOTH_LEGACY_ALIGN") {
        let result = Command::new(binary)
            .current_dir(&f.0)
            .args([
                "batch",
                "--output",
                "legacy",
                "--config",
                "config.toml",
                "--log-level",
                "error",
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        for file in [
            "consensus_features.tsv",
            "intensity_matrix.tsv",
            "qvalue_matrix.tsv",
            "decoy_intensity_matrix.tsv",
            "lfq_details.tsv",
            "feature_observations.tsv",
            "lfq_matrix.long.tsv",
        ] {
            assert_eq!(
                std::fs::read(out.join(file)).unwrap(),
                std::fs::read(f.0.join("legacy").join(file)).unwrap(),
                "legacy parity: {file}"
            );
        }
    }
}
#[test]
fn cli_rejects_conflicting_modes_and_invalid_confidence() {
    for flags in [
        vec!["--targets", "x", "--sage-psms", "y"],
        vec!["--no-mbr"],
        vec!["--targets", "x", "--max-id-qvalue", "NaN"],
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_koth_align"))
            .arg("missing-batch")
            .args(flags)
            .output()
            .unwrap();
        assert!(!result.status.success());
    }
}
