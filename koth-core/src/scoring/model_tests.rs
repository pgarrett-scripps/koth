//! Unit tests for `model.rs`. Wired in via `#[path]`, so this is a child module
//! of the production module and sees its private items.

use super::*;

use crate::scoring::elements::cache;

/// The shipped default is the peptide averagine — the model the published
/// benchmark ran under. A change here re-scores every feature koth reports.
#[test]
fn default_model_is_peptide() {
    assert_eq!(IsotopeModel::default(), PEPTIDE);
    assert_eq!(IsotopeModelSpec::default().model(), PEPTIDE);
    assert_eq!(PEPTIDE.residue_mass, 111.1254);
}

/// RNA composition at a mass a 10-mer oligonucleotide actually lands on.
/// Derived from the four chain residues (NMP − H₂O), mean C₉.₅ H₁₁.₇₅ N₃.₇₅ O₇
/// per 321.2916 Da.
#[test]
fn rna_counts_at_three_kilodaltons() {
    assert_eq!(RNA.counts(3000.0), (89, 110, 35, 65, 0));
    assert_eq!(DNA.counts(3000.0), (95, 119, 36, 58, 0));
}

/// Every model's distribution is a normalized probability vector.
#[test]
fn distributions_sum_to_one() {
    for model in [PEPTIDE, RNA, DNA] {
        for mass in [500.0, 1500.0, 3000.0, 5000.0, 9000.0] {
            let sum: f64 = model.distribution(mass).iter().sum();
            assert!((sum - 1.0).abs() < 1e-9, "{model:?} at {mass}: sum={sum}");
        }
    }
}

/// The whole point: at one mass the classes predict different envelopes. A
/// peptide packs ~133 carbons into 3 kDa where RNA packs 89, so the peptide
/// envelope sits much further right — M+1/M0 is 1.62 against RNA's 1.13. Score
/// an oligonucleotide against the peptide model and this is the error.
#[test]
fn rna_envelope_differs_from_peptide() {
    let m = 3000.0;
    let pep = PEPTIDE.distribution(m);
    let rna = RNA.distribution(m);
    let pep_ratio = pep[1] / pep[0];
    let rna_ratio = rna[1] / rna[0];
    assert!(
        pep_ratio > rna_ratio + 0.4,
        "peptide M+1/M0 {pep_ratio:.3} should exceed RNA's {rna_ratio:.3} by a wide margin"
    );
    // And the peptide's monoisotope carries less of the envelope.
    assert!(pep[0] < rna[0] - 0.05, "{} vs {}", pep[0], rna[0]);
}

/// A 5 kDa oligonucleotide needs 109 oxygens. The element cache was sized for
/// tryptic peptides (MAX_O = 100) and clamped silently above that, which is a
/// wrong pattern rather than a slow one. Pins the headroom.
#[test]
fn oxygen_rich_models_are_not_clamped() {
    let (c, h, n, o, s) = RNA.counts(5000.0);
    assert_eq!(o, 109, "RNA at 5 kDa should need 109 O");
    let real = cache().distribution(c, h, n, o, s);
    let clamped = cache().distribution(c, h, n, 100, s);
    assert!(
        real != clamped,
        "distribution at {o} O is identical to the 100-O clamp — the cache is still too small"
    );
    // And the largest oligonucleotide the model is plausibly used on.
    assert_eq!(RNA.counts(9000.0).3, 196);
    assert!(
        cache().distribution(266, 329, 105, 196, 0) != cache().distribution(266, 329, 105, 100, 0)
    );
}

/// A model with no sulfur ignores `sulfur_offsets`: scoring is one template, not
/// a max over templates that vary an atom the analyte does not contain.
#[test]
fn sulfur_free_model_ignores_offsets() {
    let obs = [5.0e7, 3.0e7, 1.5e7, 4.0e6, 8.0e5];
    let m = 3000.0;
    assert!(!RNA.has_sulfur());
    let plain = crate::scoring::averagine::bhattacharyya_score(&obs, &RNA.distribution(m));
    let (scored, n_s) =
        crate::scoring::averagine::bhattacharyya_score_best_sulfur(&obs, m, &[-1, 0, 1], &RNA);
    assert!((scored - plain).abs() < 1e-12, "{scored} vs {plain}");
    assert_eq!(n_s, 0);
    // The peptide model, for contrast, does vary it.
    assert!(PEPTIDE.has_sulfur());
}

#[derive(serde::Deserialize)]
struct Wrapper {
    isotope_model: IsotopeModelSpec,
}

/// A class name and an explicit table both parse; `s` defaults to absent.
#[test]
fn spec_parses_both_forms() {
    let by_name: Wrapper = toml::from_str(r#"isotope_model = "rna""#).unwrap();
    assert_eq!(by_name.isotope_model.model(), RNA);

    let custom: Wrapper = toml::from_str(
        r#"isotope_model = { residue_mass = 321.2916, c = 9.5, h = 11.75, n = 3.75, o = 7.0 }"#,
    )
    .unwrap();
    assert_eq!(custom.isotope_model.model(), RNA);

    let with_s: Wrapper = toml::from_str(
        r#"isotope_model = { residue_mass = 111.1254, c = 4.9384, h = 7.7583, n = 1.3577, o = 1.4773, s = 0.0417 }"#,
    )
    .unwrap();
    assert_eq!(with_s.isotope_model.model(), PEPTIDE);
}

/// An unknown class name fails the parse rather than falling back to peptide.
#[test]
fn unknown_class_is_a_parse_error() {
    assert!(toml::from_str::<Wrapper>(r#"isotope_model = "protein""#).is_err());
    assert!(toml::from_str::<Wrapper>(r#"isotope_model = { c = 9.5 }"#).is_err());
}

/// Nucleic acids are acquired in negative mode, where the adduct is lost rather
/// than gained. Getting the sign wrong mis-masses a 4-charged ion by 8 Da.
#[test]
fn negative_mode_flips_the_adduct() {
    use crate::models::{Polarity, PROTON_MASS};
    let (mz, z) = (750.0, 4u8);
    let pos = Polarity::Positive.neutral_mass(mz, z);
    let neg = Polarity::Negative.neutral_mass(mz, z);
    assert!((pos - 2_995.970_894).abs() < 1e-6, "{pos}");
    assert!((neg - 3_004.029_106).abs() < 1e-6, "{neg}");
    assert!((neg - pos - 2.0 * z as f64 * PROTON_MASS).abs() < 1e-9);
    assert_eq!(Polarity::default(), Polarity::Positive);
}

/// Polarity is set by name in `[file]`, and defaults to positive so every
/// existing peptide config is unaffected.
#[test]
fn polarity_parses_by_name() {
    use crate::config::FileConfig;
    use crate::models::Polarity;
    assert_eq!(FileConfig::default().polarity, Polarity::Positive);
    #[derive(serde::Deserialize)]
    struct W {
        polarity: Polarity,
    }
    let w: W = toml::from_str(r#"polarity = "negative""#).unwrap();
    assert_eq!(w.polarity, Polarity::Negative);
    assert!(toml::from_str::<W>(r#"polarity = "anion""#).is_err());
}
