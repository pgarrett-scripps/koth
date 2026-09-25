use crate::scoring::model::IsotopeModelSpec;
use serde::{Deserialize, Serialize};

/// Resolved form of `FeaturesConfig::cosine_anchor` — which hill the
/// chromatographic-cosine gate is measured against during isotope-chain
/// extension. Parsed (case-insensitively) from the config string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CosineAnchor {
    /// Anchor to the immediate predecessor in the chain (legacy koth default).
    Adjacent,
    /// Anchor every isotope to the monoisotope seed hill.
    Seed,
}

impl CosineAnchor {
    /// Parse a `cosine_anchor` string strictly: an unrecognised value is a
    /// config error rather than a silent fallback.
    pub fn parse(s: &str) -> Result<Self, crate::error::Error> {
        match s.to_ascii_lowercase().as_str() {
            "adjacent" => Ok(CosineAnchor::Adjacent),
            "seed" => Ok(CosineAnchor::Seed),
            other => Err(crate::error::Error::Config(format!(
                "invalid features.cosine_anchor `{other}`: expected \"adjacent\" or \"seed\""
            ))),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct FeaturesConfig {
    pub min_charge: u8,
    pub max_charge: u8,
    /// Per-extension chromatographic-cosine threshold used while building an
    /// isotope chain. The candidate hill's cosine vs the seed hill must clear
    /// this to be added; otherwise chain extension stops in that direction.
    pub min_chain_cosine: f64,
    /// Minimum intensity a candidate isotope hill must retain relative to its
    /// chain predecessor, as a fraction: a heavier isotope must be at least
    /// this fraction of its predecessor's intensity or the chain stops. Chain
    /// extension runs upward only (the seed IS the monoisotope hypothesis), so
    /// there is no downward counterpart to this knob. Default 0.01. The
    /// pre-0.3.0 name `right_max_decrease` is still accepted.
    #[serde(alias = "right_max_decrease")]
    pub min_isotope_step_ratio: f64,
    pub max_isotopes: usize,
    /// Per-extension intensity-ratio gate. After a candidate hill clears
    /// `min_chain_cosine`, also check that its apex intensity vs the
    /// predecessor's matches the averagine ratio within ± this many log2
    /// units. Catches column-bleed / contaminant hills that co-elute (so
    /// pass cosine) but have wrong intensity for an isotope. Set to a
    /// very large value to effectively disable. Typical: 1.5 (within
    /// ~2.83× of expected).
    #[serde(default = "default_max_isotope_log2_ratio")]
    pub max_isotope_log2_ratio: f64,
    /// When true (legacy behaviour), isotope-chain extension stops as
    /// soon as the averagine-*predicted* intensity of the next isotope falls
    /// below the per-run noise floor (5th-pct hill intensity × 0.8), computed
    /// from the seed intensity *before the hill is searched for*. This cheaply
    /// trims decayed chain tails, but it also blocks legitimate low-abundance
    /// monoisotopic seeds from ever pairing their M+1: a seed already near the
    /// floor has a predicted M+1 below the floor, so the chain never forms and
    /// the feature collapses to charge 0. When false (the default), this
    /// predicted-intensity break is skipped and chain extension is terminated
    /// purely by *evidence*:
    /// a missing hill (`find_neighbors` empty), `min_chain_cosine`, the
    /// `max_isotope_log2_ratio` intensity-ratio gate, `min_isotope_step_ratio`
    /// (found hill must be ≥ this fraction of the predecessor), and the
    /// averagine-template-length / `max_isotopes` caps. Setting this false
    /// recovers dim 2+/3+ features that have both isotope hills present but
    /// were never paired — recommended when feeding an MS1 search that needs
    /// feature depth. Re-benchmark FDR/quant when flipping.
    #[serde(default = "default_chain_predicted_intensity_gate")]
    pub chain_predicted_intensity_gate: bool,
    /// Minimum number of mutually-overlapping scans two hills must share
    /// before a chromatographic cosine is computed between them; below it the
    /// cosine is 0 and the isotope chain will not extend across that pair.
    /// Default 3 (the long-standing hard gate — keep for PXD003881).
    ///
    /// On fast gradients hills are only 3–5 scans wide, so a dim isotope hill
    /// frequently overlaps its monoisotope by just 1–2 scans and is rejected
    /// outright regardless of how well the shapes agree — the dominant cause of
    /// koth emitting a low-abundance monoisotope as an unpaired charge-0 feature
    /// on such data (biosaur2 requires only ≥1 shared scan). Lowering this to 2
    /// recovers those pairs; pair it with a shorter `hills.min_scans` (a 2-scan
    /// hill can never reach a 3-scan overlap) and re-benchmark quant, since
    /// short-overlap pairs are noisier.
    #[serde(default = "default_min_scan_overlap")]
    pub min_scan_overlap: usize,
    /// Minimum isotope-pattern (Bhattacharyya) score a candidate — original or
    /// truncated — must reach before it may *claim* its hills; below-bar
    /// candidates are dropped so their hills stay free for a better-fitting
    /// feature. Default 0.0 = no gate.
    #[serde(default = "default_exhaustive_min_isotope_score")]
    pub exhaustive_min_isotope_score: f64,
    /// Signal standard deviation of seed-relative log2 isotope-ratio errors.
    /// The broad noise null has fixed sigma 2; valid range is 0 < sigma < 2.
    pub isotope_evidence_ratio_sigma: f64,
    /// Shape of the Beta(shape, 1) co-elution model against a uniform null.
    /// Must be finite and > 1; larger values favor tighter co-elution.
    pub isotope_evidence_cosine_shape: f64,
    /// Legacy compatibility field; no longer affects claim ordering.
    /// The exhaustive resolver always ranks by additive isotope log evidence.
    /// Accepted so existing configuration files continue to deserialize.
    #[serde(default)]
    pub exhaustive_isotope_priority: bool,
    /// Chromatographic-cosine **anchor** for isotope-chain extension: which hill
    /// each candidate isotope's cosine gate is measured against.
    ///
    /// `"seed"` (the default): anchor *every* isotope's cosine to the
    /// monoisotope seed hill (the convention biosaur2 / AlphaPept / Dinosaur all
    /// use). The m/z step target and the intensity-ratio predecessor still step
    /// from the immediate predecessor — only the cosine reference changes.
    /// Rejects a far isotope that co-elutes with its neighbour but not with the
    /// mono. Beat `"adjacent"` on the 20-run PXD003881 cohort (+0.31 pp recall);
    /// see `default_cosine_anchor`.
    ///
    /// `"adjacent"` (legacy): anchor to the immediate predecessor in the chain —
    /// for M+1 that is the seed, for M+k≥2 the previously-claimed isotope hill.
    /// koth's former default; byte-identical to the pre-2026-07 paper output.
    ///
    /// The per-extension cosine that feeds the composite `mean_cosine` is
    /// computed against whichever reference this selects (so the composite is
    /// consistent with the gate). The reported `cosine_score` field and the
    /// exhaustive-resolver rescoring remain adjacent-style regardless.
    #[serde(default = "default_cosine_anchor")]
    pub cosine_anchor: String,
    /// Sulfur-count offsets to score each isotope chain against, relative to
    /// the **ceiling** of the averagine-expected sulfur count. The chain is
    /// scored against one averagine template per offset and the best
    /// Bhattacharyya kept — this corrects the systematic penalty on Cys/Met-rich
    /// peptides whose M+2 is elevated by ³⁴S (4.25 %, +2 Da).
    ///
    /// Default `[-1, 0, 1]`, i.e. `{0, 1, 2}` sulfurs up to 2665 Da, `{1, 2, 3}`
    /// to 5330 Da, `{2, 3, 4}` to 7995 Da. Negative results saturate at 0 and
    /// duplicates collapse, so `[-2, -1, 0, 1]` on a small peptide is `{0, 1, 2}`.
    /// Widen it (e.g. `[-2, -1, 0, 1]`) to keep the no-sulfur template on large
    /// peptides — ~32 % of 3 kDa tryptic peptides have no sulfur at all.
    ///
    /// An **empty list disables sulfur awareness**: one plain averagine template
    /// at the half-up rounded count.
    ///
    /// Scoring takes a **max over templates**, so a longer list can only raise
    /// scores — including for decoys and mis-assembled chains. Judge a change on
    /// discrimination (PSM recall / target–decoy separation), never on the
    /// isotope-score distribution alone.
    #[serde(default = "default_sulfur_offsets")]
    pub sulfur_offsets: Vec<i8>,
    /// Which analyte class's average composition the theoretical isotope
    /// pattern is built from.
    ///
    /// `"peptide"` (the default) is Senko's averagine and the only model the
    /// published benchmark exercises. `"rna"` and `"dna"` are the means of the
    /// four ribo- / deoxyribonucleotide chain residues. An explicit table
    /// overrides both:
    ///
    /// ```toml
    /// isotope_model = { residue_mass = 321.2916, c = 9.5, h = 11.75, n = 3.75, o = 7.0 }
    /// ```
    ///
    /// `s` defaults to 0 and a model without sulfur ignores `sulfur_offsets`
    /// entirely — varying a structurally absent atom would hand every candidate,
    /// decoys included, a free maximum over templates. An unknown class name is
    /// a config parse error, never a silent fallback to peptide.
    #[serde(default)]
    pub isotope_model: IsotopeModelSpec,
    /// Neutron (C13) mass in Da
    pub neutron_mass: f64,
    /// Drop features whose **isotope_score** (Bhattacharyya vs averagine) is
    /// below this. 0.0 = keep all. Default 0.5.
    pub min_isotope_score: f64,
    /// Drop features whose **cosine_score** (mean chromatographic cosine of
    /// adjacent isotope hills) is below this. 0.0 = keep all.
    pub min_cosine_score: f64,
    /// Drop features whose **combined_score** (= isotope × cosine) is below
    /// this. 0.0 = keep all. All three filters are AND-ed.
    pub min_combined_score: f64,
}

fn default_max_isotope_log2_ratio() -> f64 {
    1.5
}

fn default_chain_predicted_intensity_gate() -> bool {
    // `false` since the downward (M-1, M-2, …) chain extension was removed.
    //
    // The gate exists only in the upward loop, so while the downward walk
    // existed it could route around it: a dim monoisotope whose own upward walk
    // was killed by the predicted-intensity break was still recovered by seeding
    // on its M+1 and stepping down (the assembler's `uncontested_envelope_assembled`
    // fixture is exactly this — every upward walk there breaks at step 1, and
    // the whole 3-hill envelope used to be built by seeding the M+2 and walking
    // down twice). With one direction only there is no second route, so leaving
    // the gate on silently drops those features. Termination is now purely by
    // evidence — see the field docs.
    false
}

fn default_min_scan_overlap() -> usize {
    3
}

fn default_exhaustive_min_isotope_score() -> f64 {
    0.0
}

fn default_cosine_anchor() -> String {
    // `seed` is the default: anchoring every isotope's chromatographic cosine
    // to the monoisotope seed (as biosaur2 / AlphaPept / Dinosaur all do) beat
    // the former `adjacent` (predecessor) anchor on the full 20-run PXD003881
    // cohort, at +0.31 pp recall (+1565 covered PSMs) with no quant regression
    // (median CV, MV rate, and HUMAN FPR all flat-to-better).
    //
    // The absolute recall figures that used to be quoted here (0.7933 ->
    // 0.7964) predate `b83e1e1` and no longer describe any current run. The
    // current cohort recall at this default is 0.7943; the matching `adjacent`
    // arm has not been re-measured since, so the +0.31 pp delta above is
    // historical and should be re-derived before being cited anywhere.
    // Set to `adjacent` to reproduce the pre-2026-07 paper feature output.
    "seed".to_string()
}

pub(crate) fn default_sulfur_offsets() -> Vec<i8> {
    vec![-1, 0, 1]
}

impl Default for FeaturesConfig {
    fn default() -> Self {
        Self {
            min_charge: 2,
            max_charge: 6,
            min_chain_cosine: 0.4,
            min_isotope_step_ratio: 0.01,
            max_isotopes: 6,
            max_isotope_log2_ratio: 1.5,
            chain_predicted_intensity_gate: false,
            min_scan_overlap: 3,
            exhaustive_min_isotope_score: 0.0,
            exhaustive_isotope_priority: false,
            isotope_evidence_ratio_sigma: 0.75,
            isotope_evidence_cosine_shape: 2.0,
            cosine_anchor: default_cosine_anchor(),
            sulfur_offsets: default_sulfur_offsets(),
            isotope_model: IsotopeModelSpec::default(),
            neutron_mass: 1.003_354_835,
            min_isotope_score: 0.5,
            min_cosine_score: 0.0,
            min_combined_score: 0.0,
        }
    }
}

impl FeaturesConfig {
    /// Validate distribution parameters and strings that serde alone cannot
    /// check. Called after TOML deserialization so invalid values fail at load time
    /// rather than silently falling back at runtime.
    pub fn validate(&self) -> Result<(), crate::error::Error> {
        CosineAnchor::parse(&self.cosine_anchor)?;
        if !self.isotope_evidence_ratio_sigma.is_finite()
            || self.isotope_evidence_ratio_sigma <= 0.0
            || self.isotope_evidence_ratio_sigma >= 2.0
        {
            return Err(crate::error::Error::Config(
                "features.isotope_evidence_ratio_sigma must be finite and between 0 and 2 (exclusive)".into(),
            ));
        }
        if !self.isotope_evidence_cosine_shape.is_finite()
            || self.isotope_evidence_cosine_shape <= 1.0
        {
            return Err(crate::error::Error::Config(
                "features.isotope_evidence_cosine_shape must be finite and greater than 1".into(),
            ));
        }
        Ok(())
    }

    /// Resolve the chromatographic-cosine anchor. Infallible at runtime because
    /// `validate()` has already rejected bad values at config load; the
    /// defensive fallback preserves legacy behaviour if reached.
    pub fn cosine_anchor_mode(&self) -> CosineAnchor {
        CosineAnchor::parse(&self.cosine_anchor).unwrap_or(CosineAnchor::Adjacent)
    }
}
