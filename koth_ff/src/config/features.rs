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
    pub fn parse(s: &str) -> Result<Self, crate::error::KothError> {
        match s.to_ascii_lowercase().as_str() {
            "adjacent" => Ok(CosineAnchor::Adjacent),
            "seed" => Ok(CosineAnchor::Seed),
            other => Err(crate::error::KothError::ConfigError(format!(
                "invalid features.cosine_anchor `{other}`: expected \"adjacent\" or \"seed\""
            ))),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeaturesConfig {
    pub min_charge: u8,
    pub max_charge: u8,
    /// Per-extension chromatographic-cosine threshold used while building an
    /// isotope chain. The candidate hill's cosine vs the seed hill must clear
    /// this to be added; otherwise chain extension stops in that direction.
    pub min_chain_cosine: f64,
    /// Minimum intensity a candidate isotope hill must retain relative to its
    /// chain predecessor, as a fraction. Chain extension runs upward only
    /// (the seed IS the monoisotope hypothesis), so there is no downward
    /// counterpart to this knob.
    pub right_max_decrease: f64,
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
    /// When true (default, legacy behaviour), isotope-chain extension stops as
    /// soon as the averagine-*predicted* intensity of the next isotope falls
    /// below the per-run noise floor (5th-pct hill intensity × 0.8), computed
    /// from the seed intensity *before the hill is searched for*. This cheaply
    /// trims decayed chain tails, but it also blocks legitimate low-abundance
    /// monoisotopic seeds from ever pairing their M+1: a seed already near the
    /// floor has a predicted M+1 below the floor, so the chain never forms and
    /// the feature collapses to charge 0. When false, this predicted-intensity
    /// break is skipped and chain extension is terminated purely by *evidence*:
    /// a missing hill (`find_neighbors` empty), `min_chain_cosine`, the
    /// `max_isotope_log2_ratio` intensity-ratio gate, `right_max_decrease`
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
    /// When true, contested-hill claim priority is ordered by envelope length,
    /// then isotope-pattern score, then composite — so the best averagine fit
    /// wins a shared hill within a length class. Default false.
    #[serde(default)]
    pub exhaustive_isotope_priority: bool,
    /// Chromatographic-cosine **anchor** for isotope-chain extension: which hill
    /// each candidate isotope's cosine gate is measured against.
    ///
    /// `"adjacent"` (default, legacy): anchor to the immediate predecessor in
    /// the chain — for M+1 that is the seed, for M+k≥2 the previously-claimed
    /// isotope hill. This is koth's deliberate choice: chains drift in S/N as
    /// they extend from the seed, so anchoring far isotopes to the seed
    /// over-rejects them. Default = byte-identical to legacy koth.
    ///
    /// `"seed"`: anchor *every* isotope's cosine to the monoisotope seed hill
    /// (the convention biosaur2 / AlphaPept / Dinosaur all use). The m/z step
    /// target and the intensity-ratio predecessor still step from the immediate
    /// predecessor — only the cosine reference changes. Rejects a far isotope
    /// that co-elutes with its neighbour but not with the mono.
    ///
    /// The per-extension cosine that feeds the composite `mean_cosine` is
    /// computed against whichever reference this selects (so the composite is
    /// consistent with the gate). The reported `cosine_score` field and the
    /// exhaustive-resolver rescoring remain adjacent-style regardless.
    #[serde(default = "default_cosine_anchor")]
    pub cosine_anchor: String,
    /// Score isotope chains against multiple averagine templates that vary
    /// the sulfur atom count `{0, avg, avg+2, avg+4}` and keep the best fit.
    /// Corrects the systematic Bhattacharyya penalty on Cys/Met-rich
    /// peptides whose M+2 is elevated by ³⁴S (4.25%, +2 Da).
    /// Default true.
    #[serde(default = "default_sulfur_aware_scoring")]
    pub sulfur_aware_scoring: bool,
    /// Neutron (C13) mass in Da
    pub neutron_mass: f64,
    /// Drop features whose **isotope_score** (Bhattacharyya vs averagine) is
    /// below this. 0.0 = keep all.
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
    // cohort — recall 0.7933 -> 0.7964 (+0.31 pp, +1565 covered PSMs) with no
    // quant regression (median CV, MV rate, and HUMAN FPR all flat-to-better).
    // Set to `adjacent` to reproduce the pre-2026-07 paper feature output.
    "seed".to_string()
}

fn default_sulfur_aware_scoring() -> bool {
    true
}

impl Default for FeaturesConfig {
    fn default() -> Self {
        Self {
            min_charge: 1,
            max_charge: 7,
            min_chain_cosine: 0.5,
            right_max_decrease: 0.05,
            max_isotopes: 6,
            max_isotope_log2_ratio: 1.5,
            chain_predicted_intensity_gate: false,
            min_scan_overlap: 3,
            exhaustive_min_isotope_score: 0.0,
            exhaustive_isotope_priority: false,
            cosine_anchor: default_cosine_anchor(),
            sulfur_aware_scoring: true,
            neutron_mass: 1.003_354_835,
            min_isotope_score: 0.0,
            min_cosine_score: 0.0,
            min_combined_score: 0.0,
        }
    }
}

impl FeaturesConfig {
    /// Validate string-valued knobs that serde alone cannot check. Called after
    /// TOML deserialization so an unrecognised value is rejected at load time
    /// rather than silently falling back at runtime.
    pub fn validate(&self) -> Result<(), crate::error::KothError> {
        CosineAnchor::parse(&self.cosine_anchor)?;
        Ok(())
    }

    /// Resolve the chromatographic-cosine anchor. Infallible at runtime because
    /// `validate()` has already rejected bad values at config load; the
    /// defensive fallback preserves legacy behaviour if reached.
    pub fn cosine_anchor_mode(&self) -> CosineAnchor {
        CosineAnchor::parse(&self.cosine_anchor).unwrap_or(CosineAnchor::Adjacent)
    }
}
