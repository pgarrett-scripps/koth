//! Analyte-class isotope models: the average composition koth scores an isotope
//! envelope against.
//!
//! Hill detection and chain assembly assume only that isotopes sit at
//! charge-spaced m/z intervals. Scoring is the single place an analyte class
//! enters, through an "averagine": an average residue composition scaled to a
//! molecule's neutral mass (Senko 1995, for peptides). Carrying that
//! composition as data rather than as constants is what lets a non-peptide
//! analyte be scored against its own model.
//!
//! **Phosphorus is deliberately absent.** ³¹P is the only stable phosphorus
//! isotope, so a phosphate contributes a delta at offset 0 and cannot shift an
//! isotope *pattern*. The nucleic-acid models carry their P in `residue_mass`,
//! which is what sets the scale, and nowhere else — excluding it from the
//! convolution is exact, not an approximation.

use serde::{Deserialize, Serialize};

use super::elements::{cache, K_PATTERN};

/// Average residue composition of an analyte class, expressed as element counts
/// per `residue_mass` Da of neutral mass. Scaling these by
/// `neutral_mass / residue_mass` and rounding gives the composition whose exact
/// isotopologue distribution [`super::elements`] convolves.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IsotopeModel {
    /// Mass of one average residue, in Da. The denominator of every ratio.
    pub residue_mass: f64,
    /// Carbon atoms per residue.
    pub c: f64,
    /// Hydrogen atoms per residue.
    pub h: f64,
    /// Nitrogen atoms per residue.
    pub n: f64,
    /// Oxygen atoms per residue.
    pub o: f64,
    /// Sulfur atoms per residue. Zero for the nucleic-acid models, which
    /// disables sulfur-aware scoring (see [`IsotopeModel::has_sulfur`]).
    #[serde(default)]
    pub s: f64,
}

/// Peptide averagine — Senko's averaged amino acid, C₄.₉₃₈₄ H₇.₇₅₈₃ N₁.₃₅₇₇
/// O₁.₄₇₇₃ S₀.₀₄₁₇ per 111.1254 Da. koth's default and the only model the
/// published benchmark exercises.
pub const PEPTIDE: IsotopeModel = IsotopeModel {
    residue_mass: 111.1254,
    c: 4.9384,
    h: 7.7583,
    n: 1.3577,
    o: 1.4773,
    s: 0.0417,
};

/// RNA — the mean of the four ribonucleotide chain residues (a nucleoside
/// monophosphate less one water, the unit a phosphodiester chain repeats):
///
/// | residue | formula        | monoisotopic |
/// |---------|----------------|--------------|
/// | A       | C₁₀H₁₂N₅O₆P    | 329.05252    |
/// | C       | C₉H₁₂N₃O₇P     | 305.04129    |
/// | G       | C₁₀H₁₂N₅O₇P    | 345.04743    |
/// | U       | C₉H₁₁N₂O₈P     | 306.02530    |
///
/// Unweighted mean: C₉.₅ H₁₁.₇₅ N₃.₇₅ O₇ (P₁) per 321.2916 Da. Unweighted
/// because base composition is sequence-specific; weight it yourself with a
/// custom model if you know the organism's rRNA composition.
pub const RNA: IsotopeModel = IsotopeModel {
    residue_mass: 321.2916,
    c: 9.5,
    h: 11.75,
    n: 3.75,
    o: 7.0,
    s: 0.0,
};

/// DNA — the same construction over the four deoxyribonucleotide residues
/// (dA C₁₀H₁₂N₅O₅P, dC C₉H₁₂N₃O₆P, dG C₁₀H₁₂N₅O₆P, dT C₁₀H₁₃N₂O₇P):
/// C₉.₇₅ H₁₂.₂₅ N₃.₇₅ O₆ (P₁) per 308.8006 Da.
pub const DNA: IsotopeModel = IsotopeModel {
    residue_mass: 308.8006,
    c: 9.75,
    h: 12.25,
    n: 3.75,
    o: 6.0,
    s: 0.0,
};

impl Default for IsotopeModel {
    fn default() -> Self {
        PEPTIDE
    }
}

impl IsotopeModel {
    /// Integer element counts (C, H, N, O, S) at `neutral_mass`. Half-up
    /// rounding from the scaled floats.
    pub fn counts(&self, neutral_mass: f64) -> (u32, u32, u32, u32, u32) {
        let scale = neutral_mass.max(0.0) / self.residue_mass;
        let round = |x: f64| -> u32 { (x + 0.5).floor().max(0.0) as u32 };
        (
            round(self.c * scale),
            round(self.h * scale),
            round(self.n * scale),
            round(self.o * scale),
            round(self.s * scale),
        )
    }

    /// Theoretical isotope distribution at `neutral_mass`, normalized
    /// `[p0, …, p9]`.
    pub fn distribution(&self, neutral_mass: f64) -> [f64; K_PATTERN] {
        let (c, h, n, o, s) = self.counts(neutral_mass);
        cache().distribution(c, h, n, o, s)
    }

    /// As [`Self::distribution`] with the sulfur count overridden: C/H/N/O held
    /// at their scaled values, S set to `n_s`.
    pub fn distribution_with_sulfur(&self, neutral_mass: f64, n_s: u32) -> [f64; K_PATTERN] {
        let (c, h, n, o, _) = self.counts(neutral_mass);
        cache().distribution(c, h, n, o, n_s)
    }

    /// Expected sulfur count rounded **up** rather than half-up. Only the sulfur
    /// term changes; C/H/N/O keep [`Self::counts`]' rounding. For the peptide
    /// model (`s = 0.0417`, i.e. 3.7524e-4·M) the ceiling steps at M ≈ 2665 /
    /// 5330 / 7995 Da.
    pub fn sulfur_ceil(&self, neutral_mass: f64) -> u32 {
        let scale = neutral_mass.max(0.0) / self.residue_mass;
        (self.s * scale).ceil().max(0.0) as u32
    }

    /// Whether the model contains sulfur at all. A model that does not (RNA,
    /// DNA) is scored against one template regardless of `sulfur_offsets`:
    /// varying a count that is structurally zero would only hand every
    /// candidate, decoys included, a free maximum over templates.
    pub fn has_sulfur(&self) -> bool {
        self.s > 0.0
    }
}

/// The built-in analyte classes, selected by name in a config file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum IsotopeClass {
    /// [`PEPTIDE`] — the default.
    Peptide,
    /// [`RNA`].
    Rna,
    /// [`DNA`].
    Dna,
}

/// What a config file may say for `isotope_model`: either a class name
/// (`isotope_model = "rna"`) or an explicit composition table
/// (`isotope_model = { residue_mass = 321.2916, c = 9.5, h = 11.75, n = 3.75,
/// o = 7.0 }`). An unknown class name is a parse error, not a silent fallback.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum IsotopeModelSpec {
    /// One of the built-in classes.
    Class(IsotopeClass),
    /// An explicit composition.
    Custom(IsotopeModel),
}

impl Default for IsotopeModelSpec {
    fn default() -> Self {
        Self::Class(IsotopeClass::Peptide)
    }
}

impl IsotopeModelSpec {
    /// Resolve to the composition to score against.
    pub fn model(&self) -> IsotopeModel {
        match self {
            Self::Class(IsotopeClass::Peptide) => PEPTIDE,
            Self::Class(IsotopeClass::Rna) => RNA,
            Self::Class(IsotopeClass::Dna) => DNA,
            Self::Custom(m) => *m,
        }
    }
}

#[cfg(test)]
#[path = "model_tests.rs"]
mod tests;
