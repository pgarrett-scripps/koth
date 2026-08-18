//! Per-element isotopologue distribution cache.
//!
//! For each of the five "averagine" elements (C, H, N, O, S), pre-compute the
//! exact isotopologue distribution for atom counts `0..=MAX`. At lookup time,
//! we convolve the five cached distributions for a given `(n_C, n_H, n_N,
//! n_O, n_S)` to get the molecule's isotope pattern.
//!
//! This replaces the older per-mass averagine template table:
//!   - flexible: hold C/H/N/O fixed, vary S → sulfur-aware scoring is one knob
//!   - accurate: no 50 Da mass bins
//!   - cheap: ~82 kB cache, ~4 convolutions per lookup (~1 µs)
//!
//! The five element distributions are exact (NOT Poisson approximations) —
//! built by incrementally convolving the single-atom distribution.

use std::sync::OnceLock;

/// First slot in the isotope-pattern vector that we care about.
/// All distributions are kept at length `K_PATTERN` (= 10 isotope peaks).
pub const K_PATTERN: usize = 10;

/// Upper bounds on element counts the cache supports. Anything above gets
/// silently clamped. Sized to cover any tryptic peptide up to ~5050 Da
/// with room for cysteine-rich edge cases (high S count).
const MAX_C: usize = 300;
const MAX_H: usize = 500;
const MAX_N: usize = 100;
const MAX_O: usize = 100;
/// 20 covers everything biological — even a peptide of pure cysteines on a
/// 5000 Da chain is ~50 C and 50 S; the 20-S cap is the highest reasonable
/// bound we'd ever need for tryptic peptides. Clamped above this.
const MAX_S: usize = 20;

/// Single-atom isotope distributions (mass-offset slot → abundance).
/// Index `i` corresponds to mass shift `+i Da` from the monoisotopic peak.
///
/// Sources: IUPAC isotopic compositions for natural terrestrial elements.
const SINGLE_C: [f64; K_PATTERN] = [0.9893, 0.0107, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
const SINGLE_H: [f64; K_PATTERN] = [0.999885, 0.000115, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
const SINGLE_N: [f64; K_PATTERN] = [0.99632, 0.00368, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
const SINGLE_O: [f64; K_PATTERN] = [0.99757, 0.00038, 0.00205, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
// S has isotopes at offsets 0, +1, +2, +4 Da (no +3 — that slot is zero).
const SINGLE_S: [f64; K_PATTERN] = [0.9499, 0.0075, 0.0425, 0.0, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0];

/// Cache of per-element distributions for atom counts `0..=max_n`.
pub struct ElementCache {
    c: Vec<[f64; K_PATTERN]>,
    h: Vec<[f64; K_PATTERN]>,
    n: Vec<[f64; K_PATTERN]>,
    o: Vec<[f64; K_PATTERN]>,
    s: Vec<[f64; K_PATTERN]>,
}

impl ElementCache {
    fn build() -> Self {
        Self {
            c: build_element_table(&SINGLE_C, MAX_C),
            h: build_element_table(&SINGLE_H, MAX_H),
            n: build_element_table(&SINGLE_N, MAX_N),
            o: build_element_table(&SINGLE_O, MAX_O),
            s: build_element_table(&SINGLE_S, MAX_S),
        }
    }

    /// Convolve the five element distributions for the given atom counts.
    /// Out-of-range counts are clamped to the table's max.
    pub fn distribution(
        &self,
        n_c: u32,
        n_h: u32,
        n_n: u32,
        n_o: u32,
        n_s: u32,
    ) -> [f64; K_PATTERN] {
        let c = &self.c[(n_c as usize).min(MAX_C)];
        let h = &self.h[(n_h as usize).min(MAX_H)];
        let n = &self.n[(n_n as usize).min(MAX_N)];
        let o = &self.o[(n_o as usize).min(MAX_O)];
        let s = &self.s[(n_s as usize).min(MAX_S)];
        let mut d = convolve_k(c, h);
        d = convolve_k(&d, n);
        d = convolve_k(&d, o);
        d = convolve_k(&d, s);
        normalize_in_place(&mut d);
        d
    }
}

/// Build the per-element cache: `out[n]` = distribution for `n` atoms of the element.
fn build_element_table(single: &[f64; K_PATTERN], max_n: usize) -> Vec<[f64; K_PATTERN]> {
    let mut table = Vec::with_capacity(max_n + 1);
    let mut zero = [0.0f64; K_PATTERN];
    zero[0] = 1.0;
    table.push(zero);
    for n in 1..=max_n {
        table.push(convolve_k(&table[n - 1], single));
    }
    table
}

/// Fixed-length convolution: out[i+j] += a[i]*b[j], truncated at K_PATTERN.
fn convolve_k(a: &[f64; K_PATTERN], b: &[f64; K_PATTERN]) -> [f64; K_PATTERN] {
    let mut out = [0.0f64; K_PATTERN];
    for (i, &av) in a.iter().enumerate() {
        if av == 0.0 {
            continue;
        }
        for (j, &bv) in b.iter().enumerate() {
            if i + j >= K_PATTERN {
                break;
            }
            out[i + j] += av * bv;
        }
    }
    out
}

fn normalize_in_place(v: &mut [f64; K_PATTERN]) {
    let sum: f64 = v.iter().sum();
    if sum > 0.0 {
        for x in v.iter_mut() {
            *x /= sum;
        }
    }
}

/// Global singleton — built once on first access.
pub fn cache() -> &'static ElementCache {
    static CACHE: OnceLock<ElementCache> = OnceLock::new();
    CACHE.get_or_init(ElementCache::build)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_atom_matches_input() {
        let c = cache();
        // 1 carbon should match SINGLE_C exactly.
        let d = c.distribution(1, 0, 0, 0, 0);
        for i in 0..K_PATTERN {
            assert!(
                (d[i] - SINGLE_C[i]).abs() < 1e-12,
                "1 C: pos {i}, got {}",
                d[i]
            );
        }
    }

    #[test]
    fn zero_atoms_is_delta() {
        let c = cache();
        let d = c.distribution(0, 0, 0, 0, 0);
        assert!((d[0] - 1.0).abs() < 1e-12);
        for i in 1..K_PATTERN {
            assert_eq!(d[i], 0.0);
        }
    }

    #[test]
    fn distribution_sums_to_one() {
        let c = cache();
        // Random-ish averagine-shape composition for a ~1500 Da peptide.
        let d = c.distribution(67, 105, 18, 20, 1);
        let sum: f64 = d.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9, "sum = {sum}");
    }

    #[test]
    fn sulfur_lifts_m_plus_two() {
        let c = cache();
        // Hold C/H/N/O fixed for a ~1500 Da peptide; vary S.
        let no_s = c.distribution(67, 105, 18, 20, 0);
        let two_s = c.distribution(67, 105, 18, 20, 2);
        // S has 4.25% +2 abundance; 2 sulfurs should lift M+2 measurably
        // (~+8 pp on the raw spectrum, less after re-normalization).
        // ³⁴S is 4.25% per atom; for a ~1500 Da peptide the C13/C13 and O18
        // contributions already dominate M+2, so a 2-S override lifts M+2 by
        // ~15% (relative) — small but unambiguous and the right direction.
        assert!(
            two_s[2] > no_s[2] * 1.10,
            "2-S M+2 ({}) should be at least 10% > 0-S M+2 ({})",
            two_s[2],
            no_s[2],
        );
    }

    #[test]
    fn two_carbons_is_convolution_of_singles() {
        let c = cache();
        let one_c = c.distribution(1, 0, 0, 0, 0);
        let two_c = c.distribution(2, 0, 0, 0, 0);
        // Convolve one_c with itself manually
        let mut expected = [0.0f64; K_PATTERN];
        for i in 0..K_PATTERN {
            for j in 0..K_PATTERN - i {
                expected[i + j] += one_c[i] * one_c[j];
            }
        }
        // Normalize for comparison (in case of tiny float drift)
        let s: f64 = expected.iter().sum();
        for x in &mut expected {
            *x /= s;
        }
        for i in 0..K_PATTERN {
            assert!(
                (two_c[i] - expected[i]).abs() < 1e-12,
                "pos {i}: {} vs {}",
                two_c[i],
                expected[i]
            );
        }
    }
}
