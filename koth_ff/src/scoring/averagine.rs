/// Averagine-based theoretical isotope distribution table.
///
/// The averagine model approximates the average amino acid composition:
/// C₄.₉₃₈₄ H₇.₇₅₈₃ N₁.₃₅₇₇ O₁.₄₇₇₃ S₀.₀₄₁₇ per 111.1254 Da
///
/// For a molecule of mass M, element counts are scaled proportionally,
/// then the multinomial isotope distribution is computed via convolution.
use std::sync::OnceLock;

/// Table of (neutral_mass, [p0, p1, ..., p9]) where each entry is the
/// theoretical isotope distribution (normalized to sum=1) for the given mass.
/// Masses span 50 to 5050 in steps of 50.
static TEMPLATE_TABLE: OnceLock<Vec<(f64, [f64; 10])>> = OnceLock::new();

/// Get or build the averagine template table.
pub fn get_table() -> &'static Vec<(f64, [f64; 10])> {
    TEMPLATE_TABLE.get_or_init(build_table)
}

fn build_table() -> Vec<(f64, [f64; 10])> {
    const STEP: f64 = 50.0;
    const MIN: f64 = 50.0;
    const MAX: f64 = 5050.0;

    let mut table = Vec::new();
    let mut mass = MIN;
    while mass <= MAX {
        let dist = compute_averagine_distribution(mass);
        table.push((mass, dist));
        mass += STEP;
    }
    table
}

/// Look up the theoretical isotope template for the given neutral mass.
/// Returns the template for the nearest tabulated mass.
pub fn lookup_template(neutral_mass: f64) -> &'static [f64; 10] {
    let table = get_table();
    // Binary search for closest mass
    let idx = table.partition_point(|(m, _)| *m < neutral_mass);
    let idx = idx.min(table.len() - 1);

    // Compare with the entry before, if available
    if idx > 0 {
        let prev_dist = (table[idx - 1].0 - neutral_mass).abs();
        let curr_dist = (table[idx].0 - neutral_mass).abs();
        if prev_dist < curr_dist {
            return &table[idx - 1].1;
        }
    }
    &table[idx].1
}

/// Compute the isotope distribution for a molecule of given neutral mass
/// using the averagine model.
///
/// Returns probabilities [p(M), p(M+1), ..., p(M+9)] normalized to sum=1.
fn compute_averagine_distribution(mass: f64) -> [f64; 10] {
    // Averagine formula: per 111.1254 Da
    const AVG_MASS: f64 = 111.1254;
    let scale = mass / AVG_MASS;

    // Element counts (scaled)
    let n_c = 4.9384 * scale;
    let n_h = 7.7583 * scale;
    let n_n = 1.3577 * scale;
    let n_o = 1.4773 * scale;
    let n_s = 0.0417 * scale;

    // Isotope abundances for each element
    // C: C12, C13
    let c_dist = element_dist(n_c, &[(0, 0.9893), (1, 0.0107)]);
    // H: H1, D (H2)
    let h_dist = element_dist(n_h, &[(0, 0.999885), (1, 0.000115)]);
    // N: N14, N15
    let n_dist = element_dist(n_n, &[(0, 0.99632), (1, 0.00368)]);
    // O: O16, O17, O18
    let o_dist = element_dist(n_o, &[(0, 0.99757), (1, 0.00038), (2, 0.00205)]);
    // S: S32, S33, S34, S36
    let s_dist = element_dist(n_s, &[(0, 0.9499), (1, 0.0075), (2, 0.0425), (4, 0.0001)]);

    // Convolve all element distributions
    let mut dist = convolve(&c_dist, &h_dist);
    dist = convolve(&dist, &n_dist);
    dist = convolve(&dist, &o_dist);
    dist = convolve(&dist, &s_dist);

    // Take first 10 and normalize
    let mut result = [0.0f64; 10];
    let n = result.len().min(dist.len());
    result[..n].copy_from_slice(&dist[..n]);
    normalize(&mut result);
    result
}

/// Compute isotope distribution for `n` atoms of an element given its isotope masses/abundances.
///
/// `isotopes` is a list of (mass_offset, abundance) pairs.
/// Uses the binomial/multinomial approximation via repeated convolution.
fn element_dist(n: f64, isotopes: &[(usize, f64)]) -> Vec<f64> {
    if n < 1e-9 {
        let mut d = vec![0.0; 10];
        d[0] = 1.0;
        return d;
    }

    // Poisson approximation: for each heavy isotope with abundance p,
    // the number of heavy atoms ~ Poisson(λ = n * p). Uses the full
    // floating-point n, so fractional atom counts are handled correctly.
    let heavy: Vec<(usize, f64)> = isotopes[1..].iter().copied().collect();

    let mut result = vec![0.0f64; 10];
    result[0] = 1.0;

    for &(offset, abundance) in &heavy {
        let lambda = n * abundance;
        // Poisson PMF: P(k) = e^(-λ) * λ^k / k!
        let e_neg_lambda = (-lambda).exp();
        let mut poisson = vec![0.0f64; 10];
        let mut term = e_neg_lambda;
        poisson[0] = term;
        for k in 1..10 {
            term *= lambda / k as f64;
            poisson[k] = term;
        }
        // Spread: k heavy atoms contribute k*offset mass units, so P(k) goes at index k*offset.
        // (A plain shift was wrong — it placed P(0) at index `offset` instead of 0.)
        let mut spread = vec![0.0f64; 10];
        for (k, &p) in poisson.iter().enumerate() {
            let pos = k * offset;
            if pos < 10 {
                spread[pos] += p;
            }
        }
        result = convolve_fixed(&result, &spread);
        normalize_slice(&mut result);
    }

    result
}

/// Convolve two distributions (polynomial multiplication, take first 10 terms).
fn convolve(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0f64; 10];
    for (i, &av) in a.iter().enumerate().take(10) {
        for (j, &bv) in b.iter().enumerate().take(10 - i) {
            result[i + j] += av * bv;
        }
    }
    result
}

/// Fixed-size (10-element) convolution.
fn convolve_fixed(a: &[f64], b: &[f64]) -> Vec<f64> {
    convolve(a, b)
}

fn normalize(v: &mut [f64]) {
    let sum: f64 = v.iter().sum();
    if sum > 0.0 {
        for x in v.iter_mut() {
            *x /= sum;
        }
    }
}

fn normalize_slice(v: &mut [f64]) {
    let sum: f64 = v.iter().sum();
    if sum > 0.0 {
        for x in v.iter_mut() {
            *x /= sum;
        }
    }
}

/// Compute Bhattacharyya coefficient between observed and theoretical distributions.
///
/// BC = Σ sqrt(p_i * q_i) over active positions only.
///
/// The template is scaled to absolute experimental units (template[i] * obs_sum /
/// template_k_sum). Positions where the expected absolute intensity falls below
/// `min_intensity` are excluded from the comparison — theoretical peaks below the
/// noise floor should not be matched against observed peaks. Pass 0.0 to disable.
///
/// Returns BC in [0, 1].
pub fn bhattacharyya_score(obs: &[f64], template: &[f64; 10], min_intensity: f64) -> f64 {
    static ONCE: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    let debug = obs.iter().any(|&x| x > 1e8);
    if debug {
        ONCE.get_or_init(|| {
            let k = obs.len().min(10);
            let obs_sum: f64 = obs.iter().sum();
            let template_k_sum: f64 = template[..k].iter().sum();
            let scale = if template_k_sum > 0.0 { obs_sum / template_k_sum } else { 0.0 };
            let active: Vec<usize> = (0..k)
                .filter(|&i| min_intensity <= 0.0 || template[i] * scale >= min_intensity)
                .collect();
            eprintln!("[bc] k={k} obs_sum={obs_sum:.3e} template_k_sum={template_k_sum:.6} scale={scale:.3e} min_intensity={min_intensity:.3e}");
            eprintln!("[bc] template={:?}", &template[..k]);
            eprintln!("[bc] active={active:?}");
            if !active.is_empty() {
                let obs_active_sum: f64 = active.iter().map(|&i| obs[i]).sum();
                let theo_active_sum: f64 = active.iter().map(|&i| template[i]).sum();
                eprintln!("[bc] obs_active_sum={obs_active_sum:.3e} theo_active_sum={theo_active_sum:.6}");
                let bc: f64 = active.iter().map(|&i| {
                    let p = obs[i] / obs_active_sum;
                    let q = template[i] / theo_active_sum;
                    eprintln!("[bc]   i={i} obs={:.3e} p={p:.4} q={q:.4} sqrt(pq)={:.4}", obs[i], (p*q).sqrt());
                    (p * q).sqrt()
                }).sum();
                eprintln!("[bc] bc={bc:.4} missed_penalty={:.4} final={:.4}", 1.0-template_k_sum, bc*(template_k_sum));
            }
        });
    }

    let k = obs.len().min(10);
    if k == 0 {
        return 0.0;
    }

    let obs_sum: f64 = obs.iter().sum();
    if obs_sum <= 0.0 {
        return 0.0;
    }

    // Missed penalty: template signal beyond the observed k peaks.
    let template_k_sum: f64 = template[..k].iter().sum();
    let missed_penalty = 1.0 - template_k_sum;

    // Scale template to absolute experimental units and build active-position mask.
    let scale = if template_k_sum > 0.0 { obs_sum / template_k_sum } else { 0.0 };
    let active: Vec<usize> = (0..k)
        .filter(|&i| min_intensity <= 0.0 || template[i] * scale >= min_intensity)
        .collect();

    if active.is_empty() {
        return 0.0;
    }

    let obs_active_sum: f64 = active.iter().map(|&i| obs[i]).sum();
    let theo_active_sum: f64 = active.iter().map(|&i| template[i]).sum();
    if obs_active_sum <= 0.0 || theo_active_sum <= 0.0 {
        return 0.0;
    }

    let bc: f64 = active
        .iter()
        .map(|&i| {
            let p = obs[i] / obs_active_sum;
            let q = template[i] / theo_active_sum;
            (p * q).sqrt()
        })
        .sum();

    (bc * (1.0 - missed_penalty)).clamp(0.0, 1.0)
}
