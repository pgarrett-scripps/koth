//! Hill-splitting experimentation lab.
//!
//! Standalone — does NOT import koth internals. It re-implements the CURRENT
//! production splitter (faithful to `src/hills/split.rs`) plus several candidate
//! algorithms, then scores them against a suite of synthetic noisy profiles with
//! known ground-truth peak counts (single, dual, shoulder, spike-contaminated…).
//!
//! Run:  cargo run -p koth_ff --example hill_split_lab
//!
//! Noise is deterministic (seeded LCG) so results are reproducible run-to-run.

use std::f32::consts::PI;

// ----------------------------------------------------------------------------
// Deterministic RNG (LCG + Box-Muller) — reproducible without external crates.
// ----------------------------------------------------------------------------
struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed ^ 0x9E37_79B9_7F4A_7C15)
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }
    fn unit(&mut self) -> f32 {
        ((self.next_u64() >> 40) as f32) / ((1u64 << 24) as f32) // [0,1)
    }
    fn gauss(&mut self) -> f32 {
        let u1 = self.unit().max(1e-7);
        let u2 = self.unit();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
}

// ----------------------------------------------------------------------------
// Synthetic profile construction.
// ----------------------------------------------------------------------------
fn add_gaussian(p: &mut [f32], center: f32, sigma: f32, amp: f32) {
    for (i, v) in p.iter_mut().enumerate() {
        let d = (i as f32 - center) / sigma;
        *v += amp * (-0.5 * d * d).exp();
    }
}

/// Apply realistic, valley-filling noise:
///  - multiplicative jitter (signal-proportional)
///  - a positive additive floor (|gauss|*floor) — this is what lifts valleys
///  - occasional sharp spikes (corrupts the global max)
fn apply_noise(p: &mut [f32], rel: f32, floor: f32, spike_prob: f32, spike_mult: f32, rng: &mut Rng) {
    for v in p.iter_mut() {
        let mult = 1.0 + rel * rng.gauss();
        let add = floor * rng.gauss().abs();
        let mut x = *v * mult + add;
        if rng.unit() < spike_prob {
            x *= spike_mult;
        }
        *v = x.max(0.0);
    }
}

// ----------------------------------------------------------------------------
// Shared primitives.
// ----------------------------------------------------------------------------
fn smooth(p: &[f32], half: usize) -> Vec<f32> {
    if half == 0 || p.len() < 2 {
        return p.to_vec();
    }
    let n = p.len();
    (0..n)
        .map(|i| {
            let lo = i.saturating_sub(half);
            let hi = (i + half + 1).min(n);
            p[lo..hi].iter().sum::<f32>() / (hi - lo) as f32
        })
        .collect()
}

/// Width-3 median filter (nonlinear despike). Removes lone 1-sample spikes
/// completely while preserving genuine peak edges — something no linear
/// (moving-average) smoother can do.
fn median3(p: &[f32]) -> Vec<f32> {
    let n = p.len();
    if n < 3 {
        return p.to_vec();
    }
    let mut out = p.to_vec();
    for i in 1..n - 1 {
        let mut t = [p[i - 1], p[i], p[i + 1]];
        t.sort_by(|a, b| a.partial_cmp(b).unwrap());
        out[i] = t[1];
    }
    out
}

fn percentile(p: &[f32], q: f32) -> f32 {
    let mut s = p.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let idx = ((s.len() as f32 - 1.0) * q).round() as usize;
    s[idx]
}

/// Robust noise σ via the median absolute first-difference (MAD of the
/// derivative). Peaks contribute a few large diffs; the median ignores them.
fn noise_sigma(p: &[f32]) -> f32 {
    if p.len() < 3 {
        return 0.0;
    }
    let mut d: Vec<f32> = p.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = d[d.len() / 2];
    // successive differences inflate σ by √2; 1.4826 converts MAD→σ for Gaussian
    med * 1.4826 / std::f32::consts::SQRT_2
}

fn local_maxima(p: &[f32], min_height: f32) -> Vec<usize> {
    let n = p.len();
    let mut out = Vec::new();
    for i in 0..n {
        let l = i == 0 || p[i] >= p[i - 1];
        let r = i == n - 1 || p[i] >= p[i + 1];
        if p[i] >= min_height && l && r {
            out.push(i);
        }
    }
    out
}

/// Keep the tallest candidate within `min_distance` (faithful to production).
fn enforce_spacing(cands: &[usize], p: &[f32], min_distance: usize) -> Vec<usize> {
    let mut spaced: Vec<usize> = Vec::new();
    'outer: for &idx in cands {
        let mut remove = None;
        for (ri, &ex) in spaced.iter().enumerate() {
            if idx.abs_diff(ex) < min_distance {
                if p[idx] <= p[ex] {
                    continue 'outer;
                } else {
                    remove = Some(ri);
                    break;
                }
            }
        }
        if let Some(ri) = remove {
            spaced.remove(ri);
        }
        spaced.push(idx);
    }
    spaced.sort_unstable();
    spaced
}

fn prominence(p: &[f32], peak: usize) -> f32 {
    let h = p[peak];
    let mut lmin = h;
    for i in (0..peak).rev() {
        if p[i] > h {
            break;
        }
        lmin = lmin.min(p[i]);
    }
    let mut rmin = h;
    for i in (peak + 1)..p.len() {
        if p[i] > h {
            break;
        }
        rmin = rmin.min(p[i]);
    }
    h - lmin.max(rmin)
}

/// Lowest valley value strictly between two indices.
fn valley(p: &[f32], a: usize, b: usize) -> f32 {
    p[a..=b].iter().cloned().fold(f32::INFINITY, f32::min)
}

#[derive(Clone, Copy)]
struct Cfg {
    min_peak_distance: usize,
    min_peak_height: f32,
    min_prominence: f32,
    min_scans: usize,
}
const DEFAULT: Cfg = Cfg {
    min_peak_distance: 10,
    min_peak_height: 0.2,
    min_prominence: 0.2,
    min_scans: 3,
};

/// Count resulting segments after cutting at the valley between each adjacent
/// peak pair, dropping any segment shorter than `min_scans`.
fn enforce_segments(peaks: &[usize], p: &[f32], min_scans: usize) -> usize {
    if peaks.len() <= 1 {
        return 1;
    }
    let mut bounds = vec![0usize];
    for w in peaks.windows(2) {
        let v = (w[0]..=w[1])
            .min_by(|&a, &b| p[a].partial_cmp(&p[b]).unwrap())
            .unwrap();
        bounds.push(v);
    }
    bounds.push(p.len());
    let segs = bounds.windows(2).filter(|w| w[1] - w[0] >= min_scans).count();
    segs.max(1)
}

// ----------------------------------------------------------------------------
// ALGORITHMS — each returns the number of resulting hills (segments).
// ----------------------------------------------------------------------------

/// CURRENT production algorithm (faithful port of src/hills/split.rs).
fn algo_current(p: &[f32], c: Cfg) -> usize {
    if p.len() < c.min_scans * 2 {
        return 1;
    }
    let maxi = p.iter().cloned().fold(0.0f32, f32::max);
    if maxi == 0.0 {
        return 1;
    }
    let min_h = maxi * c.min_peak_height;
    let prom_t = maxi * c.min_prominence;
    let cands = local_maxima(p, min_h);
    let spaced = enforce_spacing(&cands, p, c.min_peak_distance);
    let peaks: Vec<usize> = spaced
        .into_iter()
        .filter(|&i| prominence(p, i) >= prom_t)
        .collect();
    if peaks.len() <= 1 {
        return 1;
    }
    enforce_segments(&peaks, p, c.min_scans)
}

/// CANDIDATE A — detect on a smoothed copy + robust (95th-pct) threshold anchor.
/// Minimal change from current: only the noise-sensitive parts are hardened.
fn algo_smooth_robust(p: &[f32], c: Cfg) -> usize {
    if p.len() < c.min_scans * 2 {
        return 1;
    }
    let s = smooth(p, 2);
    let robust_max = percentile(&s, 0.95).max(1e-6);
    let min_h = robust_max * c.min_peak_height;
    let prom_t = robust_max * c.min_prominence;
    let cands = local_maxima(&s, min_h);
    let spaced = enforce_spacing(&cands, &s, c.min_peak_distance);
    let peaks: Vec<usize> = spaced
        .into_iter()
        .filter(|&i| prominence(&s, i) >= prom_t)
        .collect();
    if peaks.len() <= 1 {
        return 1;
    }
    enforce_segments(&peaks, &s, c.min_scans)
}

/// CANDIDATE B — scale-free valley-ratio (agglomerative / persistence-style).
/// Keep a split between adjacent peaks iff the valley between them drops below
/// `ratio * min(peak_left, peak_right)`. Immune to global-max spike inflation;
/// naturally fair to unequal-height duals (compares against the *smaller* peak).
fn algo_valley_ratio(p: &[f32], c: Cfg, ratio: f32) -> usize {
    if p.len() < c.min_scans * 2 {
        return 1;
    }
    let s = smooth(p, 2);
    let robust_max = percentile(&s, 0.95).max(1e-6);
    let min_h = robust_max * c.min_peak_height;
    let cands = local_maxima(&s, min_h);
    let mut peaks = enforce_spacing(&cands, &s, c.min_peak_distance);
    // Agglomerative merge: repeatedly drop the worst-separated peak.
    loop {
        if peaks.len() <= 1 {
            break;
        }
        // Find the adjacent pair with the SHALLOWEST notch relative to its
        // smaller peak (highest valley/min-peak ratio = least separated).
        let mut worst = None;
        for (k, w) in peaks.windows(2).enumerate() {
            let v = valley(&s, w[0], w[1]);
            let smaller = s[w[0]].min(s[w[1]]);
            let r = v / smaller.max(1e-6);
            if worst.map_or(true, |(_, wr)| r > wr) {
                worst = Some((k, r));
            }
        }
        let (k, r) = worst.unwrap();
        if r <= ratio {
            break; // every notch is deep enough — done
        }
        // Merge: remove the smaller of the two peaks in the offending pair.
        let (a, b) = (peaks[k], peaks[k + 1]);
        let drop = if s[a] <= s[b] { k } else { k + 1 };
        peaks.remove(drop);
    }
    if peaks.len() <= 1 {
        return 1;
    }
    enforce_segments(&peaks, &s, c.min_scans)
}

/// CANDIDATE C — noise-aware absolute prominence. Estimate σ from the trace
/// (MAD of derivative) and require prominence > k·σ in ABSOLUTE units, plus the
/// valley to sit k·σ below both peaks. No dependence on the global max at all.
fn algo_noise_aware(p: &[f32], c: Cfg, k: f32) -> usize {
    if p.len() < c.min_scans * 2 {
        return 1;
    }
    let s = smooth(p, 2);
    let sigma = noise_sigma(p).max(1e-6);
    let floor = k * sigma;
    let robust_max = percentile(&s, 0.95).max(1e-6);
    let min_h = (robust_max * c.min_peak_height).max(floor);
    let cands = local_maxima(&s, min_h);
    let mut peaks = enforce_spacing(&cands, &s, c.min_peak_distance);
    // Merge adjacent peaks whose separating notch is < k·σ below the smaller peak.
    loop {
        if peaks.len() <= 1 {
            break;
        }
        let mut worst: Option<(usize, f32)> = None; // (k, notch_depth) smallest depth
        for (idx, w) in peaks.windows(2).enumerate() {
            let v = valley(&s, w[0], w[1]);
            let depth = s[w[0]].min(s[w[1]]) - v; // how far the notch drops below smaller peak
            if worst.map_or(true, |(_, d)| depth < d) {
                worst = Some((idx, depth));
            }
        }
        let (idx, depth) = worst.unwrap();
        if depth >= floor {
            break;
        }
        let (a, b) = (peaks[idx], peaks[idx + 1]);
        let drop = if s[a] <= s[b] { idx } else { idx + 1 };
        peaks.remove(drop);
    }
    if peaks.len() <= 1 {
        return 1;
    }
    enforce_segments(&peaks, &s, c.min_scans)
}

/// CANDIDATE D — the kitchen sink. Combines every robustness idea:
///   1. median-3 despike  (kills lone spikes — fixes the "spike" cases)
///   2. moving-average smooth (half=2)
///   3. noise σ estimated from the trace (MAD of derivative)
///   4. NOISE-RELATIVE height floor = max(0.05·robust_max, k·σ)  — lets a small
///      but real peak (10:1 duals) through while still rejecting noise bumps
///   5. NO hard min-distance eviction — close real peaks survive
///   6. agglomerative merge: drop a peak whenever its separating notch is
///      shallower than k·σ below the smaller neighbour (persistence/valley test)
fn algo_best(p: &[f32], c: Cfg, k: f32) -> usize {
    let (peaks, s) = best_peaks(p, c, k);
    if peaks.len() <= 1 {
        return 1;
    }
    enforce_segments(&peaks, &s, c.min_scans)
}

/// Core of CANDIDATE D — returns (final peak indices, smoothed profile) so the
/// debug path can inspect exactly what it found.
fn best_peaks(p: &[f32], c: Cfg, k: f32) -> (Vec<usize>, Vec<f32>) {
    if p.len() < c.min_scans * 2 {
        return (vec![0], p.to_vec());
    }
    let despiked = median3(p);
    let s = smooth(&despiked, 2);
    let sigma = noise_sigma(&despiked).max(1e-6);
    let floor = k * sigma;
    let robust_max = percentile(&s, 0.95).max(1e-6);
    // Height floor is the dominant lever: low (0.05) catches faint 10:1 minor
    // peaks but admits baseline bumps; 0.10 is the practical sweet spot.
    let min_h = (0.10 * robust_max).max(floor);

    let cands = local_maxima(&s, min_h);
    let mut peaks = cands; // deliberately NO enforce_spacing
    // A split between adjacent peaks is kept iff BOTH hold:
    //   (relative) valley / smaller_peak <= ratio   -> fair to unequal duals,
    //                                                   immune to spike inflation
    //   (absolute) smaller_peak - valley >= k·σ     -> rejects noise notches
    // The relative test kills noise-on-a-peak-top (valley sits near the apex);
    // the absolute test kills shallow wiggles on a low baseline.
    let ratio = 0.70f32;
    loop {
        if peaks.len() <= 1 {
            break;
        }
        // Merge the single worst-separated pair this round (the one whose notch
        // most badly fails either test), then re-evaluate.
        let mut merge_at: Option<(usize, f32)> = None; // (idx, badness)
        for (idx, w) in peaks.windows(2).enumerate() {
            let v = valley(&s, w[0], w[1]);
            let smaller = s[w[0]].min(s[w[1]]).max(1e-6);
            let r = v / smaller; // higher = shallower notch
            let depth = smaller - v;
            let fails = r > ratio || depth < floor;
            if fails {
                // badness: how far over the ratio, plus how far under the floor
                let badness = (r - ratio).max(0.0) + (floor - depth).max(0.0) / floor;
                if merge_at.map_or(true, |(_, b)| badness > b) {
                    merge_at = Some((idx, badness));
                }
            }
        }
        let Some((idx, _)) = merge_at else { break };
        let (a, b) = (peaks[idx], peaks[idx + 1]);
        let drop = if s[a] <= s[b] { idx } else { idx + 1 };
        peaks.remove(drop);
    }
    // Width gate: a real chromatographic peak is several scans wide at half-max.
    // Drops surviving 2-sample spikes that median-3 alone can't remove.
    let min_w = 4usize;
    peaks.retain(|&pk| {
        let half = 0.5 * s[pk];
        let lo = pk.saturating_sub(8);
        let hi = (pk + 9).min(s.len());
        (lo..hi).filter(|&i| s[i] >= half).count() >= min_w
    });
    (peaks, s)
}

// ----------------------------------------------------------------------------
// Test suite.
// ----------------------------------------------------------------------------
struct Case {
    name: &'static str,
    expected: usize,
    build: fn(&mut Rng) -> Vec<f32>,
}

fn cases() -> Vec<Case> {
    vec![
        Case {
            name: "single clean",
            expected: 1,
            build: |_r| {
                let mut p = vec![0.0; 80];
                add_gaussian(&mut p, 40.0, 6.0, 1000.0);
                p
            },
        },
        Case {
            name: "single noisy",
            expected: 1,
            build: |r| {
                let mut p = vec![0.0; 80];
                add_gaussian(&mut p, 40.0, 6.0, 1000.0);
                apply_noise(&mut p, 0.18, 60.0, 0.04, 1.6, r);
                p
            },
        },
        Case {
            name: "single + big spike",
            expected: 1,
            build: |r| {
                let mut p = vec![0.0; 80];
                add_gaussian(&mut p, 40.0, 6.0, 1000.0);
                apply_noise(&mut p, 0.12, 50.0, 0.0, 1.0, r);
                p[12] += 1800.0; // lone spike, taller than the real apex
                p
            },
        },
        Case {
            name: "dual clean equal",
            expected: 2,
            build: |_r| {
                let mut p = vec![0.0; 100];
                add_gaussian(&mut p, 35.0, 5.0, 1000.0);
                add_gaussian(&mut p, 65.0, 5.0, 1000.0);
                p
            },
        },
        Case {
            name: "dual clean 10:1",
            expected: 2,
            build: |_r| {
                let mut p = vec![0.0; 100];
                add_gaussian(&mut p, 35.0, 5.0, 1000.0);
                add_gaussian(&mut p, 65.0, 5.0, 100.0);
                p
            },
        },
        Case {
            name: "dual noisy equal",
            expected: 2,
            build: |r| {
                let mut p = vec![0.0; 100];
                add_gaussian(&mut p, 35.0, 5.0, 1000.0);
                add_gaussian(&mut p, 65.0, 5.0, 1000.0);
                apply_noise(&mut p, 0.20, 90.0, 0.04, 1.6, r);
                p
            },
        },
        Case {
            name: "dual noisy 5:1",
            expected: 2,
            build: |r| {
                let mut p = vec![0.0; 100];
                add_gaussian(&mut p, 35.0, 5.0, 1000.0);
                add_gaussian(&mut p, 65.0, 5.0, 200.0);
                apply_noise(&mut p, 0.18, 70.0, 0.03, 1.6, r);
                p
            },
        },
        Case {
            name: "dual + spike on one",
            expected: 2,
            build: |r| {
                let mut p = vec![0.0; 100];
                add_gaussian(&mut p, 35.0, 5.0, 1000.0);
                add_gaussian(&mut p, 65.0, 5.0, 1000.0);
                apply_noise(&mut p, 0.15, 70.0, 0.0, 1.0, r);
                p[35] += 1500.0; // spike inflates global max -> raises thresholds
                p
            },
        },
        Case {
            name: "shoulder (TN)",
            expected: 1,
            build: |r| {
                // bump on the flank: valley only dips to ~0.75 of the smaller -> NOT a split
                let mut p = vec![0.0; 90];
                add_gaussian(&mut p, 45.0, 7.0, 1000.0);
                add_gaussian(&mut p, 60.0, 5.0, 350.0);
                apply_noise(&mut p, 0.10, 40.0, 0.0, 1.0, r);
                p
            },
        },
        Case {
            name: "dual close (8 scans)",
            expected: 2,
            build: |r| {
                // narrow peaks (sigma=2) 8 scans apart -> deep notch (~0.27 of
                // apex): clearly TWO peaks, but closer than min_peak_distance=10.
                let mut p = vec![0.0; 80];
                add_gaussian(&mut p, 36.0, 2.0, 1000.0);
                add_gaussian(&mut p, 44.0, 2.0, 900.0);
                apply_noise(&mut p, 0.12, 50.0, 0.0, 1.0, r);
                p
            },
        },
        Case {
            name: "triple noisy",
            expected: 3,
            build: |r| {
                let mut p = vec![0.0; 130];
                add_gaussian(&mut p, 30.0, 5.0, 1000.0);
                add_gaussian(&mut p, 65.0, 5.0, 800.0);
                add_gaussian(&mut p, 100.0, 5.0, 1000.0);
                apply_noise(&mut p, 0.18, 80.0, 0.03, 1.6, r);
                p
            },
        },
    ]
}

fn main() {
    let c = DEFAULT;
    #[allow(clippy::type_complexity)]
    let algos: Vec<(&str, Box<dyn Fn(&[f32]) -> usize>)> = vec![
        ("current", Box::new(move |p: &[f32]| algo_current(p, c))),
        ("smooth+rob", Box::new(move |p: &[f32]| algo_smooth_robust(p, c))),
        ("valley", Box::new(move |p: &[f32]| algo_valley_ratio(p, c, 0.70))),
        ("noise-aware", Box::new(move |p: &[f32]| algo_noise_aware(p, c, 4.0))),
        ("BEST", Box::new(move |p: &[f32]| algo_best(p, c, 4.0))),
    ];
    let w = 12usize; // column width

    // Average over several noise seeds so the verdict isn't one unlucky draw.
    const SEEDS: u64 = 40;

    print!("\n{:<22} {:>4} ", "case", "exp");
    for (n, _) in &algos {
        print!(" {:>w$}", n);
    }
    println!();
    println!("{}", "-".repeat(28 + (w + 1) * algos.len()));

    let mut totals = vec![0u32; algos.len()];

    for case in cases() {
        let mut hits = vec![0u32; algos.len()];
        for seed in 0..SEEDS {
            let mut rng = Rng::new(seed.wrapping_mul(2654435761) ^ case.name.len() as u64);
            let prof = (case.build)(&mut rng);
            for (ai, (_n, f)) in algos.iter().enumerate() {
                if f(&prof) == case.expected {
                    hits[ai] += 1;
                }
            }
        }
        print!("{:<22} {:>4} ", case.name, case.expected);
        for (ai, h) in hits.iter().enumerate() {
            let pct = 100.0 * *h as f32 / SEEDS as f32;
            let mark = if pct >= 90.0 { "✓" } else if pct >= 50.0 { "~" } else { "✗" };
            print!(" {:>w$}", format!("{:>3.0}% {}", pct, mark));
            totals[ai] += *h;
        }
        println!();
    }

    println!("{}", "-".repeat(28 + (w + 1) * algos.len()));
    let denom = (SEEDS as u32) * cases().len() as u32;
    print!("{:<22} {:>4} ", "OVERALL", "");
    for t in &totals {
        print!(" {:>w$}", format!("{:>3.0}%", 100.0 * *t as f32 / denom as f32));
    }
    println!("\n\n  ✓ ≥90% correct   ~ ≥50%   ✗ <50%   (over {SEEDS} noise seeds)\n");

    // Optional: BEST_DEBUG=<case substring> dumps why BEST mis-counts that case.
    if let Ok(want) = std::env::var("BEST_DEBUG") {
        for case in cases().into_iter().filter(|c| c.name.contains(&want)) {
            println!("--- debug BEST on '{}' (expected {}) ---", case.name, case.expected);
            for seed in 0..SEEDS {
                let mut rng = Rng::new(seed.wrapping_mul(2654435761) ^ case.name.len() as u64);
                let prof = (case.build)(&mut rng);
                let (peaks, s) = best_peaks(&prof, c, 4.0);
                let got = algo_best(&prof, c, 4.0);
                if got != case.expected {
                    let sigma = noise_sigma(&median3(&prof));
                    let heights: Vec<String> =
                        peaks.iter().map(|&i| format!("@{i}={:.0}", s[i])).collect();
                    println!(
                        "  seed {seed:>2}: got {got} peaks {:?}  σ≈{sigma:.0}  4σ≈{:.0}",
                        heights,
                        4.0 * sigma
                    );
                }
            }
        }
    }
}
