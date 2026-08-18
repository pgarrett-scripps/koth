//! Isotope-consistency m/z recalibration (Biosaur-style, ID-free).
//!
//! Self-calibration of the m/z axis learned from the *isotope spacing* of
//! assembled features. For a feature of charge `z`, adjacent isotope hills
//! should be separated by exactly `neutron_mass / z`. The signed deviation of
//! the observed spacing from that theoretical step is a sample of the local,
//! systematic m/z error — expressed in ppm and binned over (m/z, RT).
//!
//! A pass-1 feature detection collects these residuals; the resulting median
//! offset surface is applied in pass 2 to shift the *expected* isotope
//! position during chain extension, so isotope hills are searched for at their
//! recalibrated location rather than the naive grid position.
//!
//! ## What this can and cannot correct
//!
//! Isotope spacing constrains only the *proportional* (ppm) component of the
//! m/z error, and only its per-region structure. A hypothetical constant Da
//! offset shared by every ion is unobservable without an external mass
//! reference (a lock mass or peptide identifications), neither of which exists
//! at feature-finding time. In practice the proportional, m/z- and
//! RT-dependent component is the dominant Orbitrap error mode, and that is
//! exactly what this recovers. The surface is deterministic (robust medians,
//! fixed binning) and read immutably by the parallel feature pass.

/// A 2-D (m/z × RT) surface of signed median ppm offsets, with hierarchical
/// fallback: cell median → m/z-marginal median → global median. Predicting
/// outside the trained range clamps to the edge bins.
#[derive(Debug, Clone)]
pub struct MzRecalModel {
    mz_min: f64,
    mz_inv_width: f64,
    n_mz: usize,
    rt_min: f64,
    rt_inv_width: f64,
    n_rt: usize,
    /// Per-cell median ppm (row-major, `mz * n_rt + rt`); `None` if the cell
    /// held fewer than `min_samples` residuals.
    cell: Vec<Option<f64>>,
    /// Per-m/z-bin marginal median ppm; `None` if that column was empty.
    mz_marginal: Vec<Option<f64>>,
    global: f64,
    /// Per-cell robust (MAD) spread ppm; `None` below `min_samples`. Parallel
    /// to `cell`. Used by the region-adaptive isotope-match tolerance.
    cell_sigma: Vec<Option<f64>>,
    /// Per-m/z-bin marginal spread ppm; `None` if empty.
    mz_marginal_sigma: Vec<Option<f64>>,
    /// Total residual samples that fed the surface (diagnostics).
    pub n_samples: u64,
    /// Robust (MAD-based) global spread of the residuals, in ppm (diagnostics
    /// and the final fallback for [`predict_sigma`](Self::predict_sigma)).
    pub global_sigma: f64,
}

impl MzRecalModel {
    /// Signed ppm offset predicted for a peak at `(mz, rt)`. Positive means the
    /// observed m/z runs *high* of the true value, so callers shift the
    /// expected isotope position by `(1 + predict/1e6)`.
    #[inline]
    pub fn predict(&self, mz: f64, rt: f64) -> f64 {
        let mi = self.mz_bin(mz);
        let ri = self.rt_bin(rt);
        if let Some(v) = self.cell[mi * self.n_rt + ri] {
            return v;
        }
        if let Some(v) = self.mz_marginal[mi] {
            return v;
        }
        self.global
    }

    /// Global median offset (ppm) — the fallback value and a useful summary.
    pub fn global_offset(&self) -> f64 {
        self.global
    }

    /// Robust spread (ppm) predicted for `(mz, rt)`, same hierarchical fallback
    /// as [`predict`](Self::predict): cell σ → m/z-marginal σ → global σ. Used
    /// to set a region-adaptive isotope-match tolerance.
    #[inline]
    pub fn predict_sigma(&self, mz: f64, rt: f64) -> f64 {
        let mi = self.mz_bin(mz);
        let ri = self.rt_bin(rt);
        if let Some(v) = self.cell_sigma[mi * self.n_rt + ri] {
            return v;
        }
        if let Some(v) = self.mz_marginal_sigma[mi] {
            return v;
        }
        self.global_sigma
    }

    #[inline]
    fn mz_bin(&self, mz: f64) -> usize {
        let b = ((mz - self.mz_min) * self.mz_inv_width) as isize;
        b.clamp(0, self.n_mz as isize - 1) as usize
    }

    #[inline]
    fn rt_bin(&self, rt: f64) -> usize {
        let b = ((rt - self.rt_min) * self.rt_inv_width) as isize;
        b.clamp(0, self.n_rt as isize - 1) as usize
    }
}

/// Accumulates `(m/z, RT, residual-ppm)` samples during a pass-1 feature pass,
/// then bins them into an [`MzRecalModel`] on [`finalize`](Self::finalize).
#[derive(Debug)]
pub struct MzRecalBuilder {
    samples: Vec<(f64, f64, f64)>,
    max_abs_ppm: f64,
}

impl MzRecalBuilder {
    /// `max_abs_ppm` hard-caps accepted residuals so a mis-assembled feature
    /// (e.g. a spacing that skipped an isotope) can't poison the medians. The
    /// medians are robust anyway; this mainly protects the marginal fallback.
    pub fn new(max_abs_ppm: f64) -> Self {
        Self {
            samples: Vec::new(),
            max_abs_ppm: max_abs_ppm.max(0.0),
        }
    }

    #[inline]
    pub fn add(&mut self, mz: f64, rt: f64, resid_ppm: f64) {
        if !mz.is_finite() || !rt.is_finite() || !resid_ppm.is_finite() {
            return;
        }
        if resid_ppm.abs() > self.max_abs_ppm {
            return;
        }
        self.samples.push((mz, rt, resid_ppm));
    }

    pub fn len(&self) -> usize {
        self.samples.len()
    }

    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Build the surface. Returns `None` if fewer than `min_samples` (min 1)
    /// residuals were collected — not enough to be statistically meaningful.
    /// Bin extents are derived from the observed sample range, so a run that
    /// only populated part of the m/z axis doesn't waste bins on empty space.
    pub fn finalize(self, n_mz: usize, n_rt: usize, min_samples: usize) -> Option<MzRecalModel> {
        let n_mz = n_mz.max(1);
        let n_rt = n_rt.max(1);
        let need = min_samples.max(1);
        if self.samples.len() < need {
            return None;
        }

        let (mut mz_lo, mut mz_hi) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut rt_lo, mut rt_hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for &(mz, rt, _) in &self.samples {
            mz_lo = mz_lo.min(mz);
            mz_hi = mz_hi.max(mz);
            rt_lo = rt_lo.min(rt);
            rt_hi = rt_hi.max(rt);
        }
        // Guard against a degenerate (zero-span) axis: widen by an epsilon so
        // the inverse width is finite and every sample lands in bin 0.
        let mz_span = (mz_hi - mz_lo).max(1e-9);
        let rt_span = (rt_hi - rt_lo).max(1e-9);
        // Nudge the upper edge so the max sample maps inside the last bin
        // rather than one past it (before clamping).
        let mz_inv_width = n_mz as f64 / (mz_span * (1.0 + 1e-9));
        let rt_inv_width = n_rt as f64 / (rt_span * (1.0 + 1e-9));

        let mut cell_buckets: Vec<Vec<f64>> = vec![Vec::new(); n_mz * n_rt];
        let mut mz_buckets: Vec<Vec<f64>> = vec![Vec::new(); n_mz];
        let mut all: Vec<f64> = Vec::with_capacity(self.samples.len());

        let mz_bin = |mz: f64| -> usize {
            (((mz - mz_lo) * mz_inv_width) as isize).clamp(0, n_mz as isize - 1) as usize
        };
        let rt_bin = |rt: f64| -> usize {
            (((rt - rt_lo) * rt_inv_width) as isize).clamp(0, n_rt as isize - 1) as usize
        };

        for &(mz, rt, resid) in &self.samples {
            let mi = mz_bin(mz);
            let ri = rt_bin(rt);
            cell_buckets[mi * n_rt + ri].push(resid);
            mz_buckets[mi].push(resid);
            all.push(resid);
        }

        // Median + MAD σ per bucket in one pass. Cells below `need` samples
        // yield `None` (→ fallback); marginals need only be non-empty.
        let mut cell = Vec::with_capacity(cell_buckets.len());
        let mut cell_sigma = Vec::with_capacity(cell_buckets.len());
        for mut v in cell_buckets {
            if v.len() >= need {
                let m = median(&mut v);
                cell.push(Some(m));
                cell_sigma.push(Some(mad_sigma(&mut v, m)));
            } else {
                cell.push(None);
                cell_sigma.push(None);
            }
        }
        let mut mz_marginal = Vec::with_capacity(mz_buckets.len());
        let mut mz_marginal_sigma = Vec::with_capacity(mz_buckets.len());
        for mut v in mz_buckets {
            if v.is_empty() {
                mz_marginal.push(None);
                mz_marginal_sigma.push(None);
            } else {
                let m = median(&mut v);
                mz_marginal.push(Some(m));
                mz_marginal_sigma.push(Some(mad_sigma(&mut v, m)));
            }
        }

        let global = median(&mut all);
        let global_sigma = mad_sigma(&mut all, global);

        Some(MzRecalModel {
            mz_min: mz_lo,
            mz_inv_width,
            n_mz,
            rt_min: rt_lo,
            rt_inv_width,
            n_rt,
            cell,
            mz_marginal,
            global,
            cell_sigma,
            mz_marginal_sigma,
            n_samples: self.samples.len() as u64,
            global_sigma,
        })
    }
}

/// Median of `v` (mutates by sorting). `v` must be non-empty.
fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.total_cmp(b));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// MAD-scaled robust σ estimate (1.4826 × median absolute deviation).
fn mad_sigma(v: &mut [f64], center: f64) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    let mut dev: Vec<f64> = v.iter().map(|x| (x - center).abs()).collect();
    1.4826 * median(&mut dev)
}

#[cfg(test)]
#[path = "recalibration_tests.rs"]
mod tests;
