use std::path::Path;
use std::time::Instant;

use ::koth_ff::config::{
    FeaturesConfig, HillsConfig, ImToleranceType, ScoringConfig, ToleranceType,
};
use ::koth_ff::models::{Hill, ScoredFeature};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

// ---------------------------------------------------------------------------
// Helpers: Rust structs → Python dicts
// ---------------------------------------------------------------------------

fn hill_to_dict<'py>(py: Python<'py>, h: &Hill) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("mz", h.mz)?;
    d.set_item("mz_std", h.mz_std)?;
    d.set_item("rt", h.rt)?;
    d.set_item("rt_start", h.rt_start)?;
    d.set_item("rt_end", h.rt_end)?;
    d.set_item("rt_width", h.rt_width)?;
    d.set_item("im", h.im)?;
    d.set_item("im_std", h.im_std)?;
    d.set_item("scan_start", h.scan_start)?;
    d.set_item("scan_apex", h.scan_apex)?;
    d.set_item("scan_end", h.scan_end)?;
    d.set_item("n_scans", h.n_scans)?;
    d.set_item("skipped_scans", h.skipped_scans)?;
    d.set_item("intensity_sum", h.intensity_sum)?;
    d.set_item("intensity_max", h.intensity_max)?;
    d.set_item("intensity_profile", h.intensity_profile.to_vec())?;
    Ok(d)
}

fn scored_feature_to_dict<'py>(
    py: Python<'py>,
    sf: &ScoredFeature,
) -> PyResult<Bound<'py, PyDict>> {
    let f = &sf.feature;
    let d = PyDict::new(py);
    d.set_item("mz", f.monoisotopic_mz())?;
    d.set_item("mass", f.monoisotopic_neutral_mass())?;
    d.set_item("charge", f.charge)?;
    d.set_item("rt_apex", f.rt_apex())?;
    d.set_item("rt_start", f.rt_start())?;
    d.set_item("rt_end", f.rt_end())?;
    d.set_item("im", f.im_apex())?;
    d.set_item("n_isotopes", f.hills.len())?;
    d.set_item("n_scans", f.n_scans_total())?;
    d.set_item("cosine_similarity", f.cosine_similarity)?;
    d.set_item("ppm_error", f.ppm_error)?;
    d.set_item("score", sf.score)?;
    d.set_item("neutron_offset", sf.neutron_offset)?;
    d.set_item("intensity_sum", f.total_intensity())?;
    d.set_item("intensity_apex", f.total_intensity_at_apex())?;
    d.set_item("theoretical_pattern", sf.theoretical_pattern.clone())?;
    let (_, _, elution) = f.elution_profile();
    d.set_item("elution_profile", elution)?;
    d.set_item("isotope_profile", f.isotope_profile_apex())?;
    Ok(d)
}

// ---------------------------------------------------------------------------
// Config helpers: extract config fields from Python **kwargs dict
// ---------------------------------------------------------------------------

fn hills_config_from_kwargs(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<HillsConfig> {
    let mut cfg = HillsConfig::default();
    let Some(kw) = kwargs else { return Ok(cfg) };

    macro_rules! get {
        ($key:literal, $field:expr, $T:ty) => {
            if let Some(v) = kw.get_item($key)? {
                $field = v.extract::<$T>()?;
            }
        };
    }

    get!("mz_tolerance", cfg.mz_tolerance, f64);
    get!("min_scans", cfg.min_scans, usize);
    get!("max_gap", cfg.max_gap, usize);
    get!("split_hills", cfg.split_hills, bool);
    get!("im_tolerance", cfg.im_tolerance, f64);
    get!("intensity_coverage", cfg.intensity_coverage, f64);
    get!("global_min_mz", cfg.global_min_mz, f64);
    get!("global_max_mz", cfg.global_max_mz, f64);
    get!("bruker_mz_ppm", cfg.bruker_mz_ppm, f64);
    get!("bruker_im_pct", cfg.bruker_im_pct, f64);
    if let Some(v) = kw.get_item("mz_tolerance_type")? {
        let s: String = v.extract()?;
        cfg.mz_tolerance_type = match s.as_str() {
            "ppm" => ToleranceType::Ppm,
            "da" => ToleranceType::Da,
            other => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("unknown mz_tolerance_type '{other}', expected 'ppm' or 'da'")
            )),
        };
    }
    if let Some(v) = kw.get_item("im_tolerance_type")? {
        let s: String = v.extract()?;
        cfg.im_tolerance_type = match s.as_str() {
            "relative" => ImToleranceType::Relative,
            "absolute" => ImToleranceType::Absolute,
            other => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("unknown im_tolerance_type '{other}', expected 'relative' or 'absolute'")
            )),
        };
    }
    Ok(cfg)
}

fn features_config_from_kwargs(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<FeaturesConfig> {
    let mut cfg = FeaturesConfig::default();
    let Some(kw) = kwargs else { return Ok(cfg) };

    macro_rules! get {
        ($key:literal, $field:expr, $T:ty) => {
            if let Some(v) = kw.get_item($key)? {
                $field = v.extract::<$T>()?;
            }
        };
    }

    get!("mz_tolerance", cfg.mz_tolerance, f64);
    get!("min_charge", cfg.min_charge, u8);
    get!("max_charge", cfg.max_charge, u8);
    get!("min_cosine_similarity", cfg.min_cosine_similarity, f64);
    get!("left_max_decrease", cfg.left_max_decrease, f64);
    get!("right_max_decrease", cfg.right_max_decrease, f64);
    get!("im_tolerance", cfg.im_tolerance, f64);
    get!("max_isotopes", cfg.max_isotopes, usize);
    Ok(cfg)
}

fn scoring_config_from_kwargs(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<ScoringConfig> {
    let mut cfg = ScoringConfig::default();
    let Some(kw) = kwargs else { return Ok(cfg) };

    macro_rules! get {
        ($key:literal, $field:expr, $T:ty) => {
            if let Some(v) = kw.get_item($key)? {
                $field = v.extract::<$T>()?;
            }
        };
    }

    get!("isotope_offset_min", cfg.isotope_offset_min, i8);
    get!("isotope_offset_max", cfg.isotope_offset_max, i8);
    get!("offset_zero_bonus", cfg.offset_zero_bonus, f64);
    get!("min_score_threshold", cfg.min_score_threshold, f64);
    Ok(cfg)
}

// ---------------------------------------------------------------------------
// Hill dict → Hill (needed for detect_features round-trip)
// ---------------------------------------------------------------------------

fn dict_to_hill(d: &Bound<'_, PyAny>) -> PyResult<Hill> {
    let d: &Bound<'_, PyDict> = d.cast()?;

    macro_rules! req {
        ($key:literal, $T:ty) => {
            d.get_item($key)?
                .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyKeyError, _>(
                    format!("hill dict missing key '{}'", $key)
                ))?
                .extract::<$T>()?
        };
    }

    let profile_vec: Vec<f32> = req!("intensity_profile", Vec<f32>);

    Ok(Hill {
        mz:               req!("mz", f64),
        mz_std:           req!("mz_std", f64),
        rt:               req!("rt", f64),
        rt_start:         req!("rt_start", f64),
        rt_end:           req!("rt_end", f64),
        rt_width:         req!("rt_width", f64),
        im:               req!("im", f64),
        im_std:           req!("im_std", f64),
        scan_start:       req!("scan_start", usize),
        scan_apex:        req!("scan_apex", usize),
        scan_end:         req!("scan_end", usize),
        n_scans:          req!("n_scans", usize),
        skipped_scans:    req!("skipped_scans", usize),
        intensity_sum:    req!("intensity_sum", f64),
        intensity_max:    req!("intensity_max", f64),
        intensity_profile: std::sync::Arc::from(profile_vec.as_slice()),
    })
}

// ---------------------------------------------------------------------------
// Exposed Python functions
// ---------------------------------------------------------------------------

/// Detect chromatographic hills from an mzML or Bruker .d file.
///
/// Returns a list of dicts; pass to ``pd.DataFrame(hills)`` for a table.
/// Keyword arguments override any :class:`HillsConfig` field, e.g.
/// ``mz_tolerance``, ``min_scans``, ``max_gap``, ``intensity_coverage``.
#[pyfunction]
#[pyo3(signature = (path, **kwargs))]
fn detect_hills(
    py: Python<'_>,
    path: &str,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Py<PyList>> {
    let cfg = hills_config_from_kwargs(kwargs)?;
    let hills = ::koth_ff::run_hills_streaming(Path::new(path), &cfg)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
    let list = PyList::new(
        py,
        hills
            .iter()
            .map(|h| hill_to_dict(py, h))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    Ok(list.unbind())
}

/// Detect isotope features from a list of hill dicts and score them.
///
/// ``hills`` is the list returned by :func:`detect_hills` or the ``"hills"``
/// key of :func:`run_pipeline`.  Returns a list of feature dicts.
/// Keyword arguments override :class:`FeaturesConfig` / :class:`ScoringConfig`
/// fields, e.g. ``min_charge``, ``max_charge``, ``min_cosine_similarity``.
#[pyfunction]
#[pyo3(signature = (hills, **kwargs))]
fn detect_features(
    py: Python<'_>,
    hills: Bound<'_, PyList>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Py<PyList>> {
    let features_cfg = features_config_from_kwargs(kwargs)?;
    let scoring_cfg = scoring_config_from_kwargs(kwargs)?;

    let rust_hills = hills
        .iter()
        .map(|item| dict_to_hill(&item))
        .collect::<PyResult<Vec<_>>>()?;

    let features = ::koth_ff::run_features(&rust_hills, &features_cfg)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
    let scored = ::koth_ff::run_scoring(&features, &scoring_cfg);

    let list = PyList::new(
        py,
        scored
            .iter()
            .map(|sf| scored_feature_to_dict(py, sf))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    Ok(list.unbind())
}

/// Run the complete hill → feature → scoring pipeline in one call.
///
/// Returns a dict with two keys:
///
/// - ``"hills"``    – list of hill dicts
/// - ``"features"`` – list of scored feature dicts
///
/// Example::
///
///     import koth_ff, pandas as pd
///     result      = koth_ff.run_pipeline("data.mzML", mz_tolerance=8.0)
///     hills_df    = pd.DataFrame(result["hills"])
///     features_df = pd.DataFrame(result["features"])
#[pyfunction]
#[pyo3(signature = (path, **kwargs))]
fn run_pipeline(
    py: Python<'_>,
    path: &str,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Py<PyDict>> {
    let hills_cfg = hills_config_from_kwargs(kwargs)?;
    let features_cfg = features_config_from_kwargs(kwargs)?;
    let scoring_cfg = scoring_config_from_kwargs(kwargs)?;

    let t0 = Instant::now();
    let hills = ::koth_ff::run_hills_streaming(Path::new(path), &hills_cfg)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
    let hills_s = t0.elapsed().as_secs_f64();
    eprintln!("[koth_ff] hills: {:.2}s ({} hills)", hills_s, hills.len());

    let t1 = Instant::now();
    let features = ::koth_ff::run_features(&hills, &features_cfg)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
    let features_s = t1.elapsed().as_secs_f64();
    eprintln!("[koth_ff] feature detection: {:.2}s ({} features, {} charge=0)",
        features_s, features.len(),
        features.iter().filter(|f| f.charge == 0).count());

    let t2 = Instant::now();
    let scored = ::koth_ff::run_scoring(&features, &scoring_cfg);
    let scoring_s = t2.elapsed().as_secs_f64();
    eprintln!("[koth_ff] scoring: {:.2}s ({} retained)", scoring_s, scored.len());

    let hills_py = PyList::new(
        py,
        hills
            .iter()
            .map(|h| hill_to_dict(py, h))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    let features_py = PyList::new(
        py,
        scored
            .iter()
            .map(|sf| scored_feature_to_dict(py, sf))
            .collect::<PyResult<Vec<_>>>()?,
    )?;

    let result = PyDict::new(py);
    result.set_item("hills", hills_py)?;
    result.set_item("features", features_py)?;
    Ok(result.unbind())
}

// ---------------------------------------------------------------------------
// Module registration
// ---------------------------------------------------------------------------

#[pymodule]
fn koth_ff(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Initialize env_logger so Rust log::info! output is visible when
    // RUST_LOG=koth_ff=info (or similar) is set in the environment.
    let _ = env_logger::try_init();
    m.add_function(wrap_pyfunction!(detect_hills, m)?)?;
    m.add_function(wrap_pyfunction!(detect_features, m)?)?;
    m.add_function(wrap_pyfunction!(run_pipeline, m)?)?;
    Ok(())
}
