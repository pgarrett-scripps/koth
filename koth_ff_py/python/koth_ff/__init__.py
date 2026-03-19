"""
koth_ff — high-performance LC-MS feature finder.

The native extension (built by maturin) is imported here so that users can do::

    import koth_ff
    result = koth_ff.run_pipeline("data.mzML")

Functions
---------
detect_hills(path, **kwargs) -> list[dict]
    Stage 1: detect chromatographic hills from an mzML or Bruker .d file.

detect_features(hills, **kwargs) -> list[dict]
    Stages 2–3: group hills into isotope features and score them.

run_pipeline(path, **kwargs) -> dict
    Run all three stages and return {"hills": [...], "features": [...]}.

All functions accept keyword arguments to override any config field.
See the Rust docs for available keys (mz_tolerance, min_scans, max_gap,
intensity_coverage, min_charge, max_charge, …).
"""

from .koth_ff import detect_hills, detect_features, run_pipeline

__all__ = ["detect_hills", "detect_features", "run_pipeline"]
