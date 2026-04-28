import koth_ff

MZML = "tests/data/example_dda.d"

result = koth_ff.run_pipeline(
    MZML,
    # hills
    mz_tolerance=5.0,
    mz_tolerance_type="ppm",
    min_scans=3,
    max_gap=0,
    split_hills=False,
    im_tolerance=0.05,
    im_tolerance_type="relative",
    intensity_coverage=0.99,
    # features
    min_charge=1,
    max_charge=7,
    min_cosine_similarity=0.5,
    max_isotopes=6,
    # scoring
    min_score_threshold=0.5,
    offset_zero_bonus=0.15,
)

hills    = result["hills"]
features = result["features"]

print(f"\n=== Hills ({len(hills):,} rows) ===")
print(hills.head(10))

print(f"\n=== Features ({len(features):,} rows) ===")
print(features.head(10))
