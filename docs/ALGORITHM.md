# Feature-detection algorithm

The single-run pipeline runs in three sequential stages. Multi-run analysis is
described in [alignment and LFQ](LFQ.md).

### Stage 1: Hill detection

Hills are chromatographic traces — a single m/z signal tracked across consecutive MS1 scans.

The detector consumes MS1 spectra incrementally for mzML, Bruker `.d`, and Thermo `.raw`
(`io::stream_spectra`). Bruker decodes batches of at most 32 frames in parallel and queues
at most 32 spectra; mzML and Thermo prefetch at most 64 spectra. Native readers sort scan
metadata before decoding to preserve retention-time order without retaining all peak arrays.

For each scan the detector builds a scratch buffer of `(mz_mean, hill_id)` pairs sorted by m/z,
then binary-searches it for each incoming peak. Because peaks arrive sorted by m/z within a scan,
each binary search lands in a narrow window and the inner loop is short. Compared to a HashMap bin
approach there is no hash overhead and no bin-size tuning.

Matching rules:

- Closest unmatched hill within `mz_tolerance` (ppm or Da) is selected, with a combined m/z +
  intensity log-fold-change distance score weighted by `lfc_weight`.
- If ion mobility data is present, an IM tolerance check is applied and the combined
  m/z + IM distance is used to break ties. Bruker `.d` mobility is converted to
  1/K0 with the run's acquisition calibration (`TimsCalibration`, ModelType 2)
  before hills are built; `bruker_mobility_scale = "linear"` restores the
  timsrust linear converter used up to 0.9.0.
- A matched peak extends the hill's running m/z mean (Welford online update) and appends
  the intensity to the profile.
- An unmatched peak starts a new hill.

After each peak is processed the detector periodically evicts stale hills — any hill whose
`last_scan_seen` is more than `max_gap` scans behind the current scan index. Evicted hills
shorter than `min_scans` are discarded; the rest are finalized and appended to the output list.

Intensity profiles are stored as `f32` vectors (4 bytes per position). Gap positions — scans
where no peak matched — are represented as `0.0` intensity rather than `Option<f64>`, cutting
per-element memory from 16 bytes to 4 bytes.

### Co-elution splitting

After all hills are finalized, any hill that contains multiple local intensity maxima is split
into sub-hills at the valley between them.

Split criteria (all must hold):

1. Two local maxima are at least `min_peak_distance` scans apart.
2. Each maximum reaches at least `min_peak_height` fraction of the hill's global maximum.
3. Each maximum has a prominence (peak height minus the highest valley between it and any taller
   neighbour) of at least `min_prominence` fraction of the global maximum. This prevents
   noise wiggles on a flank from triggering false splits.
4. Both resulting segments are at least `min_scans` scans long.

Splitting is controlled by `split_hills = true` in the config.

### Stage 2: Feature detection

Features are isotope envelopes — groups of hills whose m/z values are spaced by
`neutron_mass / charge` (where `neutron_mass` = 1.003354835 Da, the C13 offset).

Every hill is tried as a monoisotopic seed, at every charge state from `min_charge` to
`max_charge`. Chains extend **upward only** — seed → M+1 → M+2 → … — because the seed *is*
the monoisotopic hypothesis. There is no downward walk: a hill that one could reach is
itself a seed that builds the same envelope upward, with the averagine template indexed
from its own position.

A candidate partner must satisfy:

- m/z within `mz_tolerance` of the expected isotope position.
- Scan range overlaps with the reference hill.
- If IM data is present, IM within `im_tolerance` of the reference hill.
- `intensity_max >= ref_hill.intensity_max * min_isotope_step_ratio` (prevents linking to an
  implausibly weak signal).
- Cosine of the elution profile against the `cosine_anchor` reference hill (the seed by
  default) is at least `min_chain_cosine`.
- The apex intensity ratio against the chain predecessor matches the averagine ratio within
  ±`max_isotope_log2_ratio`.

This produces an over-complete `(seed, charge)` candidate pool, which is then resolved
**non-destructively**: each candidate competes using its best-scoring valid
monoisotope-anchored prefix. The score sums per-isotope log evidence from the
seed-relative averagine intensity ratio and chromatographic co-elution against
explicit noise models. Good isotopes add evidence and poor ones subtract;
envelope length never overrides the score. After a conflict, the resolver selects
the best remaining free prefix and requeues it, so its priority cannot increase.
Feature detection also records per-adjacent-pair cosine similarities and ppm
errors for downstream quality filtering.

### Stage 3: Scoring

The default peptide model scales an average residue composition to the inferred
neutral mass. Cached per-element isotope distributions are convolved to form
the theoretical envelope, with configurable sulfur-count alternatives for
peptides. The model is not a fixed 50-Da mass lookup table. See the
[configuration reference](CONFIGURATION.md) for model and sulfur settings.

Scoring compares the observed isotope profile with the theoretical pattern
using a Bhattacharyya coefficient and a penalty for missing theoretical
intensity. If `isotope_offset_enabled` is true, offsets from −1 through +1 are
tested; otherwise only zero is considered. `offset_zero_bonus` favors retaining
the original assignment when scores are close, and
`min_isotope_score_for_offset` controls the fallback to zero offset.

The output separates `isotope_score`, `cosine_score`, and `combined_score`.
These are pattern and co-elution quality measures, not identification
probabilities. Consult [output columns](OUTPUTS.md) for their exported names.

## Memory design

With default settings, MS1 hill detection uses bounded spectrum buffers for every input format.
It updates active hills and discards each consumed spectrum, retaining completed hills for
isotope assembly. Native readers also retain scan metadata and vendor/library state; total
process memory therefore includes more than the spectrum buffers. The `bruker_streaming`
setting selects the dnoise preprocessing mode; both modes stream spectra into detection.

Explicit decoy shuffling and optional TIC normalization still collect spectra. Gzipped mzML
is decompressed as it is parsed, so peak memory does not include the decompressed file;
multi-member (bgzip-style) files are read to the end. The collecting
`read_spectra` API remains available for callers that need a complete spectrum vector.

Once hill detection is complete the hill list is passed to feature detection. Hills hold their
intensity profiles as `Arc<[f32]>`, so features can reference the same profile data without
copying. After feature detection the hill `Vec` is explicitly dropped; the `Arc` reference
counts keep the profiles alive inside each `Feature`.

`koth_align` applies the same discipline across runs: it holds only the per-run *features* (small)
in memory for the whole run, and **streams each run's hills one at a time** during LFQ — loading a
run's hills, quantifying every consensus feature against them, then dropping them before the next
run loads. Hill-profile memory is bounded by one run at a time, while the feature,
consensus, and result structures still grow with the cohort size.
