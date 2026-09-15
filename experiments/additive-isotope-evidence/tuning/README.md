# Additive isotope evidence: tuning and held-out validation

The selected model uses best-prefix selection, a signal ratio-error sigma of **0.75 log2 units**, and a **Beta(2, 1)** co-elution model. Selection was frozen on four training files before any held-out evaluation. The noise models remain Normal(0, 2²) for ratio residuals and Uniform(0, 1) for cosine.

## Implementation

Each seed/charge candidate accumulates evidence for every monoisotope-anchored prefix under each allowed sulfur template. It enters the priority queue with its best finite-scoring prefix that passes the existing Bhattacharyya claim gate. After a conflict, it chooses again among the remaining free prefixes. Since the feasible set only shrinks, the new queue score cannot exceed the old score. This fixes the initial pilot’s possibility of hiding strong prefixes behind negative-evidence tails.

The isotope contribution is `ln(2 / 0.75) - 0.5*r²*(1 / 0.75² - 1 / 2²) + ln(2) + ln(c)`, where `r` is the seed-relative log2 apex-ratio error and `c` is seed-anchored chromatographic cosine. The mono seed supplies no evidence against itself. Sulfur selection maximizes the sum for one whole prefix; it never switches templates within a prefix. Strong positive evidence can outweigh an interior penalty, so the search does not simply stop at the first negative term.

Envelope length and the old isotope-priority option do not influence heap ordering. Ties between candidates use mono index and charge; equal-scoring prefixes retain the smaller hill set. Candidate generation, the downstream output scores and thresholds, alignment and LFQ are otherwise unchanged. Negative total evidence is not a new rejection gate: finite negative-score candidates can still be considered if no better candidate claims their hills.

## Predefined experiment

`protocol.json` was written before the sweep. It specifies exactly nine combinations: signal sigma **0.4, 0.5, 0.75** crossed with cosine shape **2, 4, 8**. The noise sigma was fixed at 2. The original four files—A1/A2 (B03_10/B03_11), E1/E2 (B03_05/B03_06)—were used for selection. The four held-out files are A3/A4 (B03_20/B03_21), E3/E4 (B03_15/B03_16). No parameters are changed in response to held-out results.

The sweep reuses the original baseline hill Parquets through `koth_ff/examples/replay_features.rs`. A separate raw-file run at the reference 0.5/4 setting produced a byte-identical feature Parquet to cached-hill replay on B03_10. This checks that replay preserves the production path. All runs use eight threads; the original pilot baseline binary was built from unmodified commit `0908b50`.

Candidates must improve or preserve recall in every training file, increase total features by no more than 5%, and remain within practical tolerances: species median absolute log2 error +0.01, median CV +0.005, fraction with absolute log2 error >0.5 +0.01, and each assignment/chance-match proxy +0.002. These are predefined acceptance tolerances, not significance tests. Among passing candidates, select highest pooled recall, with candidates within 0.05 percentage points favoring proximity to the original 0.5/4 setting.

## Training results

| Ratio sigma | Cosine shape | Recall | Accepted by guardrails |
| ---: | ---: | ---: | --- |
| 0.4 | 2 | 84.3226% | Yes |
| 0.4 | 4 | 84.4033% | wrong_charge_only_rate |
| 0.4 | 8 | 83.8996% | per-file recall, wrong_charge_only_rate |
| 0.5 | 2 | 84.4879% | Yes |
| 0.5 | 4 | 84.4525% | Yes |
| 0.5 | 8 | 83.8189% | per-file recall, wrong_charge_only_rate |
| 0.75 | 2 | 84.5873% | Yes |
| 0.75 | 4 | 84.3079% | Yes |
| 0.75 | 8 | 83.5178% | per-file recall, wrong_charge_only_rate |

The winner increased matched IDs from **85,327 to 85,972**: recall **83.9527% → 84.5873%**, a net **+645** (1,572 gained, 927 lost). Total features decreased from **538,024 to 537,623**. The unchanged 0.5/4 model with the prefix refinement alone reached **84.4525%**; the initial additive pilot without best-prefix selection reached **84.2803%**.

The selected parameters are the best tested combination, not a claimed global optimum. The stronger Beta(8, 1) settings reduced recall and increased the wrong-charge-only proxy. The original tighter parameters were not supported by this bounded comparison.

Quantification uses one fixed cohort of 12372 peptide-charge pairs complete in baseline and all nine variants. Species summaries exclude mixed/other proteins. On this cohort, the selected model leaves E. coli median absolute log2 ratio error unchanged, slightly lowers human error, and changes median CV by less than 0.02 percentage points. The wrong-charge-only proxy increases from 0.8058% to 0.8727% (+0.0669 percentage points), within the predefined 0.2-point limit; the isotope-shift-only proxy decreases from 9.0125% to 8.1368%. The shifted-mass proxy has only 11 baseline versus 9 selected matches and therefore weak statistical power.

## What the checks measure

Anchors are existing Sage target PSMs at peptide-q ≤0.01, deduplicated per run/peptide/charge using the lowest peptide-q and stable tie handling. Correct matches require charge agreement, 10 ppm mass agreement, and PSM RT inside the feature interval. An apex ±0.5-minute check is also saved. Quantification uses the largest matching feature intensity, within-level replicate medians, and expected A/E ratios 1 for human and 1/3 for E. coli; it is raw feature quantification, not a full normalized alignment/LFQ rerun.

Three secondary diagnostics limit obvious assignment regressions. Wrong-charge-only means an anchor lacks a correct match but has a same-m/z/RT feature at another charge. Isotope-shift-only means no correct match but a same-charge feature appears at ±1 or ±2 isotope spacings. The chance-match control shifts every anchor by a fixed non-isotopic neutral mass of 17.331281 Da. These are imperfect proxies: coexisting ions can explain a hit, and unmatched MS/MS anchors do not label all MS1 features. They are not measured feature FDR or confirmed false positives.

The likelihood distributions remain working models, with dependence induced by the shared seed and preselection by generation gates. The tuning/validation split is between runs in the same peptide Orbitrap cohort; it does not establish calibration or performance across instruments or analyte classes.

## Held-out result and decision

**All predefined guardrails passed. Adopt the frozen 0.75/2 model with best-prefix selection.** No parameter was changed after evaluation of the held-out files. This closes the bounded experiment; additional tuning is not required for this change.

| Held-out file | Baseline matches | Selected matches | Recall change (percentage points) |
| --- | ---: | ---: | ---: |
| B03_20 | 20,470 | 20,596 | +0.513 |
| B03_21 | 20,282 | 20,445 | +0.670 |
| B03_15 | 20,967 | 21,114 | +0.583 |
| B03_16 | 20,960 | 21,145 | +0.731 |

| Metric | Baseline | Selected |
| --- | ---: | ---: |
| Matched IDs / 99,435 anchors | 82,679 | 83,300 |
| Pooled interval recall | 83.1488% | 83.7733% |
| Apex-window recall | 80.0754% | 80.5853% |
| Features | 528,314 | 528,151 |
| Features per matched anchor | 1.001826 | 1.001813 |
| HUMAN median absolute log2 error | 0.153069 | 0.152827 |
| HUMAN median CV (%) | 4.668452 | 4.669672 |
| HUMAN absolute log2 error >0.5 (%) | 7.410498 | 7.393808 |
| ECOLI median absolute log2 error | 0.248554 | 0.247969 |
| ECOLI median CV (%) | 4.478586 | 4.477078 |
| ECOLI absolute log2 error >0.5 (%) | 22.776573 | 22.668113 |
| wrong_charge_only_rate (%) | 0.737165 | 0.807563 |
| isotope_shift_only_rate (%) | 9.440338 | 8.547292 |
| shifted_mass_match_rate (%) | 0.014080 | 0.012068 |
| End-to-end wall time, four files | 89.178 s | 88.223 s |

The net **+621** matches comprise **1,465 gained and 844 lost**. Every held-out file improved. Quantification compares 12,910 identical complete peptide-charge pairs, including 11,983 human and 922 E. coli pairs. The wrong-charge-only proxy rose by **0.0704 percentage points** (70 anchors); the isotope-shift-only proxy fell by **0.8930 points**. Thus the result is not an improvement in every diagnostic, but all changes remain inside the predefined limits. Shifted-mass matches fell from 14 to 12; those low counts cannot establish feature FDR.

Timing comes from sequential baseline/selected raw-file runs on each held-out file with eight threads. The selected version was about 1.1% faster in this pass; treat that as essentially unchanged timing rather than a proven speedup. The hill outputs are byte-identical between baseline and selected for all four held-out inputs.

## Reproduction and artifacts

Run `tune.py prepare`, `tune.py training`, `tune.py select`, then `tune.py heldout` from the experiment worktree using the analysis Python environment. `prepare` snapshots the baseline, candidate and replay executables; build the candidate and example with `cargo build --release --locked -p koth_ff --example replay_features --bin koth_ff`. Preparation and per-file execution refuse to overwrite saved outputs. The selected TOML and its SHA-256 are frozen in `selection.json` before held-out execution.

The small reviewable records are `protocol.json`, `training_selection_audit.json`, `selection.json`, the training/held-out summary JSONs, per-file CSVs, and `heldout_decision.json`. Executables, full feature data, matched-anchor Parquets, exact commands, timings, and equivalence checks remain in the ignored `../artifacts/tuning/` directory. They stay in the original worktree at `/home/patrick-garrett/Repos/koth_rust/.worktrees/additive-isotope-evidence/experiments/additive-isotope-evidence/artifacts/` after the source change is merged.

The final defaults are `features.isotope_evidence_ratio_sigma = 0.75` and `features.isotope_evidence_cosine_shape = 2.0`. Old configs pick them up through serde defaults. The legacy `exhaustive_isotope_priority` field remains accepted but no longer changes ranking. The manuscript-pinned build, paper results, Sage searches, alignment and LFQ were not updated.

## Final verification

The production defaults reproduce the frozen selected configuration's held-out
feature Parquet byte for byte on B03_20. Default and no-default workspace tests,
Clippy on all targets with warnings denied, rustdoc with warnings denied, optional
Thermo compilation, formatting and diff checks all pass. Four existing large
Bruker-fixture tests remain explicitly ignored by the default suite.

Regression coverage includes evidence signs, clean extensions, contaminated
envelopes, initial best-prefix selection, monotone conflict rescoring across
sulfur templates, claim-gate fallback to valid alternatives, interior penalties,
deterministic ties, invalid numerical inputs/configs, and actual resolver hill
ownership. Detailed checks and equivalence hashes are saved in `validation.json`.
