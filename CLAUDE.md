# koth_rust

## What this is

Two Rust binaries built together (`cargo build --release`):

- **`koth_ff`** — per-run MS1 feature finder. `raw → [file] read/tolerances →
  [hills] chromatographic traces → [features] isotope chains (charge + averagine
  score) → [scoring] optional → hills + features output`. mzML, Bruker `.d`
  (needs the default `tdf` feature, backed by the published `dnoise` crate), and optional Thermo
  `.raw` (`thermo` feature).
- **`koth_align`** — cross-run alignment + LFQ. `koth_ff run dirs → [alignment]
  RANSAC RT warp + mass/IM drift → [lfq.consensus] group the same peptide across
  runs → [lfq] XIC grid + integrate + target-decoy q-values → intensity matrix`.

Source layout: `koth_ff/src/{config,hills,features,scoring,alignment,lfq,io,output}/`.
The active q-value path is the QDA rescorer in `lfq/rescore.rs` (not the simpler
ranker in `lfq/tdc.rs`).

## Configuration

**[`docs/CONFIGURATION.md`](docs/CONFIGURATION.md) is the complete, authoritative
reference for every config setting** — type, default, what it does, what to set
it to, and Orbitrap-vs-Bruker platform notes. Read it before touching config or
answering config questions; do not guess field names or defaults from memory.

- Configs are TOML parsed with `#[serde(deny_unknown_fields)]` — a stale/renamed
  key is a hard build/parse error. A build-time test parses this repo's
  `example_config*.toml` against the structs, so field names can't silently
  drift. If you add/rename a config field, update the struct, both templates,
  **and** `docs/CONFIGURATION.md` (the doc is not build-checked).
- The benchmark and the manuscript now live in a separate repository,
  `tacular-omics/koth-paper`. Its `analysis/config/*.toml` are the tuned
  per-platform configs; they are NOT build-checked from here, so a renamed field
  breaks them at benchmark run time rather than at compile time. Update them
  there when you rename a config field.
- Experimental knobs are **default-off and byte-safe** by convention. Their
  tested verdicts (validated-win / dud / neutral) are recorded in
  `docs/CONFIGURATION.md` §5 — check there before re-running a settled
  experiment. Notably `averagine_projection` is a measured dud on Orbitrap.

## Testing conventions

Large unit-test modules live in a **sibling file** next to the module they test,
wired in via `#[path]` rather than inlined at the bottom of the source file:

```rust
// in foo.rs
#[cfg(test)]
#[path = "foo_tests.rs"]
mod tests;
```

- A file named `<stem>_tests.rs` (e.g. `assemble_tests.rs`, `warp_tests.rs`) holds
  the unit tests for its same-named production module (`assemble.rs`, `warp.rs`).
  Because `#[path]` makes it a **child module** of the production module, its tests
  see all private items exactly as an inline `mod tests` would (they start with
  `use super::*;`). These are unit tests, **not** integration tests — do not move
  them to `tests/` and do not "helpfully" re-inline them back into the source file.
- Only large test blocks (roughly >= 80 lines) are split out this way; small test
  modules are kept inline. When adding tests to a module that already has a
  `<stem>_tests.rs`, put them in that sibling file.
