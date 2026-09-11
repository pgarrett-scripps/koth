# Release checklist

The tag workflow publishes a GitHub release automatically. Do not create or push
a `v*` tag until every manual item below is complete.

## Automated gates

Run these from the repository root. CI runs the same checks on Rust 1.88.0.

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked
cargo test --workspace --locked
cargo test --workspace --locked --no-default-features
cargo package -p koth_ff --locked
cargo build --release --locked -p koth_ff --bins
target/release/koth_ff --version
target/release/koth_align --version
```

Confirm the GitHub Actions matrix passes on Linux, macOS, and Windows. The
release workflow will refuse a tag whose name is not exactly `v<version from
Cargo.toml>`, and it only runs for a tag that has been **pushed** — a local tag
publishes nothing (0.3.0 sat unpushed for three days).

## Manual gates

- Decide the release version and set it in the workspace `Cargo.toml`. Pre-1.0,
  a new config key or a changed default is a minor bump; the released line so far
  is 0.1.0 (2026-08-18), 0.2.0 (2026-09-03), 0.3.0 (2026-09-08).
- Update both example TOMLs and `docs/CONFIGURATION.md` for every config change.
- Move every `## [Unreleased]` entry into the new dated version section. Entries
  left there describe code that is already tagged and read as unreleased.
- `git push origin master` **and** `git push origin vX.Y.Z`, then confirm the
  release workflow ran and the GitHub release exists.
- Resolve or explicitly accept the changelog warning that absolute QDA q-value
  calibration has not yet been independently validated.
- Run the four ignored Bruker real-data integration tests with their fixtures.
- Build and smoke-test `--features thermo` on a machine with .NET 8 and a native
  Thermo `.raw` fixture.
- Confirm the tuned configs in the separate `koth-paper` repository still parse
  and reproduce the intended benchmark outputs.
- Confirm the repository is public if the release is intended to be public, and
  configure branch protection so CI is required on `master`.
- Decide whether to publish `koth_ff` to crates.io in addition to GitHub binaries.
  `cargo package` is already gated, but crate publication is deliberately not
  automated.
- Review generated release notes and archive metadata/citation text after the
  paper citation is final.

## Release

Only after all gates pass:

```bash
git tag -s vX.Y.Z -m "koth_ff vX.Y.Z"
git push origin vX.Y.Z
```

The workflow builds both binaries for Linux x86_64, macOS x86_64 and arm64, and
Windows x86_64; packages README, changelog, and license files; emits SHA-256
checksums; and creates one GitHub release after every target succeeds.
