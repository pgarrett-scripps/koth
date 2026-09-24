# Contributing

Report bugs and request features through [GitHub Issues](https://github.com/pgarrett-scripps/koth/issues).
Include the tool version, input format, command, relevant configuration, and
expected versus observed behavior. Use a small shareable example where possible.

For changes, open a pull request against `master`. Describe the problem and
validation. Update the configuration reference and example TOMLs when changing
settings, and add an entry under `Unreleased` in the changelog.

Run the checks used by CI:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked
cargo test --workspace --locked
cargo test --workspace --locked --no-default-features
```

Install `scripts/release-requirements.txt` in a Python environment and run
`python scripts/check_release.py` for documentation and metadata validation.
See [development](docs/DEVELOPMENT.md) for library examples and
[RELEASE.md](RELEASE.md) for publishing setup. Normal pushes and pull requests
never publish packages or create releases.
