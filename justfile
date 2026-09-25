default:
    @just --list

# Build in debug mode (all workspace members)
build:
    cargo build --workspace

# Build optimized release binary
release:
    cargo build --release -p koth-ms

# Check for compile errors without building
check:
    cargo check --workspace

# Copy Cargo.toml's version and a release date (default today) into CITATION.cff
cite-sync DATE=`date +%F`:
    #!/usr/bin/env bash
    set -euo pipefail
    v=$(sed -n 's/^version *= *"\(.*\)"/\1/p' Cargo.toml | head -1)
    sed -i -E "s/^version: .*/version: \"$v\"/; s/^date-released: .*/date-released: \"{{DATE}}\"/" CITATION.cff
    grep -E '^(version|date-released):' CITATION.cff

# Fail if CITATION.cff's version or date-released disagrees with Cargo.toml
cite-check:
    uv run --no-project --with-requirements scripts/release-requirements.txt python scripts/check_release.py

# Run all tests
test:
    cargo test --workspace

# Run on an mzML file (usage: just run-mzml path/to/file.mzML)
run-mzml FILE OUTPUT="./out":
    cargo run --release -p koth-ms -- {{FILE}} --output {{OUTPUT}}

# Run on a Bruker .d folder (usage: just run-bruker path/to/data.d)
run-bruker FILE OUTPUT="./out":
    cargo run --release -p koth-ms -- {{FILE}} --output {{OUTPUT}}

# Run with a custom config TOML
run-config FILE CONFIG OUTPUT="./out":
    cargo run --release -p koth-ms -- {{FILE}} --config {{CONFIG}} --output {{OUTPUT}}

# Run with debug logging
run-debug FILE OUTPUT="./out":
    cargo run --release -p koth-ms -- {{FILE}} --output {{OUTPUT}} --log-level debug

# Skip scoring stage (faster, no isotope pattern scoring)
run-no-score FILE OUTPUT="./out":
    cargo run --release -p koth-ms -- {{FILE}} --output {{OUTPUT}} --no-scoring

# Quick test against the zenith_feature_finder mzML sample
test-sample OUTPUT="/tmp/koth_test":
    cargo run --release -p koth-ms -- \
        ../zenith_feature_finder/20250806_ArgC_DDA_HCD-FT_01.mzML \
        --output {{OUTPUT}} \
        --log-level info

# Print the first few lines of hills output
peek-hills DIR="./out":
    head -3 {{DIR}}/*/hills.tsv

# Print the first few lines of features output
peek-features DIR="./out":
    head -3 {{DIR}}/*/features.tsv

# Count output rows
count DIR="./out":
    @echo "Hills: $(tail -n +2 {{DIR}}/*/hills.tsv | wc -l)"
    @echo "Features: $(tail -n +2 {{DIR}}/*/features.tsv | wc -l)"

# Clean build artifacts
clean:
    cargo clean

# Format code
fmt:
    cargo fmt --all

# Fail if koth-core gained a file-I/O, CLI or native dependency
core-deps:
    .github/scripts/check-core-deps.sh

# Run clippy linter
lint:
    cargo clippy --workspace -- -D warnings

# Show binary size
size:
    ls -lh target/release/koth_ff 2>/dev/null || echo "Run 'just release' first"
