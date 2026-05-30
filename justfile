default:
    @just --list

# Build in debug mode (all workspace members)
build:
    cargo build --workspace

# Build optimized release binary
release:
    cargo build --release -p koth_ff

# Check for compile errors without building
check:
    cargo check --workspace

# Run all tests
test:
    cargo test --workspace

# Run on an mzML file (usage: just run-mzml path/to/file.mzML)
run-mzml FILE OUTPUT="./out":
    cargo run --release -p koth_ff -- {{FILE}} --output {{OUTPUT}}

# Run on a Bruker .d folder (usage: just run-bruker path/to/data.d)
run-bruker FILE OUTPUT="./out":
    cargo run --release -p koth_ff -- {{FILE}} --output {{OUTPUT}}

# Run with a custom config TOML
run-config FILE CONFIG OUTPUT="./out":
    cargo run --release -p koth_ff -- {{FILE}} --config {{CONFIG}} --output {{OUTPUT}}

# Run with debug logging
run-debug FILE OUTPUT="./out":
    cargo run --release -p koth_ff -- {{FILE}} --output {{OUTPUT}} --log-level debug

# Skip scoring stage (faster, no isotope pattern scoring)
run-no-score FILE OUTPUT="./out":
    cargo run --release -p koth_ff -- {{FILE}} --output {{OUTPUT}} --no-scoring

# Quick test against the zenith_feature_finder mzML sample
test-sample OUTPUT="/tmp/koth_test":
    cargo run --release -p koth_ff -- \
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

# Run clippy linter
lint:
    cargo clippy --workspace -- -D warnings

# Show binary size
size:
    ls -lh target/release/koth_ff 2>/dev/null || echo "Run 'just release' first"
