#!/usr/bin/env bash
# Fail if koth-core's normal (non-dev, non-build) dependency tree pulls in file
# I/O, a CLI, or native code: mzdata, timsrust*, rusqlite, libsqlite3-sys,
# arrow*, parquet, csv, flate2, opentfraw, dnoise*, clap, env_logger, or any
# crate with a `links` key other than rayon-core (pure Rust; its key only stops
# two rayon-core major versions being linked together). sage-plus and other
# embedders link koth-core next to their own versions of those crates.
set -euo pipefail
cd "$(dirname "$0")/../.."

tree=$(cargo tree -p koth-core -e normal --prefix none --locked)
echo "$tree"

forbidden='^(mzdata|timsrust[a-z0-9_-]*|rusqlite|libsqlite3-sys|arrow[a-z0-9_-]*|parquet|csv|flate2|opentfraw|dnoise[a-z0-9_-]*|clap|env_logger|zstd[a-z0-9_-]*) '
if echo "$tree" | grep -Eq "$forbidden"; then
  echo "error: forbidden dependency in koth-core:" >&2
  echo "$tree" | grep -E "$forbidden" >&2
  exit 1
fi

# Packages in the tree (name + version) that declare a `links` key.
deps_file=$(mktemp)
trap 'rm -f "$deps_file"' EXIT
echo "$tree" | awk '{print $1, $2}' | sort -u > "$deps_file"
cargo metadata --format-version 1 --locked | python3 -c '
import json, sys
allowed = {"rayon-core"}
deps = {tuple(l.split()) for l in open(sys.argv[1]) if l.strip()}
meta = json.load(sys.stdin)
bad = sorted(p["name"] for p in meta["packages"]
             if (p["name"], "v" + p["version"]) in deps and p.get("links")
             and p["name"] not in allowed)
if bad:
    sys.exit("error: koth-core depends on crates with a links key: " + ", ".join(bad))
print(f"koth-core: {len(deps)} normal dependencies, none forbidden, "
      "no links key except rayon-core")
' "$deps_file"
