#!/usr/bin/env python3
"""Package both executables and documentation; never upload or publish."""
import argparse
import hashlib
from pathlib import Path
import shutil
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parents[1]
TARGETS = ("x86_64-unknown-linux-gnu", "x86_64-apple-darwin",
           "aarch64-apple-darwin", "x86_64-pc-windows-msvc")


def version():
    return tomllib.loads((ROOT / "Cargo.toml").read_text())["workspace"]["package"]["version"]


def archive_name(target):
    suffix = ".zip" if "windows" in target else ".tar.gz"
    return f"koth_ff-v{version()}-{target}{suffix}"


def verify_dist(directory):
    expected = {archive_name(t) for t in TARGETS}
    expected |= {x + ".sha256" for x in list(expected)}
    actual = {p.name for p in directory.iterdir()}
    if actual != expected:
        raise ValueError(f"Unexpected release artifacts: missing={expected-actual}, extra={actual-expected}")
    for target in TARGETS:
        name = archive_name(target)
        archive = directory / name
        wanted = hashlib.sha256(archive.read_bytes()).hexdigest() + "  " + name
        if not archive.stat().st_size or (directory / (name + ".sha256")).read_text().strip() != wanted:
            raise ValueError(f"Checksum mismatch: {name}")


def package(target, output):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    name = archive_name(target)
    extension = ".exe" if "windows" in target else ""
    with tempfile.TemporaryDirectory() as tmp:
        stage = Path(tmp) / name.removesuffix(".tar.gz").removesuffix(".zip")
        stage.mkdir()
        for binary in ("koth_ff", "koth_align"):
            shutil.copy2(ROOT / "target" / target / "release" / (binary + extension), stage)
        for source in ("README.md", "CHANGELOG.md", "LICENSE", "CITATION.cff", ".zenodo.json",
                       "CONTRIBUTING.md", "RELEASE.md", "example_config.toml", "example_config_align.toml"):
            shutil.copy2(ROOT / source, stage)
        shutil.copytree(ROOT / "docs", stage / "docs")
        fmt = "zip" if extension else "gztar"
        shutil.make_archive(str(output / stage.name), fmt, tmp, stage.name)
    archive = output / name
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    (output / (name + ".sha256")).write_text(f"{digest}  {name}\n")
    print(name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=TARGETS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        verify_dist(args.output)
    elif args.target:
        package(args.target, args.output)
    else:
        parser.error("--target or --verify is required")


if __name__ == "__main__":
    main()
