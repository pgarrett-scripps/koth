# Releasing `koth-ff`

The Python package (`koth-ff`) is built from `koth_ff_py/` with maturin and
published to PyPI by `.github/workflows/publish.yml` whenever a `v*` tag is
pushed. The CLI binary (`koth_ff`) is published to GitHub Releases by
`.github/workflows/release.yml` on the same tag.

## Pre-release checklist

1. **Bump the version in two places** (they must match exactly):
   - `Cargo.toml` → `[workspace.package].version`
   - `koth_ff_py/pyproject.toml` → `[project].version`

   The publish workflow has a tag-vs-pyproject guard and will fail-fast if
   these disagree with the tag, so a mismatch is loud rather than silent.

2. **Refresh `Cargo.lock`** so the version bump is committed atomically:
   ```bash
   cargo check --workspace
   ```

3. **Update `CHANGELOG.md`** with a new section for the version.

4. **Commit and push** the bump on a normal feature branch / PR. Don't tag
   from the bump commit until it's merged to `main`.

## Cutting the release

1. Make sure `main` is at the commit you want to release.
2. Create and push the tag:
   ```bash
   git checkout main && git pull
   git tag v0.1.0
   git push origin v0.1.0
   ```
3. Watch the `Publish Python wheels to PyPI` workflow on Actions. It will:
   - run `cargo test --workspace`
   - build wheels for Linux x86_64 / aarch64, macOS x86_64 / aarch64, Windows
     x86_64 (single `cp310-abi3` wheel per platform — covers Python 3.10+)
   - build an sdist
   - verify the tag matches `koth_ff_py/pyproject.toml`
   - upload to PyPI via OIDC trusted publishing

4. Verify on PyPI: <https://pypi.org/project/koth-ff/>

## Smoke test the published package

```bash
uv venv /tmp/koth-smoke && source /tmp/koth-smoke/bin/activate
pip install koth-ff==0.1.0
python -c "import koth_ff; print(koth_ff.__all__)"
```

## If something goes wrong

- **Tag pushed but PyPI rejected the upload** (e.g. version already exists):
  delete the tag, bump the version, push the new tag. PyPI never lets you
  re-upload the same `name==version` — that's the trap that the tag-version
  guard exists to prevent.
  ```bash
  git tag -d v0.1.0
  git push --delete origin v0.1.0
  ```

- **aarch64 Linux wheel build fails**: timsrust pulls in a vendored sqlite C
  build via `rusqlite` `bundled`. If the cross-compile fails, the fallback is
  to make `tdf` target-conditional in `koth_ff_py/Cargo.toml` (see the
  release prep plan for the snippet) so the aarch64 Linux wheel ships without
  Bruker timsTOF support.

- **Wheel test in CI passes but `import koth_ff` fails on a user's machine**:
  check that they're on Python ≥ 3.10 (wheel ABI floor) and that pip resolved
  to the abi3 wheel rather than the sdist.
