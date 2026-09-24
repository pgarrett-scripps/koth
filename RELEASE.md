# Releases

Publishing a GitHub Release is the only publication trigger. A branch push or
tag push does not publish binaries or a crate. Manual runs of the Release
workflow perform CI and packaging without publishing anything.

## One-time setup

1. Add a crates.io API token as the repository Actions secret
   `CARGO_REGISTRY_TOKEN`. Scope it to publishing `koth-ms`; for the first
   publication it must also permit creation of that crate. No token belongs in
   source control. This follows the existing dnoise release setup and supports
   the first crate publication without requiring an already registered crate.
2. When ready for public access, make the repository public and connect it in
   [Zenodo's GitHub settings](https://zenodo.org/account/settings/github/).
   Enable `pgarrett-scripps/koth` before publishing the release. The repository
   remains private during preparation; the workflow refuses public package
   publication from a private repository.
3. Configure required CI checks on `master` when repository/account settings
   allow it. Run the Release workflow manually to test all target builds and
   the crate dry run before publishing a release.

Zenodo's native GitHub integration archives source from the release and mints
a DOI. It is an independent listener to the same published-release event:
**it does not wait for GitHub Actions or crates.io publication**. CI gates the
binaries and crate; a failed CI run does not undo a Zenodo archive or the GitHub
release. Run the manual checks before publishing, and check Zenodo's archival
status afterward. This integration does not promise to include binary assets
attached by a later Actions job. Do not add a second Zenodo uploader, which
could create duplicate records.

## Prepare a release

1. Choose the version in the workspace `Cargo.toml`, refresh `Cargo.lock`,
   run `just cite-sync [YYYY-MM-DD]` (date defaults to today) to write that
   version and the release date into `CITATION.cff`, move applicable changelog
   entries into a dated release section, and update changed configuration
   examples and documentation. Do not edit a published tag.
2. Confirm the default/no-default CI matrix, lint, rustdoc, and the
   `--no-default-features --features tdf` check pass. Complete relevant
   real-data smoke tests, including the ignored Bruker tests and the native
   Thermo `.raw` tests (`KOTH_MS1_RAW`, `KOTH_DIA_RAW`); ignored tests are not
   covered by the ordinary CI pass. Check the paper's
   pinned configurations if changing feature-detection behavior.
3. Commit and push the reviewed changes. In GitHub, run the Release workflow
   manually on that commit/branch. Inspect all four binary archives and the
   crate dry run; this operation publishes nothing.
4. Create and publish a GitHub Release for exactly that commit, using a tag
   `v<workspace version>`. Publishing a draft is the trigger. Use the GitHub UI
   or your own authenticated `gh release create`; releases created by another
   workflow's `GITHUB_TOKEN` do not trigger this workflow.
5. Verify the release workflow, crates.io version, and Zenodo DOI. Publication
   success in one service does not imply success in the others.

## Automated order

```text
Published GitHub Release
  ├─ verify tag/version + metadata + public repository + crate credential
  │    → full reusable CI
  │    → four binary builds + crate publish dry run
  │    → verify archive checksums → attach downloads → cargo publish
  └─ Zenodo GitHub integration → source archive and DOI
```

The CI workflow is reused on the actual release commit, so a previous green
branch build cannot substitute for release checks. Release builds include both
executables, configuration examples, documentation, license, and citation files.
The crate is built with `--locked` and published only after all target builds
and package verification succeed. There is no automatic version bump or tag
creation.

## Citation metadata policy

Keep `.zenodo.json` version-independent: no `version`, publication date,
funding, funders, or grants. The GitHub release supplies the Zenodo
version/date. `CITATION.cff` carries `version` and `date-released`, written only
by `just cite-sync` and never by hand; `just cite-check` and CI fail if the
version differs from `Cargo.toml` or the date is missing. `CITATION.cff` never
carries funding. `cff-version` identifies the file format and is not the
software version. `CITATION.cff` and the README carry only the concept DOI;
version-specific DOIs belong in papers that cite a particular release.

Keep title, description, author order, affiliations, ORCIDs, keywords, and
license synchronized. Zenodo gives `.zenodo.json` precedence over CFF when
both exist. CI validates CFF and checks this shared metadata policy.

## Failure and retry

If CI or packaging fails, no crate is uploaded and no new binary assets are
attached. Fix the problem on a new commit and use a new release tag; do not
move a tag Zenodo may already have archived. The GitHub release itself already
exists because publishing it triggered the workflow.

If uploading assets fails before crate publication, rerun failed jobs once the
cause is fixed. If `cargo publish` returns an uncertain result, check crates.io
before retrying: a version cannot be overwritten. If the crate already exists,
verify it belongs to the intended release rather than trying to republish it.
Zenodo failures are shown in its GitHub settings and are recovered there;
re-running Actions does not re-trigger Zenodo.

## References

[GitHub release events](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#release),
[Cargo publishing](https://doc.rust-lang.org/cargo/reference/publishing.html),
[Zenodo GitHub integration](https://help.zenodo.org/docs/github/), and
[Zenodo metadata precedence](https://help.zenodo.org/docs/github/describe-software/zenodo-json/).
