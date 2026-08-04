//! Embed the git commit into the binary at compile time.
//!
//! Without this, both binaries reported a bare `koth_ff 0.1.0` -- a
//! hand-maintained Cargo version that had not moved across 173 changelog lines
//! of behaviour changes. Every build since 0.1.0 identified itself identically,
//! so no output file could be traced back to the code that produced it, and the
//! manuscript's "producing commit" could only ever be filled in from memory.
//!
//! The benchmark in the koth-paper repository pins a specific commit and refuses
//! to run against anything else (`just koth-check`). That check reads the string
//! this file produces, so this is the thing that makes the pin enforceable
//! rather than decorative.
//!
//! Degrades quietly. A build from a source tarball with no `.git` (a Zenodo
//! archive, a `cargo package`) reports `unknown` rather than failing, because a
//! missing git directory is not a build error.

use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    // Rebuild when HEAD moves, so a checkout or commit refreshes the stamp
    // instead of leaving a stale SHA baked into an otherwise-fresh binary.
    for p in ["../.git/HEAD", "../.git/index"] {
        if std::path::Path::new(p).exists() {
            println!("cargo:rerun-if-changed={p}");
        }
    }
    // Escape hatch for reproducible/offline builds that want to state the commit
    // explicitly (packaging, CI from a tarball, `SOURCE_DATE_EPOCH`-style flows).
    println!("cargo:rerun-if-env-changed=KOTH_BUILD_GIT_SHA");

    let sha = std::env::var("KOTH_BUILD_GIT_SHA")
        .ok()
        .filter(|s| !s.is_empty())
        .or_else(git_describe)
        .unwrap_or_else(|| "unknown".to_string());

    println!("cargo:rustc-env=KOTH_GIT_SHA={sha}");
}

/// Short SHA, suffixed `-dirty` when the working tree has uncommitted changes.
///
/// The dirty flag is the point: a number produced by an edited tree is not
/// reproducible from any commit, and that has to be visible in the output rather
/// than discovered later.
fn git_describe() -> Option<String> {
    let sha = run(&["rev-parse", "--short=12", "HEAD"])?;
    // `--quiet` makes this exit non-zero when there ARE differences, which is
    // how we detect a dirty tree without parsing anything.
    let clean = Command::new("git")
        .args(["diff", "--quiet", "HEAD"])
        .status()
        .map(|s| s.success())
        .unwrap_or(true);
    Some(if clean { sha } else { format!("{sha}-dirty") })
}

fn run(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8(out.stdout).ok()?.trim().to_string();
    (!s.is_empty()).then_some(s)
}
