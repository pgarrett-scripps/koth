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

use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    // Emitting ANY rerun-if-changed replaces cargo's default, which is "rerun
    // when any file in the package changed". So everything that can move the
    // stamp has to be listed explicitly:
    //
    //   src, Cargo.toml   editing a source file makes the tree dirty, and the
    //                     `-dirty` suffix has to appear. Listing only build.rs
    //                     here was a bug: cargo recompiled the crate on a source
    //                     edit without re-running this file, so a clean SHA
    //                     stayed baked into a binary built from an edited tree.
    //   HEAD              a checkout or branch switch changes the commit.
    //   the branch ref    `git commit` on the current branch moves the ref
    //                     HEAD points to, not HEAD itself.
    //   packed-refs       where that ref lives after `git pack-refs`/`gc`.
    //   index             `git commit` clears dirty without touching any source
    //                     file, so index is what catches dirty -> clean.
    //
    // These are resolved through the real git directory. In a linked worktree
    // `.git` is a FILE (`gitdir: ...`), and the old hard-coded `../.git/HEAD`
    // and `../.git/index` did not exist there, so the stamp never re-ran after
    // a commit: a worktree build reported the SHA of whatever commit first
    // compiled it, with a stale `-dirty`.
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=src");
    println!("cargo:rerun-if-changed=Cargo.toml");
    for p in git_watch_paths() {
        println!("cargo:rerun-if-changed={}", p.display());
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

/// Files whose change must re-run this script: the worktree's HEAD and index,
/// the ref HEAD points to, and packed-refs. Empty when no git directory is found
/// (a source tarball), which leaves the `src`/`Cargo.toml` triggers in charge.
fn git_watch_paths() -> Vec<PathBuf> {
    let Some(git_dir) = find_git_dir() else {
        return Vec::new();
    };
    // A linked worktree keeps HEAD and index in its own gitdir but shares refs
    // with the main repository, named by the `commondir` file.
    let common_dir = std::fs::read_to_string(git_dir.join("commondir"))
        .ok()
        .map(|c| git_dir.join(c.trim()))
        .unwrap_or_else(|| git_dir.clone());

    let mut paths = vec![git_dir.join("HEAD"), git_dir.join("index")];
    if let Ok(head) = std::fs::read_to_string(git_dir.join("HEAD")) {
        if let Some(r) = head.trim().strip_prefix("ref:") {
            let r = r.trim();
            // Per-worktree refs (bisect, worktree/) live in the gitdir; branch
            // refs live in the common dir. Watch whichever exists, and always
            // packed-refs, where a loose ref goes after packing.
            let local = git_dir.join(r);
            paths.push(if local.exists() {
                local
            } else {
                common_dir.join(r)
            });
        }
    }
    paths.push(common_dir.join("packed-refs"));
    // A path that does not exist yet (no packed-refs) still has to be listed:
    // cargo re-runs when it appears. Canonicalise what exists so the output
    // does not carry `../..` segments.
    paths
        .into_iter()
        .map(|p| p.canonicalize().unwrap_or(p))
        .collect()
}

/// The git directory for this checkout: `.git` itself when it is a directory,
/// or the target of its `gitdir:` line when it is a file (linked worktree,
/// submodule). Searches upward from the crate so the workspace layout is not
/// hard-coded.
fn find_git_dir() -> Option<PathBuf> {
    let start = PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR")?);
    let mut dir: Option<&Path> = Some(&start);
    while let Some(d) = dir {
        let dot_git = d.join(".git");
        if dot_git.is_dir() {
            return Some(dot_git);
        }
        if dot_git.is_file() {
            let text = std::fs::read_to_string(&dot_git).ok()?;
            let target = text.lines().find_map(|l| l.strip_prefix("gitdir:"))?.trim();
            return Some(d.join(target));
        }
        dir = d.parent();
    }
    None
}
