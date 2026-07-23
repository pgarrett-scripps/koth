# koth_rust

## Testing conventions

Large unit-test modules live in a **sibling file** next to the module they test,
wired in via `#[path]` rather than inlined at the bottom of the source file:

```rust
// in foo.rs
#[cfg(test)]
#[path = "foo_tests.rs"]
mod tests;
```

- A file named `<stem>_tests.rs` (e.g. `assemble_tests.rs`, `warp_tests.rs`) holds
  the unit tests for its same-named production module (`assemble.rs`, `warp.rs`).
  Because `#[path]` makes it a **child module** of the production module, its tests
  see all private items exactly as an inline `mod tests` would (they start with
  `use super::*;`). These are unit tests, **not** integration tests — do not move
  them to `tests/` and do not "helpfully" re-inline them back into the source file.
- Only large test blocks (roughly >= 80 lines) are split out this way; small test
  modules are kept inline. When adding tests to a module that already has a
  `<stem>_tests.rs`, put them in that sibling file.
