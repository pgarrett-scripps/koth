//! The error type of the in-memory stages and the detection configuration.

use thiserror::Error;

/// Convenience alias for results produced by this crate.
pub type Result<T, E = Error> = std::result::Result<T, E>;

/// Anything that can go wrong in `koth-core`: reading or validating a
/// configuration. The messages match the `koth-ms` crate's `KothError` variants
/// of the same meaning (`Io`, `ConfigError`), which this converts into one to one.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum Error {
    /// Reading a configuration file failed.
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// An invalid or unparseable configuration value.
    #[error("Config error: {0}")]
    Config(String),
}
