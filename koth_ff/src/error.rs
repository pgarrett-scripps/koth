use thiserror::Error;

#[derive(Error, Debug)]
pub enum KothError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Unsupported input format: {0}")]
    UnsupportedFormat(String),

    #[error("mzML parse error: {0}")]
    MzmlError(String),

    #[cfg(feature = "tdf")]
    #[error("Bruker TDF error: {0}")]
    TdfError(String),

    // Ungated: the `thermo` feature does not enable `tdf`, so the Thermo reader
    // needs an error variant that exists in a thermo-only build.
    #[error("Thermo .raw error: {0}")]
    ThermoError(String),

    #[error("CSV write error: {0}")]
    CsvError(#[from] csv::Error),

    #[error("JSON error: {0}")]
    JsonError(#[from] serde_json::Error),

    #[error("Config error: {0}")]
    ConfigError(String),

    #[error("No spectra found in input file")]
    NoSpectra,

    #[error("Parquet error: {0}")]
    ParquetError(String),
}

/// `koth-core` errors convert one-to-one: its `Io` becomes [`KothError::Io`] and
/// its `Config` becomes [`KothError::ConfigError`], with the same message.
impl From<koth_core::Error> for KothError {
    fn from(e: koth_core::Error) -> Self {
        match e {
            koth_core::Error::Io(e) => KothError::Io(e),
            koth_core::Error::Config(msg) => KothError::ConfigError(msg),
            other => KothError::ConfigError(other.to_string()),
        }
    }
}
