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

    #[error("CSV write error: {0}")]
    CsvError(#[from] csv::Error),

    #[error("JSON error: {0}")]
    JsonError(#[from] serde_json::Error),

    #[error("Config error: {0}")]
    ConfigError(String),

    #[error("No spectra found in input file")]
    NoSpectra,
}
