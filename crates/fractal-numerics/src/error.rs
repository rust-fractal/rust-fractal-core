#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MalachiteError {
    InvalidPrecision,
    Native(String),
}

impl core::fmt::Display for MalachiteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidPrecision => f.write_str("Malachite precision is invalid"),
            Self::Native(error) => f.write_str(error),
        }
    }
}

impl std::error::Error for MalachiteError {}
