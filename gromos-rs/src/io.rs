//! I/O module for reading and writing GROMOS file formats
//!
//! Supports:
//! - .topo/.top - Topology files
//! - .conf/.cnf - Coordinate files
//! - .trc - Trajectory files
//! - .tre - Energy files

pub mod coordinate;
pub mod topology;

use std::io;

/// Common error type for I/O operations
#[derive(Debug)]
pub enum IoError {
    FileNotFound(String),
    ParseError(String),
    FormatError(String),
    Io(io::Error),
}

impl From<io::Error> for IoError {
    fn from(err: io::Error) -> Self {
        IoError::Io(err)
    }
}

impl std::fmt::Display for IoError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            IoError::FileNotFound(path) => write!(f, "File not found: {}", path),
            IoError::ParseError(msg) => write!(f, "Parse error: {}", msg),
            IoError::FormatError(msg) => write!(f, "Format error: {}", msg),
            IoError::Io(err) => write!(f, "I/O error: {}", err),
        }
    }
}

impl std::error::Error for IoError {}
