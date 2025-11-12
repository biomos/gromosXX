//! I/O module for reading and writing GROMOS file formats
//!
//! Supports:
//! - .topo/.top - Topology files
//! - .conf/.cnf - Coordinate files
//! - .imd - Input parameter files (simulation settings)
//! - .trc/.trj - Trajectory files (coordinates over time)
//! - .tre - Energy files (ENE/ENA format)
//! - .trf - Force files
//! - .ptp - Perturbation topology files (FEP)

pub mod coordinate;
pub mod topology;
pub mod imd;
pub mod trajectory;
pub mod energy;
pub mod force;
pub mod pdb;
pub mod g96;
pub mod ptp;
pub mod dlg;

// Re-export commonly used types
pub use imd::{ImdParameters, TempBathParameters, PressureParameters};
pub use trajectory::TrajectoryWriter;
pub use energy::{EnergyWriter, EnergyFrame, EnergyBlock};
pub use force::ForceWriter;
pub use ptp::PtpWriter;
pub use dlg::{DlgWriter, LambdaDerivativeFrame};

use std::io;

/// Common error type for I/O operations
#[derive(Debug)]
pub enum IoError {
    FileNotFound(String),
    ParseError(String),
    FormatError(String),
    WriteError(String),
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
            IoError::WriteError(msg) => write!(f, "Write error: {}", msg),
            IoError::Io(err) => write!(f, "I/O error: {}", err),
        }
    }
}

impl std::error::Error for IoError {}
