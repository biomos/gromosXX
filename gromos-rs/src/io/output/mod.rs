/// Output file writers for GROMOS formats
///
/// This module contains writers for various output files generated
/// during GROMOS MD simulations, particularly for advanced sampling methods.

pub mod gamd_stats;

pub use gamd_stats::{GamdStatsWriter, GamdBoostWriter};
