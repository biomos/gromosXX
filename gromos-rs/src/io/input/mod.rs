/// Input file parsers for GROMOS .imd format
///
/// This module contains parsers for various input file blocks
/// used in GROMOS MD simulations.

pub mod gamd_block;

pub use gamd_block::GamdBlock;
