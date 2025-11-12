/// Input file parsers for GROMOS .imd format
///
/// This module contains parsers for various input file blocks
/// used in GROMOS MD simulations.

pub mod gamd_block;
pub mod eds_block;
pub mod replica_block;

pub use gamd_block::GamdBlock;
pub use eds_block::EdsBlock;
pub use replica_block::ReplicaBlock;
