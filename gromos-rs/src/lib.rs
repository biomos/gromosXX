//! GROMOS-RS: High-performance Rust kernels for molecular dynamics
//!
//! This library provides optimized implementations of performance-critical
//! components of the GROMOS molecular dynamics package, leveraging Rust's
//! zero-cost abstractions, SIMD capabilities, and fearless concurrency.

pub mod math;
pub mod interaction;
pub mod ffi;

// Re-export main types
pub use math::{Vec3, Periodicity, BoundaryCondition};
pub use interaction::nonbonded::{LJParameters, lj_crf_interaction};

#[cfg(feature = "mimalloc")]
use mimalloc::MiMalloc;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;
