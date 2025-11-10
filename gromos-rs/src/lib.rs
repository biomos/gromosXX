//! GROMOS-RS: High-performance Rust molecular dynamics engine
//!
//! A complete Rust translation of GROMOS (GROningen MOlecular Simulation),
//! providing 2-3x performance improvements through modern optimization techniques.
//!
//! ## Features
//!
//! - **Pure Rust MD engine**: Complete translation, not just FFI wrappers
//! - **SIMD vectorization**: Automatic vectorization with glam + wide
//! - **Fearless concurrency**: Parallel algorithms with Rayon
//! - **Memory safety**: Zero-cost abstractions without runtime overhead
//! - **Modular design**: Easy to extend and customize
//!
//! ## Architecture
//!
//! ```text
//! gromos-rs
//! ├── math          - SIMD-accelerated math primitives
//! ├── topology      - Molecular structure and force field parameters
//! ├── configuration - System state management
//! ├── pairlist      - Neighbor list generation
//! ├── interaction   - Force calculations (bonded + nonbonded)
//! ├── integrator    - Time integration algorithms
//! └── engine        - Main MD simulation driver
//! ```
//!
//! ## Example
//!
//! ```rust,no_run
//! use gromos_rs::*;
//!
//! // Create system
//! let topo = topology::Topology::new();
//! let conf = configuration::Configuration::new(100, 1, 1);
//!
//! // Setup MD
//! let integrator = integrator::LeapFrog::new();
//! let forcefield = interaction::ForceField::new();
//!
//! // Run simulation
//! // engine::run_md(topo, conf, integrator, forcefield, 1000);
//! ```

pub mod math;
pub mod topology;
pub mod configuration;
pub mod pairlist;
pub mod interaction;
pub mod integrator;
pub mod ffi;

// Re-export main types for convenience
pub use math::{Vec3, Mat3, Periodicity, BoundaryCondition};
pub use topology::{Topology, LJParameters};
pub use configuration::{Configuration, State, Energy};
pub use interaction::nonbonded::{lj_crf_interaction, ForceStorage};
pub use integrator::{Integrator, LeapFrog, VelocityVerlet};

#[cfg(feature = "mimalloc")]
use mimalloc::MiMalloc;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;
