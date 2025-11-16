//! Molecular dynamics algorithms
//!
//! This module contains various MD algorithms:
//! - Constraints (SHAKE, SETTLE, M-SHAKE)
//! - Thermostats (Berendsen, Nosé-Hoover, Andersen)
//! - Barostats (Berendsen, Parrinello-Rahman)

pub mod constraints;
pub mod thermostats;
pub mod barostats;

// Re-export commonly used items
pub use constraints::{
    ShakeParameters,
    ConstraintResult,
    shake,
    m_shake,
    settle,
};

pub use thermostats::{
    BerendsenThermostatParameters,
    NoseHooverThermostatParameters,
    AndersenThermostatParameters,
    berendsen_thermostat,
    nose_hoover_thermostat,
    andersen_thermostat,
};

pub use barostats::{
    BerendsenBarostatParameters,
    ParrinelloRahmanBarostatParameters,
    berendsen_barostat,
    parrinello_rahman_barostat,
    calculate_virial,
};
