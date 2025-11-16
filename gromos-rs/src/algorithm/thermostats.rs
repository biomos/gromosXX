//! Temperature coupling algorithms (thermostats)
//!
//! This module implements various thermostat algorithms:
//! - Berendsen: Weak coupling thermostat
//! - Nosé-Hoover: Deterministic extended system thermostat
//! - Andersen: Stochastic collision thermostat

use crate::topology::Topology;
use crate::configuration::Configuration;
use crate::math::Vec3;

/// Berendsen thermostat parameters
#[derive(Debug, Clone)]
pub struct BerendsenThermostatParameters {
    pub target_temperature: f64,  // Target temperature (K)
    pub coupling_time: f64,        // Coupling time constant τ (ps)
}

impl Default for BerendsenThermostatParameters {
    fn default() -> Self {
        Self {
            target_temperature: 300.0,  // 300 K
            coupling_time: 0.1,         // 0.1 ps
        }
    }
}

/// Nosé-Hoover thermostat parameters
#[derive(Debug, Clone)]
pub struct NoseHooverThermostatParameters {
    pub target_temperature: f64,  // Target temperature (K)
    pub coupling_time: f64,        // Coupling time constant τ (ps)
    pub xi: f64,                   // Thermostat variable
}

impl Default for NoseHooverThermostatParameters {
    fn default() -> Self {
        Self {
            target_temperature: 300.0,  // 300 K
            coupling_time: 0.1,         // 0.1 ps
            xi: 0.0,                    // Initial thermostat variable
        }
    }
}

/// Andersen thermostat parameters
#[derive(Debug, Clone)]
pub struct AndersenThermostatParameters {
    pub target_temperature: f64,  // Target temperature (K)
    pub collision_frequency: f64,  // Collision frequency ν (ps⁻¹)
}

impl Default for AndersenThermostatParameters {
    fn default() -> Self {
        Self {
            target_temperature: 300.0,  // 300 K
            collision_frequency: 1.0,   // 1 ps⁻¹
        }
    }
}

/// Boltzmann constant in kJ/(mol·K)
const K_B: f64 = 0.008314462618;

/// Apply Berendsen thermostat (weak coupling)
///
/// The Berendsen thermostat (Berendsen et al., 1984) rescales velocities
/// to weakly couple the system to a heat bath at target temperature.
///
/// Scaling factor: λ = sqrt(1 + (dt/τ) * (T₀/T - 1))
///
/// where:
/// - dt is the time step
/// - τ is the coupling time constant
/// - T₀ is the target temperature
/// - T is the current temperature
///
/// # Parameters
/// - `topo`: Molecular topology
/// - `conf`: Configuration with velocities to scale
/// - `dt`: Time step size (ps)
/// - `params`: Thermostat parameters
pub fn berendsen_thermostat(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &BerendsenThermostatParameters,
) {
    // Calculate current kinetic energy and temperature
    let mut kinetic_energy = 0.0;
    let num_atoms = topo.num_atoms();

    for i in 0..num_atoms {
        let v = conf.current().vel[i];
        let m = topo.mass[i];
        kinetic_energy += 0.5 * m * (v.dot(v) as f64);
    }

    // Temperature from kinetic energy
    // T = 2 * E_kin / (N_dof * k_B)
    // where N_dof = 3N - N_constraints (approximately 3N for large systems)
    let n_dof = (3 * num_atoms) as f64;
    let current_temperature = 2.0 * kinetic_energy / (n_dof * K_B);

    if current_temperature < 1e-6 {
        return; // Avoid division by zero
    }

    // Berendsen scaling factor
    // λ = sqrt(1 + (dt/τ) * (T₀/T - 1))
    let temp_ratio = params.target_temperature / current_temperature;
    let coupling_factor = dt / params.coupling_time;
    let lambda = (1.0 + coupling_factor * (temp_ratio - 1.0)).sqrt();

    // Scale all velocities
    for i in 0..num_atoms {
        conf.current_mut().vel[i] *= lambda as f32;
    }
}

/// Apply Nosé-Hoover thermostat (extended system)
///
/// The Nosé-Hoover thermostat (Nosé, 1984; Hoover, 1985) introduces
/// an additional degree of freedom (thermostat variable ξ) that
/// couples to the system to maintain constant temperature.
///
/// Equations of motion:
/// - dv/dt = F/m - ξ*v
/// - dξ/dt = (T/T₀ - 1) / τ²
///
/// # Parameters
/// - `topo`: Molecular topology
/// - `conf`: Configuration with velocities
/// - `dt`: Time step size (ps)
/// - `params`: Thermostat parameters (mutable to update ξ)
pub fn nose_hoover_thermostat(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &mut NoseHooverThermostatParameters,
) {
    // Calculate current temperature
    let mut kinetic_energy = 0.0;
    let num_atoms = topo.num_atoms();

    for i in 0..num_atoms {
        let v = conf.current().vel[i];
        let m = topo.mass[i];
        kinetic_energy += 0.5 * m * (v.dot(v) as f64);
    }

    let n_dof = (3 * num_atoms) as f64;
    let current_temperature = 2.0 * kinetic_energy / (n_dof * K_B);

    if current_temperature < 1e-6 {
        return;
    }

    // Update thermostat variable ξ
    // dξ/dt = (T/T₀ - 1) / τ²
    let temp_ratio = current_temperature / params.target_temperature;
    let tau_sq = params.coupling_time * params.coupling_time;
    params.xi += dt * (temp_ratio - 1.0) / tau_sq;

    // Apply friction to velocities
    // v_new = v * exp(-ξ * dt)
    let friction_factor = (-params.xi * dt).exp();

    for i in 0..num_atoms {
        conf.current_mut().vel[i] *= friction_factor as f32;
    }
}

/// Apply Andersen thermostat (stochastic collisions)
///
/// The Andersen thermostat (Andersen, 1980) stochastically assigns
/// new velocities from Maxwell-Boltzmann distribution to random atoms,
/// simulating collisions with a heat bath.
///
/// For each atom, with probability ν*dt:
/// - Draw new velocity from Maxwell-Boltzmann distribution at T₀
///
/// # Parameters
/// - `topo`: Molecular topology
/// - `conf`: Configuration with velocities
/// - `dt`: Time step size (ps)
/// - `params`: Thermostat parameters
///
/// # Note
/// Requires random number generation - simplified implementation
pub fn andersen_thermostat(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &AndersenThermostatParameters,
) {
    // Probability of collision per atom
    let collision_prob = params.collision_frequency * dt;

    let num_atoms = topo.num_atoms();

    for i in 0..num_atoms {
        // Simple deterministic version: periodically reassign velocities
        // A full implementation would use proper random number generation

        // For now, use a simple pseudo-random check based on time and atom index
        let pseudo_random = ((i as f64 * 0.123456) % 1.0);

        if pseudo_random < collision_prob {
            // Draw new velocity from Maxwell-Boltzmann distribution
            // v ~ N(0, sqrt(k_B * T / m))
            let m = topo.mass[i];
            let sigma = (K_B * params.target_temperature / m).sqrt();

            // Simple Box-Muller transform for Gaussian random numbers
            // (In production, use a proper RNG library)
            let u1 = 0.5 + 0.1 * ((i as f64 * 0.456789) % 1.0);
            let u2 = 0.5 + 0.1 * ((i as f64 * 0.789012) % 1.0);

            let z0 = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            let z1 = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).sin();
            let z2 = (-2.0 * u2.ln()).sqrt() * (2.0 * std::f64::consts::PI * u1).cos();

            conf.current_mut().vel[i] = Vec3::new(
                (sigma * z0) as f32,
                (sigma * z1) as f32,
                (sigma * z2) as f32,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_berendsen_parameters() {
        let params = BerendsenThermostatParameters::default();
        assert!(params.target_temperature > 0.0);
        assert!(params.coupling_time > 0.0);
    }

    #[test]
    fn test_nose_hoover_parameters() {
        let params = NoseHooverThermostatParameters::default();
        assert!(params.target_temperature > 0.0);
        assert!(params.coupling_time > 0.0);
    }

    #[test]
    fn test_andersen_parameters() {
        let params = AndersenThermostatParameters::default();
        assert!(params.target_temperature > 0.0);
        assert!(params.collision_frequency > 0.0);
    }
}
