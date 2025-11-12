/// Replica Exchange Molecular Dynamics - Individual Replica
///
/// This module defines the data structures and methods for a single replica
/// in a replica exchange simulation. Each replica maintains its own configuration
/// and simulation parameters (temperature, lambda, etc.).

use crate::configuration::Configuration;
use crate::topology::Topology;
use crate::integrator::Integrator;
use crate::algorithm::thermostats::{BerendsenThermostatParameters, berendsen_thermostat};
use crate::algorithm::constraints::{shake, ShakeParameters};
use crate::interaction::nonbonded::{lj_crf_innerloop, lj_crf_innerloop_parallel, CRFParameters, ForceStorage};
use crate::interaction::bonded::calculate_bonded_forces;
use crate::pairlist::{PairlistContainer, StandardPairlistAlgorithm};
use crate::fep::LambdaController;
use crate::math::Rectangular;
use std::sync::Arc;

/// Unique identifier for a replica
pub type ReplicaId = usize;

/// Replica information and state
#[derive(Debug, Clone)]
pub struct ReplicaInfo {
    /// Unique replica ID
    pub id: ReplicaId,

    /// Target temperature (K)
    pub temperature: f64,

    /// Lambda value for FEP (0.0 to 1.0)
    pub lambda: f64,

    /// Timestep (ps)
    pub dt: f64,

    /// Current potential energy
    pub potential_energy: f64,

    /// Partner replica ID for exchange
    pub partner_id: Option<ReplicaId>,

    /// Exchange probability with partner
    pub exchange_probability: f64,

    /// Whether exchange was accepted
    pub exchange_accepted: bool,

    /// Current run/step number
    pub run: usize,

    /// Box volume (for NPT)
    pub volume: f64,
}

impl ReplicaInfo {
    /// Create new replica info
    pub fn new(id: ReplicaId, temperature: f64, lambda: f64, dt: f64) -> Self {
        ReplicaInfo {
            id,
            temperature,
            lambda,
            dt,
            potential_energy: 0.0,
            partner_id: None,
            exchange_probability: 0.0,
            exchange_accepted: false,
            run: 0,
            volume: 0.0,
        }
    }

    /// Calculate inverse temperature (beta = 1/kT)
    /// Uses kB = 0.008314462618 kJ/(mol·K) (GROMOS unit)
    pub fn beta(&self) -> f64 {
        const KB: f64 = 0.008314462618; // kJ/(mol·K)
        1.0 / (KB * self.temperature)
    }
}

/// Individual replica containing configuration and metadata
pub struct Replica {
    /// Replica metadata
    pub info: ReplicaInfo,

    /// System configuration (positions, velocities, forces, energies)
    pub configuration: Configuration,

    /// Lambda controller for FEP
    pub lambda_controller: Option<LambdaController>,

    /// Pairlist for nonbonded interactions
    pairlist: PairlistContainer,

    /// Pairlist algorithm
    pairlist_algorithm: StandardPairlistAlgorithm,

    /// Use parallel force calculation
    parallel_forces: bool,

    /// CRF parameters for nonbonded interactions
    crf_params: CRFParameters,
}

impl Replica {
    /// Create a new replica
    pub fn new(
        info: ReplicaInfo,
        configuration: Configuration,
        lambda_controller: Option<LambdaController>,
        cutoff: f64,
    ) -> Self {
        let pairlist = PairlistContainer::new(cutoff, cutoff, 0.2);
        let pairlist_algorithm = StandardPairlistAlgorithm::new(true);

        // Calculate CRF parameters
        // crf = (eps_rf - eps) * (2*eps + eps_rf) / (eps_rf + 2*eps) / cutoff^3
        let eps = 1.0;
        let eps_rf = 61.0;
        let crf = (eps_rf - eps) * (2.0 * eps + eps_rf) / (eps_rf + 2.0 * eps) / (cutoff * cutoff * cutoff);

        let crf_params = CRFParameters {
            crf_cut: cutoff,
            crf_2cut3i: crf / (2.0 * cutoff * cutoff * cutoff),
            crf_cut3i: (1.0 - crf / 2.0) / cutoff,
        };

        Replica {
            info,
            configuration,
            lambda_controller,
            pairlist,
            pairlist_algorithm,
            parallel_forces: true,
            crf_params,
        }
    }

    /// Update pairlist if needed
    pub fn update_pairlist_if_needed(&mut self, topology: &Topology) {
        if self.pairlist.needs_update() {
            let box_dims = self.configuration.current().box_config.dimensions();
            let periodicity = Rectangular::new(box_dims);

            self.pairlist_algorithm.update(
                topology,
                &self.configuration,
                &mut self.pairlist,
                &periodicity,
            );
            self.pairlist.reset_counter();
        } else {
            self.pairlist.step();
        }
    }

    /// Calculate all forces for this replica
    pub fn calculate_forces(&mut self, topology: &Topology) {
        // Zero forces and energies
        self.configuration.current_mut().clear_forces();
        self.configuration.current_mut().energies.clear();

        // Bonded forces (bonds + angles + dihedrals)
        let bonded_result = calculate_bonded_forces(topology, &self.configuration, true);

        // Accumulate bonded forces and energies
        for (i, &force) in bonded_result.forces.iter().enumerate() {
            self.configuration.current_mut().force[i] += force;
        }
        self.configuration.current_mut().energies.bond_total = bonded_result.energy;

        // Create periodicity from box dimensions
        let box_dims = self.configuration.current().box_config.dimensions();
        let periodicity = Rectangular::new(box_dims);

        // Convert LJ parameters to nonbonded format
        use crate::interaction::nonbonded::LJParameters as NBLJParams;
        let lj_params: Vec<Vec<NBLJParams>> = topology.lj_parameters.iter()
            .map(|row| {
                row.iter()
                    .map(|params| NBLJParams {
                        c6: params.c6,
                        c12: params.c12,
                    })
                    .collect()
            })
            .collect();

        // Nonbonded forces (LJ + CRF)
        let mut nonbonded_storage = if self.parallel_forces {
            // Convert pairlist to (u32, u32) format
            let pairlist_short: Vec<(u32, u32)> = self.pairlist.solute_short.iter()
                .map(|&(i, j)| (i as u32, j as u32))
                .collect();

            // Convert charge from Vec<f64> to Vec<f32> for compatibility
            let charges_f32: Vec<f32> = topology.charge.iter().map(|&q| q as f32).collect();

            // Convert iac from Vec<usize> to Vec<u32> for compatibility
            let iac_u32: Vec<u32> = topology.iac.iter().map(|&i| i as u32).collect();

            lj_crf_innerloop_parallel(
                &self.configuration.current().pos,
                &charges_f32,
                &iac_u32,
                &pairlist_short,
                &lj_params,
                &self.crf_params,
                &periodicity,
                topology.num_atoms(),
            )
        } else {
            // Convert pairlist to (u32, u32) format
            let pairlist_short: Vec<(u32, u32)> = self.pairlist.solute_short.iter()
                .map(|&(i, j)| (i as u32, j as u32))
                .collect();

            // Convert charge from Vec<f64> to Vec<f32> for compatibility
            let charges_f32: Vec<f32> = topology.charge.iter().map(|&q| q as f32).collect();

            // Convert iac from Vec<usize> to Vec<u32> for compatibility
            let iac_u32: Vec<u32> = topology.iac.iter().map(|&i| i as u32).collect();

            let mut storage = ForceStorage::new(topology.num_atoms());
            lj_crf_innerloop(
                &self.configuration.current().pos,
                &charges_f32,
                &iac_u32,
                &pairlist_short,
                &lj_params,
                &self.crf_params,
                &periodicity,
                &mut storage,
            );
            storage
        };

        // Accumulate nonbonded forces and energies
        for (i, &force) in nonbonded_storage.forces.iter().enumerate() {
            self.configuration.current_mut().force[i] += force;
        }
        self.configuration.current_mut().energies.lj_total = nonbonded_storage.e_lj;
        self.configuration.current_mut().energies.crf_total = nonbonded_storage.e_crf;

        // Update potential energy total
        self.configuration.current_mut().energies.update_potential_total();
        self.info.potential_energy = self.configuration.current().energies.potential_total;
    }

    /// Perform one MD step
    pub fn md_step<I: Integrator>(
        &mut self,
        topology: &Topology,
        integrator: &mut I,
        thermostat_params: Option<&BerendsenThermostatParameters>,
        shake_params: Option<&ShakeParameters>,
    ) {
        // Update pairlist if needed
        self.update_pairlist_if_needed(topology);

        // Calculate forces
        self.calculate_forces(topology);

        // Integrate
        integrator.step(self.info.dt, topology, &mut self.configuration);

        // Apply constraints (SHAKE)
        if let Some(shake_p) = shake_params {
            shake(topology, &mut self.configuration, self.info.dt, shake_p);
        }

        // Calculate kinetic energy and temperature
        self.configuration.current_mut().calculate_kinetic_energy(&topology.mass);

        // Apply thermostat
        if let Some(params) = thermostat_params {
            berendsen_thermostat(topology, &mut self.configuration, self.info.dt, params);
        }

        // Update lambda if FEP is active
        if let Some(ref mut lambda_ctrl) = self.lambda_controller {
            lambda_ctrl.update(self.info.dt);
            self.info.lambda = lambda_ctrl.lambda;
        }

        // Update volume
        self.info.volume = self.configuration.current().box_config.volume();

        // Increment run counter
        self.info.run += 1;
    }

    /// Run multiple MD steps
    pub fn run_md_steps<I: Integrator>(
        &mut self,
        topology: &Topology,
        integrator: &mut I,
        n_steps: usize,
        thermostat_params: Option<&BerendsenThermostatParameters>,
        shake_params: Option<&ShakeParameters>,
    ) {
        for _ in 0..n_steps {
            self.md_step(topology, integrator, thermostat_params, shake_params);
        }
    }

    /// Get current potential energy
    pub fn potential_energy(&self) -> f64 {
        self.configuration.current().energies.potential_total
    }

    /// Get current kinetic energy
    pub fn kinetic_energy(&self) -> f64 {
        self.configuration.current().energies.kinetic_total
    }

    /// Get current total energy
    pub fn total_energy(&self) -> f64 {
        self.configuration.current().energies.total()
    }

    /// Get current temperature (assuming 3N degrees of freedom)
    pub fn temperature(&self, topology: &Topology) -> f64 {
        let n_dof = 3 * topology.num_atoms();
        self.configuration.current().temperature(n_dof)
    }

    /// Exchange configuration with another replica
    pub fn exchange_configuration(&mut self, other: &mut Replica) {
        std::mem::swap(&mut self.configuration, &mut other.configuration);
    }

    /// Rescale velocities to match target temperature
    pub fn rescale_velocities(&mut self, topology: &Topology, target_temp: f64) {
        let n_dof = 3 * topology.num_atoms();
        let current_temp = self.configuration.current().temperature(n_dof);
        if current_temp > 0.0 {
            let scale_factor = (target_temp / current_temp).sqrt();

            // Scale velocities
            for vel in self.configuration.current_mut().vel.iter_mut() {
                *vel = *vel * scale_factor as f32;
            }

            // Recalculate kinetic energy
            self.configuration.current_mut().calculate_kinetic_energy(&topology.mass);
        }
    }

    /// Set parallel force calculation
    pub fn set_parallel_forces(&mut self, parallel: bool) {
        self.parallel_forces = parallel;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_replica_info_beta() {
        let info = ReplicaInfo::new(0, 300.0, 0.0, 0.002);
        let beta = info.beta();

        // beta = 1/(kB*T) = 1/(0.008314462618 * 300) ≈ 0.4006
        assert!((beta - 0.4006).abs() < 0.0001);
    }

    #[test]
    fn test_replica_info_creation() {
        let info = ReplicaInfo::new(5, 350.0, 0.5, 0.001);

        assert_eq!(info.id, 5);
        assert_eq!(info.temperature, 350.0);
        assert_eq!(info.lambda, 0.5);
        assert_eq!(info.dt, 0.001);
        assert_eq!(info.potential_energy, 0.0);
        assert_eq!(info.partner_id, None);
        assert!(!info.exchange_accepted);
    }
}
