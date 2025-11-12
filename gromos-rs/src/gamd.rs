/// Gaussian Accelerated Molecular Dynamics (GaMD)
///
/// GaMD enhances sampling by adding a boost potential that smoothens
/// the energy landscape, accelerating transitions over barriers without
/// requiring prior knowledge of transition pathways.
///
/// ## Theory
///
/// The boost potential is applied when the system potential energy
/// is below a threshold E:
///
/// ```text
/// ΔV(r) = ½K(E - V(r))²  when V(r) < E
/// ΔV(r) = 0              when V(r) ≥ E
/// ```
///
/// Where:
/// - K = force constant (computed from statistics)
/// - E = energy threshold (Vmax for lower bound, derived for upper bound)
/// - V(r) = system potential energy
///
/// The boosted force becomes:
/// ```text
/// F_boost = F_orig × [1 + K(E - V)]
/// ```
///
/// ## Boost Modes
///
/// - **Dihedral boost**: Accelerate only dihedral energies
/// - **Total boost**: Accelerate entire potential (including dihedrals)
/// - **Dual boost**: Separate acceleration for dihedrals and non-dihedrals
///
/// ## Parameter Computation
///
/// ### Force Constant K:
/// ```text
/// K = k0 / (Vmax - Vmin)
/// ```
///
/// ### k0 (threshold factor, 0 < k0 ≤ 1):
/// **Lower bound:** k0 = min(1.0, (σ₀/σᵥ) × (Vmax - Vmin) / (Vmax - Vmean))
/// **Upper bound:** k0 = (1 - σ₀/σᵥ) × (Vmax - Vmin) / (Vmean - Vmin)
///
/// ### Energy Threshold E:
/// **Lower bound:** E = Vmax
/// **Upper bound:** E = Vmin + (Vmax - Vmin) / k0
///
/// ## Workflow
///
/// 1. **CMD Search Phase**: Classical MD, collect statistics (no boost)
/// 2. **GaMD Search Phase**: Accelerated MD with adaptive parameters
/// 3. **Production Phase**: Fixed parameters for production runs
///
/// ## References
/// - Miao et al., J. Chem. Theory Comput. 2015, 11, 3584-3595
/// - Pang et al., J. Chem. Theory Comput. 2017, 13, 9-19

use crate::configuration::Configuration;
use crate::topology::Topology;
use crate::math::Vec3;
use crate::integrator::Integrator;
use crate::algorithm::thermostats::{BerendsenThermostatParameters, berendsen_thermostat};
use crate::algorithm::constraints::{shake, ShakeParameters};
use crate::interaction::bonded::calculate_bonded_forces;
use crate::interaction::nonbonded::{lj_crf_innerloop, lj_crf_innerloop_parallel, CRFParameters, ForceStorage, LJParameters as NBLJParams};
use crate::pairlist::{PairlistContainer, StandardPairlistAlgorithm};
use crate::math::Rectangular;

/// GaMD search mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchMode {
    /// Production: use fixed parameters, no statistics update
    NoSearch,
    /// CMD search: classical MD, collect statistics without boost
    CmdSearch,
    /// GaMD search: accelerated MD with adaptive parameters
    GamdSearch,
}

/// GaMD boost form
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoostForm {
    /// Boost only dihedral potential
    DihedralBoost,
    /// Boost total potential (including dihedrals)
    TotalBoost,
    /// Separate boost for dihedral and non-dihedral
    DualBoost,
}

/// Energy threshold type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThresholdType {
    /// Lower bound: E = Vmax
    LowerBound,
    /// Upper bound: E = Vmin + (Vmax - Vmin) / k0
    UpperBound,
}

/// GaMD statistics for a single energy component
#[derive(Debug, Clone)]
pub struct GamdStatistics {
    /// Maximum energy observed
    pub v_max: f64,
    /// Minimum energy observed
    pub v_min: f64,
    /// Mean energy
    pub v_mean: f64,
    /// Standard deviation of energy
    pub sigma_v: f64,
    /// Welford M2 accumulator for variance calculation
    pub m2: f64,
    /// Number of steps collected
    pub steps: usize,
}

impl GamdStatistics {
    /// Create new statistics tracker
    pub fn new() -> Self {
        GamdStatistics {
            v_max: f64::NEG_INFINITY,
            v_min: f64::INFINITY,
            v_mean: 0.0,
            sigma_v: 0.0,
            m2: 0.0,
            steps: 0,
        }
    }

    /// Update statistics with new energy value using Welford's algorithm
    ///
    /// This is a numerically stable online algorithm for computing
    /// mean and variance without storing all samples.
    pub fn update(&mut self, energy: f64) {
        self.steps += 1;

        // Update max/min
        if energy > self.v_max {
            self.v_max = energy;
        }
        if energy < self.v_min {
            self.v_min = energy;
        }

        // Welford's online algorithm for mean and variance
        let delta = energy - self.v_mean;
        self.v_mean += delta / self.steps as f64;
        let delta2 = energy - self.v_mean;
        self.m2 += delta * delta2;

        // Update standard deviation
        if self.steps > 1 {
            self.sigma_v = (self.m2 / self.steps as f64).sqrt();
        }
    }

    /// Reset statistics
    pub fn reset(&mut self) {
        self.v_max = f64::NEG_INFINITY;
        self.v_min = f64::INFINITY;
        self.v_mean = 0.0;
        self.sigma_v = 0.0;
        self.m2 = 0.0;
        self.steps = 0;
    }

    /// Check if statistics are ready for parameter calculation
    pub fn is_ready(&self) -> bool {
        self.steps > 10 && self.v_max > self.v_min && self.sigma_v > 0.0
    }
}

impl Default for GamdStatistics {
    fn default() -> Self {
        Self::new()
    }
}

/// GaMD acceleration parameters
#[derive(Debug, Clone)]
pub struct GamdParameters {
    /// Search mode
    pub search_mode: SearchMode,

    /// Boost form
    pub boost_form: BoostForm,

    /// Threshold type
    pub threshold_type: ThresholdType,

    /// Dihedral statistics
    pub dih_stats: GamdStatistics,

    /// Total potential statistics
    pub tot_stats: GamdStatistics,

    /// Dihedral: maximum standard deviation (sigma0)
    pub sigma0_dih: f64,

    /// Total potential: maximum standard deviation (sigma0)
    pub sigma0_tot: f64,

    /// Dihedral: k0 parameter (0 < k0 <= 1)
    pub k0_dih: f64,

    /// Total: k0 parameter (0 < k0 <= 1)
    pub k0_tot: f64,

    /// Dihedral: force constant K
    pub k_dih: f64,

    /// Total: force constant K
    pub k_tot: f64,

    /// Dihedral: energy threshold E
    pub e_dih: f64,

    /// Total: energy threshold E
    pub e_tot: f64,

    /// Current boost potential (for output)
    pub boost_potential: f64,

    /// Equilibration steps (statistics reset after equilibration)
    pub equilibration_steps: usize,

    /// Statistics window size (0 = use all data, >0 = sliding window)
    pub window_size: usize,

    /// Current step counter
    pub current_step: usize,
}

impl GamdParameters {
    /// Create new GaMD parameters
    pub fn new(
        search_mode: SearchMode,
        boost_form: BoostForm,
        threshold_type: ThresholdType,
        sigma0_dih: f64,
        sigma0_tot: f64,
    ) -> Self {
        GamdParameters {
            search_mode,
            boost_form,
            threshold_type,
            dih_stats: GamdStatistics::new(),
            tot_stats: GamdStatistics::new(),
            sigma0_dih,
            sigma0_tot,
            k0_dih: 0.0,
            k0_tot: 0.0,
            k_dih: 0.0,
            k_tot: 0.0,
            e_dih: 0.0,
            e_tot: 0.0,
            boost_potential: 0.0,
            equilibration_steps: 0,
            window_size: 0,
            current_step: 0,
        }
    }

    /// Calculate GaMD parameters from statistics
    ///
    /// Returns (k0, k, E) for the given statistics and sigma0
    fn calculate_parameters(
        &self,
        stats: &GamdStatistics,
        sigma0: f64,
    ) -> Option<(f64, f64, f64)> {
        if !stats.is_ready() {
            return None;
        }

        let v_max = stats.v_max;
        let v_min = stats.v_min;
        let v_mean = stats.v_mean;
        let sigma_v = stats.sigma_v;

        // Avoid division by zero
        if (v_max - v_min).abs() < 1e-10 || sigma_v < 1e-10 {
            return None;
        }

        let (k0, e) = match self.threshold_type {
            ThresholdType::LowerBound => {
                // Lower bound: E = Vmax
                // k0 = (σ₀/σᵥ) × (Vmax - Vmin) / (Vmax - Vmean)
                let k0 = (sigma0 / sigma_v) * (v_max - v_min) / (v_max - v_mean);
                let k0 = k0.min(1.0); // Clamp to 1.0
                (k0, v_max)
            }
            ThresholdType::UpperBound => {
                // Upper bound: E = Vmin + (Vmax - Vmin) / k0
                // k0 = (1 - σ₀/σᵥ) × (Vmax - Vmin) / (Vmean - Vmin)
                let mut k0 = (1.0 - sigma0 / sigma_v) * (v_max - v_min) / (v_mean - v_min);

                // Check if k0 is valid (0 < k0 <= 1)
                if k0 <= 0.0 || k0 > 1.0 {
                    // Fall back to lower bound
                    k0 = (sigma0 / sigma_v) * (v_max - v_min) / (v_max - v_mean);
                    k0 = k0.min(1.0);
                    (k0, v_max)
                } else {
                    let e = v_min + (v_max - v_min) / k0;
                    (k0, e)
                }
            }
        };

        // Calculate force constant K = k0 / (Vmax - Vmin)
        let k = k0 / (v_max - v_min);

        Some((k0, k, e))
    }

    /// Update dihedral parameters
    pub fn update_dihedral_parameters(&mut self) {
        if let Some((k0, k, e)) = self.calculate_parameters(&self.dih_stats, self.sigma0_dih) {
            self.k0_dih = k0;
            self.k_dih = k;
            self.e_dih = e;
        }
    }

    /// Update total potential parameters
    pub fn update_total_parameters(&mut self) {
        if let Some((k0, k, e)) = self.calculate_parameters(&self.tot_stats, self.sigma0_tot) {
            self.k0_tot = k0;
            self.k_tot = k;
            self.e_tot = e;
        }
    }

    /// Update all parameters
    pub fn update_all_parameters(&mut self) {
        self.update_dihedral_parameters();
        self.update_total_parameters();
    }

    /// Update statistics with current energies
    pub fn update_statistics(&mut self, dihedral_energy: f64, total_energy: f64) {
        self.current_step += 1;

        // Check for equilibration reset
        if self.equilibration_steps > 0 && self.current_step == self.equilibration_steps {
            self.dih_stats.reset();
            self.tot_stats.reset();
        }

        // Check for window reset
        if self.window_size > 0 && self.dih_stats.steps % self.window_size == 0 {
            self.dih_stats.reset();
            self.tot_stats.reset();
        }

        // Update statistics
        if self.boost_form == BoostForm::DihedralBoost || self.boost_form == BoostForm::DualBoost {
            self.dih_stats.update(dihedral_energy);
        }

        if self.boost_form == BoostForm::TotalBoost || self.boost_form == BoostForm::DualBoost {
            let non_dihedral_energy = total_energy - dihedral_energy;
            self.tot_stats.update(non_dihedral_energy);
        }

        // For TotalBoost, use the full potential
        if self.boost_form == BoostForm::TotalBoost {
            self.tot_stats.update(total_energy);
        }
    }

    /// Apply GaMD boost to forces
    ///
    /// Returns the total boost potential added
    pub fn apply_boost(
        &mut self,
        configuration: &mut Configuration,
        dihedral_energy: f64,
        dihedral_forces: &[Vec3],
        total_energy: f64,
        total_forces: &[Vec3],
    ) -> f64 {
        self.boost_potential = 0.0;

        // Only apply boost in GaMD search or production modes
        if self.search_mode == SearchMode::CmdSearch {
            return 0.0;
        }

        let n_atoms = configuration.current().pos.len();

        match self.boost_form {
            BoostForm::DihedralBoost => {
                // Boost only dihedral term
                let v_e = self.e_dih - dihedral_energy;
                if v_e > 0.0 {
                    let prefactor = self.k_dih * v_e;
                    self.boost_potential = 0.5 * prefactor * v_e;

                    // Apply boost: F_boost = F_orig * k(E - V)
                    for i in 0..n_atoms {
                        configuration.current_mut().force[i] += dihedral_forces[i] * prefactor as f32;
                    }
                }
            }
            BoostForm::TotalBoost => {
                // Boost entire potential (including dihedrals)
                let v_e = self.e_tot - total_energy;
                if v_e > 0.0 {
                    let prefactor = self.k_tot * v_e;
                    self.boost_potential = 0.5 * prefactor * v_e;

                    // Apply boost to all forces
                    for i in 0..n_atoms {
                        configuration.current_mut().force[i] +=
                            (dihedral_forces[i] + total_forces[i]) * prefactor as f32;
                    }
                }
            }
            BoostForm::DualBoost => {
                // Boost dihedral and non-dihedral independently
                let non_dihedral_energy = total_energy - dihedral_energy;

                // Non-dihedral boost
                let v_e_tot = self.e_tot - non_dihedral_energy;
                if v_e_tot > 0.0 {
                    let prefactor_tot = self.k_tot * v_e_tot;
                    self.boost_potential += 0.5 * prefactor_tot * v_e_tot;

                    for i in 0..n_atoms {
                        configuration.current_mut().force[i] += total_forces[i] * prefactor_tot as f32;
                    }
                }

                // Dihedral boost
                let v_e_dih = self.e_dih - dihedral_energy;
                if v_e_dih > 0.0 {
                    let prefactor_dih = self.k_dih * v_e_dih;
                    self.boost_potential += 0.5 * prefactor_dih * v_e_dih;

                    for i in 0..n_atoms {
                        configuration.current_mut().force[i] += dihedral_forces[i] * prefactor_dih as f32;
                    }
                }
            }
        }

        self.boost_potential
    }
}

/// GaMD MD Runner
///
/// Manages GaMD molecular dynamics simulation with separate
/// tracking of dihedral and non-dihedral forces
pub struct GamdRunner {
    /// GaMD parameters
    pub gamd: GamdParameters,

    /// Pairlist for nonbonded interactions
    pairlist: PairlistContainer,

    /// Pairlist algorithm
    pairlist_algorithm: StandardPairlistAlgorithm,

    /// Use parallel force calculation
    parallel_forces: bool,

    /// CRF parameters for nonbonded interactions
    crf_params: CRFParameters,

    /// Current MD step
    current_step: usize,

    /// Stored dihedral forces
    dihedral_forces: Vec<Vec3>,

    /// Stored non-dihedral forces
    other_forces: Vec<Vec3>,

    /// Dihedral energy
    dihedral_energy: f64,

    /// Non-dihedral energy
    other_energy: f64,
}

impl GamdRunner {
    /// Create a new GaMD runner
    pub fn new(
        gamd: GamdParameters,
        cutoff: f64,
        parallel_forces: bool,
    ) -> Self {
        let pairlist = PairlistContainer::new(cutoff, cutoff, 0.2);
        let pairlist_algorithm = StandardPairlistAlgorithm::new(true);

        // Calculate CRF parameters
        let eps = 1.0;
        let eps_rf = 61.0;
        let crf = (eps_rf - eps) * (2.0 * eps + eps_rf) / (eps_rf + 2.0 * eps) / (cutoff * cutoff * cutoff);

        let crf_params = CRFParameters {
            crf_cut: cutoff,
            crf_2cut3i: crf / (2.0 * cutoff * cutoff * cutoff),
            crf_cut3i: (1.0 - crf / 2.0) / cutoff,
        };

        GamdRunner {
            gamd,
            pairlist,
            pairlist_algorithm,
            parallel_forces,
            crf_params,
            current_step: 0,
            dihedral_forces: Vec::new(),
            other_forces: Vec::new(),
            dihedral_energy: 0.0,
            other_energy: 0.0,
        }
    }

    /// Update pairlist if needed
    fn update_pairlist_if_needed(&mut self, topology: &Topology, configuration: &Configuration) {
        if self.pairlist.needs_update() {
            let box_dims = configuration.current().box_config.dimensions();
            let periodicity = Rectangular::new(box_dims);

            self.pairlist_algorithm.update(
                topology,
                configuration,
                &mut self.pairlist,
                &periodicity,
            );
            self.pairlist.reset_counter();
        } else {
            self.pairlist.step();
        }
    }

    /// Calculate all forces and energies
    fn calculate_forces(
        &mut self,
        topology: &Topology,
        configuration: &mut Configuration,
    ) {
        let n_atoms = topology.num_atoms();

        // Initialize force storage
        self.dihedral_forces = vec![Vec3::ZERO; n_atoms];
        self.other_forces = vec![Vec3::ZERO; n_atoms];
        self.dihedral_energy = 0.0;
        self.other_energy = 0.0;

        // Update pairlist
        self.update_pairlist_if_needed(topology, configuration);

        // Clear configuration forces
        for force in configuration.current_mut().force.iter_mut() {
            *force = Vec3::ZERO;
        }

        // Calculate bonded forces (bonds, angles, dihedrals)
        let bonded_result = calculate_bonded_forces(topology, configuration, true);

        // For simplicity, assume all bonded energy is non-dihedral
        // (In a full implementation, you'd separate dihedral contributions)
        // TODO: Separate dihedral forces from bonded calculation
        for (i, &force) in bonded_result.forces.iter().enumerate() {
            self.other_forces[i] += force;
        }
        self.other_energy += bonded_result.energy;

        // Calculate nonbonded forces
        let box_dims = configuration.current().box_config.dimensions();
        let periodicity = Rectangular::new(box_dims);

        // Convert LJ parameters
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

        // Convert pairlist to (u32, u32) format
        let pairlist_short: Vec<(u32, u32)> = self.pairlist.solute_short.iter()
            .map(|&(i, j)| (i as u32, j as u32))
            .collect();

        // Convert charge and iac
        let charges_f32: Vec<f32> = topology.charge.iter().map(|&q| q as f32).collect();
        let iac_u32: Vec<u32> = topology.iac.iter().map(|&i| i as u32).collect();

        let nonbonded_storage = if self.parallel_forces {
            lj_crf_innerloop_parallel(
                &configuration.current().pos,
                &charges_f32,
                &iac_u32,
                &pairlist_short,
                &lj_params,
                &self.crf_params,
                &periodicity,
                n_atoms,
            )
        } else {
            let mut storage = ForceStorage::new(n_atoms);
            lj_crf_innerloop(
                &configuration.current().pos,
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

        // Add nonbonded forces to other_forces
        for (i, &force) in nonbonded_storage.forces.iter().enumerate() {
            self.other_forces[i] += force;
        }
        self.other_energy += nonbonded_storage.e_lj + nonbonded_storage.e_crf;

        // Copy forces to configuration (will be modified by GaMD boost)
        for i in 0..n_atoms {
            configuration.current_mut().force[i] = self.dihedral_forces[i] + self.other_forces[i];
        }

        // Update configuration energies
        configuration.current_mut().energies.bond_total = bonded_result.energy;
        configuration.current_mut().energies.lj_total = nonbonded_storage.e_lj;
        configuration.current_mut().energies.crf_total = nonbonded_storage.e_crf;
        configuration.current_mut().energies.potential_total = self.dihedral_energy + self.other_energy;
    }

    /// Perform one GaMD MD step
    pub fn md_step<I: Integrator>(
        &mut self,
        topology: &Topology,
        configuration: &mut Configuration,
        integrator: &mut I,
        dt: f64,
        thermostat_params: Option<&BerendsenThermostatParameters>,
        shake_params: Option<&ShakeParameters>,
    ) {
        // Calculate forces
        self.calculate_forces(topology, configuration);

        let total_energy = self.dihedral_energy + self.other_energy;

        // Update statistics (in all search modes)
        self.gamd.update_statistics(self.dihedral_energy, total_energy);

        // Update parameters (in GaMD search mode)
        if self.gamd.search_mode == SearchMode::GamdSearch {
            self.gamd.update_all_parameters();
        }

        // Apply GaMD boost
        let boost = self.gamd.apply_boost(
            configuration,
            self.dihedral_energy,
            &self.dihedral_forces,
            total_energy,
            &self.other_forces,
        );

        // Add boost to potential energy
        configuration.current_mut().energies.potential_total += boost;

        // Integrate
        integrator.step(dt, topology, configuration);

        // Apply constraints (SHAKE)
        if let Some(shake_p) = shake_params {
            shake(topology, configuration, dt, shake_p);
        }

        // Calculate kinetic energy
        configuration.current_mut().calculate_kinetic_energy(&topology.mass);

        // Apply thermostat
        if let Some(params) = thermostat_params {
            berendsen_thermostat(topology, configuration, dt, params);
        }

        self.current_step += 1;
    }

    /// Run multiple GaMD MD steps
    pub fn run_md_steps<I: Integrator>(
        &mut self,
        topology: &Topology,
        configuration: &mut Configuration,
        integrator: &mut I,
        n_steps: usize,
        dt: f64,
        thermostat_params: Option<&BerendsenThermostatParameters>,
        shake_params: Option<&ShakeParameters>,
    ) {
        for _ in 0..n_steps {
            self.md_step(
                topology,
                configuration,
                integrator,
                dt,
                thermostat_params,
                shake_params,
            );
        }
    }

    /// Get current boost potential
    pub fn boost_potential(&self) -> f64 {
        self.gamd.boost_potential
    }

    /// Get dihedral statistics
    pub fn dih_stats(&self) -> &GamdStatistics {
        &self.gamd.dih_stats
    }

    /// Get total potential statistics
    pub fn tot_stats(&self) -> &GamdStatistics {
        &self.gamd.tot_stats
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_statistics_update() {
        let mut stats = GamdStatistics::new();

        // Add some samples
        stats.update(100.0);
        stats.update(110.0);
        stats.update(90.0);
        stats.update(105.0);

        assert_eq!(stats.steps, 4);
        assert_eq!(stats.v_max, 110.0);
        assert_eq!(stats.v_min, 90.0);
        assert!((stats.v_mean - 101.25).abs() < 1e-10);
        assert!(stats.sigma_v > 0.0);
    }

    #[test]
    fn test_parameter_calculation_lower_bound() {
        let params = GamdParameters::new(
            SearchMode::GamdSearch,
            BoostForm::TotalBoost,
            ThresholdType::LowerBound,
            1.0,
            1.0,
        );

        let mut stats = GamdStatistics::new();
        for i in 0..100 {
            stats.update(-100.0 + i as f64);
        }

        if let Some((k0, k, e)) = params.calculate_parameters(&stats, 1.0) {
            assert!(k0 > 0.0 && k0 <= 1.0);
            assert!(k > 0.0);
            assert!((e - stats.v_max).abs() < 1e-10); // Lower bound: E = Vmax
        } else {
            panic!("Parameter calculation failed");
        }
    }

    #[test]
    fn test_parameter_calculation_upper_bound() {
        let params = GamdParameters::new(
            SearchMode::GamdSearch,
            BoostForm::TotalBoost,
            ThresholdType::UpperBound,
            1.0,
            1.0,
        );

        let mut stats = GamdStatistics::new();
        for i in 0..100 {
            stats.update(-100.0 + i as f64);
        }

        if let Some((k0, k, e)) = params.calculate_parameters(&stats, 1.0) {
            assert!(k0 > 0.0 && k0 <= 1.0);
            assert!(k > 0.0);
            assert!(e > stats.v_min && e < stats.v_max); // Upper bound: Vmin < E < Vmax
        } else {
            panic!("Parameter calculation failed");
        }
    }

    #[test]
    fn test_boost_potential() {
        let mut params = GamdParameters::new(
            SearchMode::GamdSearch,
            BoostForm::TotalBoost,
            ThresholdType::LowerBound,
            1.0,
            1.0,
        );

        // Set up parameters manually
        params.k_tot = 0.01;
        params.e_tot = 100.0;

        // Create dummy configuration
        let mut config = Configuration::new(10, 3, 1);
        for i in 0..10 {
            config.current_mut().force[i] = Vec3::new(1.0, 0.0, 0.0);
        }

        let dih_forces = vec![Vec3::ZERO; 10];
        let tot_forces = vec![Vec3::new(1.0, 0.0, 0.0); 10];

        // Apply boost with V = 90 < E = 100
        let boost = params.apply_boost(&mut config, 0.0, &dih_forces, 90.0, &tot_forces);

        // Check boost potential: ΔV = 0.5 * k * (E - V)^2
        let expected = 0.5 * 0.01 * (100.0 - 90.0) * (100.0 - 90.0);
        assert!((boost - expected).abs() < 1e-10);

        // Check forces were modified
        let prefactor = 0.01 * (100.0 - 90.0);
        let expected_force = 1.0 + 1.0 * prefactor as f32;
        assert!((config.current().force[0].x - expected_force).abs() < 1e-5);
    }

    #[test]
    fn test_cmd_search_no_boost() {
        let mut params = GamdParameters::new(
            SearchMode::CmdSearch,
            BoostForm::TotalBoost,
            ThresholdType::LowerBound,
            1.0,
            1.0,
        );

        params.k_tot = 0.01;
        params.e_tot = 100.0;

        let mut config = Configuration::new(10, 3, 1);
        let dih_forces = vec![Vec3::ZERO; 10];
        let tot_forces = vec![Vec3::new(1.0, 0.0, 0.0); 10];

        // CMD search should not apply boost
        let boost = params.apply_boost(&mut config, 0.0, &dih_forces, 90.0, &tot_forces);
        assert_eq!(boost, 0.0);
    }
}
