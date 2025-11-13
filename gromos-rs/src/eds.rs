/// Enveloping Distribution Sampling (EDS)
///
/// This module implements EDS for enhanced sampling across multiple potential energy surfaces.
///
/// ## Theory
///
/// EDS uses a smooth reference Hamiltonian that envelopes multiple end-state potentials:
///
/// ```text
/// V_R(R) = -(1/βs) * ln[Σ_i exp(-βs(V_i(R) - E_i^R))]
/// ```
///
/// Where:
/// - `V_R`: Reference state energy (smooth Hamiltonian)
/// - `β = 1/(k_B·T)`: Inverse temperature
/// - `s`: Smoothing parameter (0 < s ≤ 1)
/// - `V_i`: Potential energy of state i
/// - `E_i^R`: Energy offset for state i
///
/// ## References
/// - Hansen & Hünenberger, J. Chem. Theory Comput. 2014, 10, 2632-2647
/// - Christ & van Gunsteren, J. Chem. Phys. 2007, 126, 184110

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

/// Boltzmann constant in kJ/(mol·K)
const KB: f64 = 0.008314462618;

/// EDS functional form
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EDSForm {
    /// Single smoothing parameter for all state pairs
    SingleS,
    /// Pairwise smoothing parameters (N(N-1)/2 values)
    MultiS,
    /// Specified pairs with dedicated s values
    PairS,
}

/// EDS state information
#[derive(Debug, Clone)]
pub struct EDSState {
    /// State ID
    pub id: usize,

    /// Individual state energy V_i
    pub energy: f64,

    /// Energy offset E_i^R
    pub offset: f64,

    /// Forces for this state
    pub forces: Vec<Vec3>,

    /// Virial tensor for this state (for pressure calculation)
    pub virial: [[f64; 3]; 3],

    /// Number of times this state was visited
    pub visit_count: usize,

    /// Was this state visited in current round trip?
    pub visited: bool,
}

impl EDSState {
    /// Create a new EDS state
    pub fn new(id: usize, offset: f64, n_atoms: usize) -> Self {
        EDSState {
            id,
            energy: 0.0,
            offset,
            forces: vec![Vec3::ZERO; n_atoms],
            virial: [[0.0; 3]; 3],
            visit_count: 0,
            visited: false,
        }
    }

    /// Clear forces and virial for next step
    pub fn clear(&mut self) {
        for f in &mut self.forces {
            *f = Vec3::ZERO;
        }
        self.virial = [[0.0; 3]; 3];
        self.energy = 0.0;
    }

    /// Get adjusted energy (V_i - E_i^R)
    pub fn adjusted_energy(&self) -> f64 {
        self.energy - self.offset
    }
}

/// EDS parameters and state
pub struct EDSParameters {
    /// Functional form
    pub form: EDSForm,

    /// Number of states
    pub num_states: usize,

    /// Smoothing parameter(s)
    /// - SingleS: 1 value
    /// - MultiS: N(N-1)/2 values
    /// - PairS: variable
    pub s_values: Vec<f64>,

    /// Temperature (K)
    pub temperature: f64,

    /// All EDS states
    pub states: Vec<EDSState>,

    /// Current reference state energy V_R
    pub reference_energy: f64,

    /// Mixed state energy (before acceleration)
    pub mixed_energy: f64,

    /// Current state ID (state with minimum adjusted energy)
    pub current_state_id: usize,

    /// Round-trip counter
    pub round_trips: usize,
}

impl EDSParameters {
    /// Create new EDS parameters
    pub fn new(
        form: EDSForm,
        s_values: Vec<f64>,
        offsets: Vec<f64>,
        temperature: f64,
        n_atoms: usize,
    ) -> Result<Self, String> {
        let num_states = offsets.len();

        if num_states < 2 {
            return Err("EDS requires at least 2 states".to_string());
        }

        // Validate s_values
        match form {
            EDSForm::SingleS => {
                if s_values.len() != 1 {
                    return Err("SingleS requires exactly 1 s value".to_string());
                }
            }
            EDSForm::MultiS => {
                let expected = (num_states * (num_states - 1)) / 2;
                if s_values.len() != expected {
                    return Err(format!(
                        "MultiS requires {} s values for {} states",
                        expected, num_states
                    ));
                }
            }
            EDSForm::PairS => {
                // Variable number, depends on specified pairs
            }
        }

        // Create states
        let states: Vec<EDSState> = offsets
            .iter()
            .enumerate()
            .map(|(id, &offset)| EDSState::new(id, offset, n_atoms))
            .collect();

        Ok(EDSParameters {
            form,
            num_states,
            s_values,
            temperature,
            states,
            reference_energy: 0.0,
            mixed_energy: 0.0,
            current_state_id: 0,
            round_trips: 0,
        })
    }

    /// Get beta = 1/(k_B * T)
    pub fn beta(&self) -> f64 {
        1.0 / (KB * self.temperature)
    }

    /// Calculate reference energy using EDS Hamiltonian (single-s form)
    ///
    /// Uses numerically stable log-sum-exp trick:
    /// log(exp(a) + exp(b)) = max(a,b) + log(1 + exp(min(a,b) - max(a,b)))
    pub fn calculate_reference_energy_single_s(&mut self) {
        let beta = self.beta();
        let s = self.s_values[0];

        if self.num_states == 0 {
            self.reference_energy = 0.0;
            return;
        }

        // Calculate -β*s*(V_i - E_i^R) for each state
        let mut prefactors: Vec<f64> = self.states
            .iter()
            .map(|state| -beta * s * (state.energy - state.offset))
            .collect();

        // Log-sum-exp with numerical stability
        let mut sum_prefactors = prefactors[0];

        for &pf in &prefactors[1..] {
            sum_prefactors = sum_prefactors.max(pf)
                + (1.0 + ((sum_prefactors.min(pf) - sum_prefactors.max(pf)).exp())).ln();
        }

        // V_R = -(1/(β*s)) * log-sum-exp
        self.mixed_energy = -1.0 / (beta * s) * sum_prefactors;
        self.reference_energy = self.mixed_energy;
    }

    /// Calculate reference energy using multi-s form
    ///
    /// For multi-s, we use pairwise smoothing parameters s_ij
    /// The Hamiltonian becomes more complex with N(N-1)/2 parameters
    pub fn calculate_reference_energy_multi_s(&mut self) {
        let beta = self.beta();

        if self.num_states == 0 {
            self.reference_energy = 0.0;
            return;
        }

        // For multi-s, we calculate a weighted average based on pairwise s values
        // This is a simplified implementation - full multi-s is more complex
        // Using average s for now
        let avg_s: f64 = self.s_values.iter().sum::<f64>() / self.s_values.len() as f64;

        // Calculate using average s (similar to single-s)
        let prefactors: Vec<f64> = self.states
            .iter()
            .map(|state| -beta * avg_s * (state.energy - state.offset))
            .collect();

        let max_prefactor = prefactors.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let sum_exp: f64 = prefactors
            .iter()
            .map(|&pf| (pf - max_prefactor).exp())
            .sum();

        self.mixed_energy = -1.0 / (beta * avg_s) * (max_prefactor + sum_exp.ln());
        self.reference_energy = self.mixed_energy;
    }

    /// Calculate reference energy based on form
    pub fn calculate_reference_energy(&mut self) {
        match self.form {
            EDSForm::SingleS => self.calculate_reference_energy_single_s(),
            EDSForm::MultiS => self.calculate_reference_energy_multi_s(),
            EDSForm::PairS => self.calculate_reference_energy_single_s(), // Use single-s for now
        }
    }

    /// Get smoothing parameter for a specific state pair (for multi-s/pair-s)
    fn get_s_for_pair(&self, i: usize, j: usize) -> f64 {
        match self.form {
            EDSForm::SingleS => self.s_values[0],
            EDSForm::MultiS => {
                // Map (i,j) to index in s_values array
                // For i < j: index = i*n - i*(i+1)/2 + (j-i-1)
                if i < j {
                    let n = self.num_states;
                    let idx = i * n - i * (i + 1) / 2 + (j - i - 1);
                    self.s_values[idx]
                } else if i > j {
                    let n = self.num_states;
                    let idx = j * n - j * (j + 1) / 2 + (i - j - 1);
                    self.s_values[idx]
                } else {
                    self.s_values[0] // Same state, use first s
                }
            }
            EDSForm::PairS => {
                // For pair-s, we'd need a mapping structure
                // For now, use first s value
                self.s_values[0]
            }
        }
    }

    /// Calculate state probabilities
    pub fn calculate_state_probabilities(&self) -> Vec<f64> {
        let beta = self.beta();
        let s = self.s_values[0];

        // Calculate prefactors
        let prefactors: Vec<f64> = self.states
            .iter()
            .map(|state| -beta * s * (state.energy - state.offset))
            .collect();

        // Find max for numerical stability
        let max_prefactor = prefactors.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        // Calculate sum of exp(prefactor - max)
        let sum_exp: f64 = prefactors
            .iter()
            .map(|&pf| (pf - max_prefactor).exp())
            .sum();

        // Calculate probabilities
        prefactors
            .iter()
            .map(|&pf| (pf - max_prefactor).exp() / sum_exp)
            .collect()
    }

    /// Apply EDS forces to configuration
    ///
    /// F = Σ_i p_i * F_i
    /// where p_i = exp(-βs(V_i - E_i^R)) / Z
    pub fn apply_forces(&self, configuration: &mut Configuration) {
        let probabilities = self.calculate_state_probabilities();

        // Clear current forces
        for force in configuration.current_mut().force.iter_mut() {
            *force = Vec3::ZERO;
        }

        // Apply weighted forces from all states
        for (state, &prob) in self.states.iter().zip(probabilities.iter()) {
            for (i, &state_force) in state.forces.iter().enumerate() {
                configuration.current_mut().force[i] += state_force * prob as f32;
            }
        }
    }

    /// Apply EDS virial to configuration
    pub fn apply_virial(&self, configuration: &mut Configuration) {
        let probabilities = self.calculate_state_probabilities();

        // Accumulate weighted virial from all states
        let mut virial_sum = [[0.0_f64; 3]; 3];

        for (state, &prob) in self.states.iter().zip(probabilities.iter()) {
            for i in 0..3 {
                for j in 0..3 {
                    virial_sum[i][j] += state.virial[i][j] * prob;
                }
            }
        }

        // Convert to Mat3 and assign (Mat3 is column-major)
        configuration.current_mut().virial_tensor = crate::math::Mat3::from_cols(
            Vec3::new(virial_sum[0][0] as f32, virial_sum[1][0] as f32, virial_sum[2][0] as f32),
            Vec3::new(virial_sum[0][1] as f32, virial_sum[1][1] as f32, virial_sum[2][1] as f32),
            Vec3::new(virial_sum[0][2] as f32, virial_sum[1][2] as f32, virial_sum[2][2] as f32),
        );
    }

    /// Determine current state (state with minimum adjusted energy)
    pub fn determine_current_state(&mut self) {
        let mut min_energy = f64::INFINITY;
        let mut min_state = 0;

        for (i, state) in self.states.iter().enumerate() {
            let adjusted = state.adjusted_energy();
            if adjusted < min_energy {
                min_energy = adjusted;
                min_state = i;
            }
        }

        // Mark state as visited
        self.states[min_state].visited = true;
        self.states[min_state].visit_count += 1;
        self.current_state_id = min_state;
    }

    /// Check if a round trip is complete (all states visited)
    pub fn check_round_trip(&mut self) -> bool {
        let all_visited = self.states.iter().all(|s| s.visited);

        if all_visited {
            self.round_trips += 1;

            // Reset visited flags
            for state in &mut self.states {
                state.visited = false;
            }

            true
        } else {
            false
        }
    }

    /// Clear all state forces and energies
    pub fn clear_states(&mut self) {
        for state in &mut self.states {
            state.clear();
        }
    }
}

/// AEDS (Accelerated EDS) parameters
pub struct AEDSParameters {
    /// Base EDS parameters
    pub eds: EDSParameters,

    /// Maximum transition energy (acceleration region upper bound)
    pub e_max: f64,

    /// Minimum transition energy (acceleration region lower bound)
    pub e_min: f64,

    /// Average energy per state
    pub avg_energy: Vec<f64>,

    /// Standard deviation of energy per state
    pub std_energy: Vec<f64>,

    /// Free energy differences relative to state 0
    pub free_energy: Vec<f64>,

    /// Enable parameter search
    pub search_enabled: bool,
}

impl AEDSParameters {
    /// Create new AEDS parameters
    pub fn new(
        eds: EDSParameters,
        e_max: f64,
        e_min: f64,
        search_enabled: bool,
    ) -> Self {
        let n = eds.num_states;

        AEDSParameters {
            eds,
            e_max,
            e_min,
            avg_energy: vec![0.0; n],
            std_energy: vec![0.0; n],
            free_energy: vec![0.0; n],
            search_enabled,
        }
    }

    /// Calculate acceleration factor
    ///
    /// fkfac = 1 - (V_mix - E_min) / (E_max - E_min)
    /// fkfac = 0 at E_max, fkfac = 1 at E_min
    pub fn acceleration_factor(&self) -> f64 {
        let v_mix = self.eds.mixed_energy;

        if v_mix <= self.e_min {
            return 1.0; // No acceleration below E_min
        } else if v_mix >= self.e_max {
            return 1.0; // No acceleration above E_max
        } else {
            // Linear interpolation in acceleration region
            let de_mix = v_mix - self.e_min;
            let kfac = 1.0 / (self.e_max - self.e_min);
            1.0 - kfac * de_mix
        }
    }

    /// Apply AEDS forces (with acceleration)
    pub fn apply_forces(&self, configuration: &mut Configuration) {
        let factor = self.acceleration_factor();

        // Get base EDS forces
        self.eds.apply_forces(configuration);

        // Scale by acceleration factor if in acceleration region
        if factor < 1.0 {
            for force in configuration.current_mut().force.iter_mut() {
                *force = *force * factor as f32;
            }
        }
    }

    /// Update energy statistics for AEDS parameter search
    ///
    /// Accumulates running statistics for each state to estimate
    /// average energies and standard deviations
    pub fn update_statistics(&mut self, alpha: f64) {
        // Exponential moving average with decay factor alpha
        // E_avg(t+1) = alpha * E(t) + (1-alpha) * E_avg(t)
        // E_std(t+1) = sqrt(alpha * (E(t) - E_avg(t))^2 + (1-alpha) * E_std(t)^2)

        for (i, state) in self.eds.states.iter().enumerate() {
            let energy = state.energy;

            // Update average
            let old_avg = self.avg_energy[i];
            self.avg_energy[i] = alpha * energy + (1.0 - alpha) * old_avg;

            // Update standard deviation
            let diff = energy - old_avg;
            let old_var = self.std_energy[i] * self.std_energy[i];
            let new_var = alpha * diff * diff + (1.0 - alpha) * old_var;
            self.std_energy[i] = new_var.sqrt();
        }
    }

    /// Estimate free energy differences using thermodynamic integration
    ///
    /// Uses the relationship: ΔF_ij ≈ -kT * ln(<exp(-β(V_j - V_i))>_i)
    /// Simplified version uses energy differences and visit counts
    pub fn estimate_free_energies(&mut self) {
        let beta = self.eds.beta();

        // Calculate free energy relative to state 0
        self.free_energy[0] = 0.0;

        for i in 1..self.eds.num_states {
            // Simple estimate based on average energy difference
            let de = self.avg_energy[i] - self.avg_energy[0];

            // Free energy difference (simplified)
            // ΔF ≈ ΔE - T*ΔS, approximating ΔS from visit ratio
            let visit_ratio = if self.eds.states[0].visit_count > 0 {
                self.eds.states[i].visit_count as f64 / self.eds.states[0].visit_count as f64
            } else {
                1.0
            };

            // Entropy contribution from visit counts
            let ds = if visit_ratio > 0.0 {
                -visit_ratio.ln()
            } else {
                0.0
            };

            self.free_energy[i] = de - ds / beta;
        }
    }

    /// Optimize energy offsets to improve sampling efficiency
    ///
    /// Adjusts E_i^R values to flatten the free energy landscape
    /// Goal: make all states equally probable
    pub fn optimize_offsets(&mut self, learning_rate: f64) {
        // Calculate current state probabilities
        let probabilities = self.eds.calculate_state_probabilities();

        // Target probability (uniform distribution)
        let target_prob = 1.0 / self.eds.num_states as f64;

        // Update offsets based on probability deviation
        for i in 0..self.eds.num_states {
            let prob_error = probabilities[i] - target_prob;

            // Adjust offset: if state is over-visited, increase offset to disfavor it
            // ΔE_i^R = -η * ln(p_i / p_target)
            let adjustment = if probabilities[i] > 1e-10 {
                -learning_rate * (probabilities[i] / target_prob).ln()
            } else {
                0.0
            };

            self.eds.states[i].offset += adjustment;
        }
    }

    /// Optimize smoothing parameter to balance exploration/exploitation
    ///
    /// Adjusts s based on round-trip efficiency and energy overlap
    pub fn optimize_smoothing(&mut self, learning_rate: f64) {
        if self.eds.form != EDSForm::SingleS {
            return; // Only optimize single-s for now
        }

        // Calculate energy overlap between states
        let mut total_overlap = 0.0;
        let mut n_pairs = 0;

        for i in 0..self.eds.num_states {
            for j in (i + 1)..self.eds.num_states {
                // Overlap metric: exp(-|ΔE|/(kT))
                let de = (self.avg_energy[i] - self.avg_energy[j]).abs();
                let beta = self.eds.beta();
                let overlap = (-beta * de).exp();
                total_overlap += overlap;
                n_pairs += 1;
            }
        }

        let avg_overlap = if n_pairs > 0 {
            total_overlap / n_pairs as f64
        } else {
            0.5
        };

        // Target overlap: ~0.2-0.3 for good sampling
        let target_overlap = 0.25;
        let overlap_error = avg_overlap - target_overlap;

        // Adjust s: if overlap too high, increase s (sharper peaks)
        // if overlap too low, decrease s (smoother landscape)
        let adjustment = learning_rate * overlap_error;
        let new_s = (self.eds.s_values[0] + adjustment).clamp(0.1, 1.0);
        self.eds.s_values[0] = new_s;
    }

    /// Full AEDS parameter search iteration
    ///
    /// Performs one iteration of adaptive parameter optimization
    pub fn search_iteration(&mut self, learning_rate: f64, update_alpha: f64) {
        if !self.search_enabled {
            return;
        }

        // Update statistics
        self.update_statistics(update_alpha);

        // Every N steps, optimize parameters
        if self.eds.states[0].visit_count % 100 == 0 {
            self.estimate_free_energies();
            self.optimize_offsets(learning_rate);
            self.optimize_smoothing(learning_rate * 0.5); // More conservative for s
        }
    }
}

/// EDS MD Runner
///
/// Manages the EDS molecular dynamics simulation loop,
/// coordinating force calculations for multiple states
pub struct EDSRunner {
    /// EDS or AEDS parameters
    pub aeds: AEDSParameters,

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
}

impl EDSRunner {
    /// Create a new EDS runner
    pub fn new(
        aeds: AEDSParameters,
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

        EDSRunner {
            aeds,
            pairlist,
            pairlist_algorithm,
            parallel_forces,
            crf_params,
            current_step: 0,
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

    /// Calculate forces for a specific EDS state
    fn calculate_state_forces(
        &self,
        state_id: usize,
        topology: &Topology,
        configuration: &Configuration,
    ) -> (f64, Vec<Vec3>, [[f64; 3]; 3]) {
        let n_atoms = topology.num_atoms();
        let mut forces = vec![Vec3::ZERO; n_atoms];
        let mut virial = [[0.0; 3]; 3];
        let mut total_energy = 0.0;

        // Bonded forces
        let bonded_result = calculate_bonded_forces(topology, configuration, true);
        for (i, &force) in bonded_result.forces.iter().enumerate() {
            forces[i] += force;
        }
        total_energy += bonded_result.energy;

        // Nonbonded forces
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

        // Accumulate nonbonded forces
        for (i, &force) in nonbonded_storage.forces.iter().enumerate() {
            forces[i] += force;
        }
        total_energy += nonbonded_storage.e_lj + nonbonded_storage.e_crf;

        // TODO: Extract virial from nonbonded calculation
        // For now, virial is not properly accumulated

        (total_energy, forces, virial)
    }

    /// Calculate forces for all EDS states
    pub fn calculate_all_state_forces(
        &mut self,
        topology: &Topology,
        configuration: &mut Configuration,
    ) {
        // Update pairlist
        self.update_pairlist_if_needed(topology, configuration);

        // Clear all state forces
        self.aeds.eds.clear_states();

        // Calculate forces for each state
        for state_id in 0..self.aeds.eds.num_states {
            let (energy, forces, virial) = self.calculate_state_forces(
                state_id,
                topology,
                configuration,
            );

            // Store state energy and forces
            self.aeds.eds.states[state_id].energy = energy;
            self.aeds.eds.states[state_id].forces = forces;
            self.aeds.eds.states[state_id].virial = virial;
        }

        // Calculate EDS reference energy
        self.aeds.eds.calculate_reference_energy();

        // Determine current state
        self.aeds.eds.determine_current_state();

        // Check for round trips
        self.aeds.eds.check_round_trip();

        // Apply AEDS forces to configuration
        self.aeds.apply_forces(configuration);

        // Apply AEDS virial
        self.aeds.eds.apply_virial(configuration);

        // Update configuration energy with reference energy
        configuration.current_mut().energies.potential_total = self.aeds.eds.reference_energy;
    }

    /// Perform one EDS MD step
    pub fn md_step<I: Integrator>(
        &mut self,
        topology: &Topology,
        configuration: &mut Configuration,
        integrator: &mut I,
        dt: f64,
        thermostat_params: Option<&BerendsenThermostatParameters>,
        shake_params: Option<&ShakeParameters>,
    ) {
        // Calculate forces for all states and apply EDS forces
        self.calculate_all_state_forces(topology, configuration);

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

        // AEDS parameter search
        if self.aeds.search_enabled {
            self.aeds.search_iteration(0.01, 0.05);
        }

        self.current_step += 1;
    }

    /// Run multiple EDS MD steps
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

    /// Get current EDS state ID
    pub fn current_state(&self) -> usize {
        self.aeds.eds.current_state_id
    }

    /// Get reference energy
    pub fn reference_energy(&self) -> f64 {
        self.aeds.eds.reference_energy
    }

    /// Get number of round trips
    pub fn round_trips(&self) -> usize {
        self.aeds.eds.round_trips
    }

    /// Get state visit counts
    pub fn visit_counts(&self) -> Vec<usize> {
        self.aeds.eds.states.iter().map(|s| s.visit_count).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eds_creation() {
        let eds = EDSParameters::new(
            EDSForm::SingleS,
            vec![0.5],
            vec![0.0, 10.0, 20.0],
            300.0,
            100,
        );

        assert!(eds.is_ok());
        let eds = eds.unwrap();
        assert_eq!(eds.num_states, 3);
        assert_eq!(eds.states.len(), 3);
    }

    #[test]
    fn test_beta_calculation() {
        let eds = EDSParameters::new(
            EDSForm::SingleS,
            vec![0.5],
            vec![0.0, 10.0],
            300.0,
            100,
        ).unwrap();

        let beta = eds.beta();
        // β = 1/(k_B * T) = 1/(0.008314462618 * 300) ≈ 0.4006
        assert!((beta - 0.4006).abs() < 0.001);
    }

    #[test]
    fn test_state_probabilities() {
        let mut eds = EDSParameters::new(
            EDSForm::SingleS,
            vec![0.5],
            vec![0.0, 0.0],
            300.0,
            100,
        ).unwrap();

        // Set equal energies
        eds.states[0].energy = 100.0;
        eds.states[1].energy = 100.0;

        let probs = eds.calculate_state_probabilities();

        // Should be approximately equal
        assert!((probs[0] - 0.5).abs() < 0.01);
        assert!((probs[1] - 0.5).abs() < 0.01);

        // Should sum to 1
        let sum: f64 = probs.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_round_trip_detection() {
        let mut eds = EDSParameters::new(
            EDSForm::SingleS,
            vec![0.5],
            vec![0.0, 10.0, 20.0],
            300.0,
            100,
        ).unwrap();

        // Mark first two states as visited
        eds.states[0].visited = true;
        eds.states[1].visited = true;

        // Not complete yet
        assert!(!eds.check_round_trip());

        // Mark last state
        eds.states[2].visited = true;

        // Now complete
        assert!(eds.check_round_trip());
        assert_eq!(eds.round_trips, 1);

        // Flags should be reset
        assert!(!eds.states[0].visited);
        assert!(!eds.states[1].visited);
        assert!(!eds.states[2].visited);
    }

    #[test]
    fn test_aeds_acceleration_factor() {
        let eds = EDSParameters::new(
            EDSForm::SingleS,
            vec![0.5],
            vec![0.0, 10.0],
            300.0,
            100,
        ).unwrap();

        let mut aeds = AEDSParameters::new(eds, 10.0, -50.0, false);

        // Below E_min: no acceleration
        aeds.eds.mixed_energy = -60.0;
        assert_eq!(aeds.acceleration_factor(), 1.0);

        // Above E_max: no acceleration
        aeds.eds.mixed_energy = 15.0;
        assert_eq!(aeds.acceleration_factor(), 1.0);

        // At E_min: full force
        aeds.eds.mixed_energy = -50.0;
        assert_eq!(aeds.acceleration_factor(), 1.0);

        // In middle: partial acceleration
        aeds.eds.mixed_energy = -20.0; // Halfway between -50 and 10
        let factor = aeds.acceleration_factor();
        assert!(factor > 0.0 && factor < 1.0);
    }
}
