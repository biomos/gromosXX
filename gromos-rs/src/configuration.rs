//! Configuration module - system state and energy management
//!
//! Direct Rust translation of GROMOS configuration structures from:
//! - md++/src/configuration/configuration.h
//! - md++/src/configuration/energy.h

use crate::math::{Vec3, Mat3};

/// Simulation box representation
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoxType {
    Vacuum,
    Rectangular,
    Triclinic,
    TruncatedOctahedral,
}

/// Simulation box with periodic boundary conditions
#[derive(Debug, Clone)]
pub struct Box {
    pub box_type: BoxType,
    pub vectors: Mat3,      // Three column vectors defining the box
    pub inv_vectors: Mat3,  // Inverse for wrapping
}

impl Box {
    /// Create vacuum box (no periodicity)
    pub fn vacuum() -> Self {
        Self {
            box_type: BoxType::Vacuum,
            vectors: Mat3::IDENTITY,
            inv_vectors: Mat3::IDENTITY,
        }
    }

    /// Create rectangular box with dimensions (lx, ly, lz)
    pub fn rectangular(lx: f32, ly: f32, lz: f32) -> Self {
        let vectors = Mat3::from_cols(
            Vec3::new(lx, 0.0, 0.0),
            Vec3::new(0.0, ly, 0.0),
            Vec3::new(0.0, 0.0, lz),
        );

        let inv_vectors = Mat3::from_cols(
            Vec3::new(1.0 / lx, 0.0, 0.0),
            Vec3::new(0.0, 1.0 / ly, 0.0),
            Vec3::new(0.0, 0.0, 1.0 / lz),
        );

        Self {
            box_type: BoxType::Rectangular,
            vectors,
            inv_vectors,
        }
    }

    /// Create triclinic box from box vectors
    pub fn triclinic(vectors: Mat3) -> Self {
        Self {
            box_type: BoxType::Triclinic,
            vectors,
            inv_vectors: vectors.inverse(),
        }
    }

    /// Get box volume
    pub fn volume(&self) -> f64 {
        // Volume = determinant of box matrix
        self.vectors.determinant() as f64
    }

    /// Get box dimensions (for rectangular box)
    pub fn dimensions(&self) -> Vec3 {
        Vec3::new(
            self.vectors.x_axis.x,
            self.vectors.y_axis.y,
            self.vectors.z_axis.z,
        )
    }
}

/// Energy storage for the system
///
/// Stores both total energies and per-group energies for detailed accounting
#[derive(Debug, Clone)]
pub struct Energy {
    // Total energies
    pub kinetic_total: f64,
    pub potential_total: f64,

    // Bonded energies
    pub bond_total: f64,
    pub angle_total: f64,
    pub dihedral_total: f64,
    pub improper_total: f64,
    pub cross_dihedral_total: f64,

    // Nonbonded energies
    pub lj_total: f64,
    pub crf_total: f64,  // Coulomb reaction field
    pub ls_total: f64,   // Lattice sum (Ewald, P3M, etc.)
    pub ls_realspace_total: f64,
    pub ls_kspace_total: f64,

    // Special interactions
    pub special_total: f64,
    pub sasa_total: f64,  // Solvent accessible surface area
    pub constraint_total: f64,

    // Per-group energies (for temperature/energy groups)
    pub kinetic_energy: Vec<f64>,  // One per temperature group

    // Per-group pair energies (energy_group x energy_group matrix)
    pub lj_energy: Vec<Vec<f64>>,
    pub crf_energy: Vec<Vec<f64>>,

    // Virial and pressure
    pub virial_total: f64,
}

impl Energy {
    pub fn new(num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            kinetic_total: 0.0,
            potential_total: 0.0,
            bond_total: 0.0,
            angle_total: 0.0,
            dihedral_total: 0.0,
            improper_total: 0.0,
            cross_dihedral_total: 0.0,
            lj_total: 0.0,
            crf_total: 0.0,
            ls_total: 0.0,
            ls_realspace_total: 0.0,
            ls_kspace_total: 0.0,
            special_total: 0.0,
            sasa_total: 0.0,
            constraint_total: 0.0,
            kinetic_energy: vec![0.0; num_temperature_groups],
            lj_energy: vec![vec![0.0; num_energy_groups]; num_energy_groups],
            crf_energy: vec![vec![0.0; num_energy_groups]; num_energy_groups],
            virial_total: 0.0,
        }
    }

    /// Total energy (kinetic + potential)
    pub fn total(&self) -> f64 {
        self.kinetic_total + self.potential_total
    }

    /// Update potential total from components
    pub fn update_potential_total(&mut self) {
        self.potential_total =
            self.bond_total +
            self.angle_total +
            self.dihedral_total +
            self.improper_total +
            self.cross_dihedral_total +
            self.lj_total +
            self.crf_total +
            self.ls_total +
            self.special_total +
            self.sasa_total;
    }

    /// Clear all energies
    pub fn clear(&mut self) {
        self.kinetic_total = 0.0;
        self.potential_total = 0.0;
        self.bond_total = 0.0;
        self.angle_total = 0.0;
        self.dihedral_total = 0.0;
        self.improper_total = 0.0;
        self.cross_dihedral_total = 0.0;
        self.lj_total = 0.0;
        self.crf_total = 0.0;
        self.ls_total = 0.0;
        self.ls_realspace_total = 0.0;
        self.ls_kspace_total = 0.0;
        self.special_total = 0.0;
        self.sasa_total = 0.0;
        self.constraint_total = 0.0;
        self.virial_total = 0.0;

        for ke in &mut self.kinetic_energy {
            *ke = 0.0;
        }

        for row in &mut self.lj_energy {
            for e in row {
                *e = 0.0;
            }
        }

        for row in &mut self.crf_energy {
            for e in row {
                *e = 0.0;
            }
        }
    }
}

/// State of the system at a single time point
///
/// Contains positions, velocities, forces, and all associated properties
#[derive(Debug, Clone)]
pub struct State {
    pub pos: Vec<Vec3>,              // Positions
    pub vel: Vec<Vec3>,              // Velocities
    pub force: Vec<Vec3>,            // Forces
    pub constraint_force: Vec<Vec3>, // Forces from constraints
    pub box_config: Box,             // Simulation box
    pub virial_tensor: Mat3,         // Virial tensor (for pressure calculation)
    pub energies: Energy,            // Energy storage
}

impl State {
    pub fn new(n_atoms: usize, num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            pos: vec![Vec3::ZERO; n_atoms],
            vel: vec![Vec3::ZERO; n_atoms],
            force: vec![Vec3::ZERO; n_atoms],
            constraint_force: vec![Vec3::ZERO; n_atoms],
            box_config: Box::vacuum(),
            virial_tensor: Mat3::ZERO,
            energies: Energy::new(num_temperature_groups, num_energy_groups),
        }
    }

    /// Clear forces and virial
    pub fn clear_forces(&mut self) {
        for f in &mut self.force {
            *f = Vec3::ZERO;
        }
        self.virial_tensor = Mat3::ZERO;
    }

    /// Clear constraint forces
    pub fn clear_constraint_forces(&mut self) {
        for f in &mut self.constraint_force {
            *f = Vec3::ZERO;
        }
    }

    /// Calculate kinetic energy
    ///
    /// KE = 0.5 * sum_i(m_i * v_i^2)
    pub fn calculate_kinetic_energy(&mut self, masses: &[f64]) {
        self.energies.kinetic_total = self.vel.iter()
            .zip(masses.iter())
            .map(|(v, &m)| 0.5 * m * (v.length_squared() as f64))
            .sum();
    }

    /// Calculate kinetic energy per group
    pub fn calculate_kinetic_energy_groups(
        &mut self,
        masses: &[f64],
        temperature_groups: &[Vec<usize>],
    ) {
        for (group_idx, group) in temperature_groups.iter().enumerate() {
            self.energies.kinetic_energy[group_idx] = group.iter()
                .map(|&i| 0.5 * masses[i] * (self.vel[i].length_squared() as f64))
                .sum();
        }

        // Total kinetic energy is sum of all groups
        self.energies.kinetic_total = self.energies.kinetic_energy.iter().sum();
    }

    /// Calculate temperature from kinetic energy
    ///
    /// T = 2 * KE / (k_B * N_dof)
    /// where N_dof = 3N - N_constraints
    pub fn temperature(&self, num_degrees_of_freedom: usize) -> f64 {
        const BOLTZMANN: f64 = 0.00831446; // kJ/(mol·K)
        2.0 * self.energies.kinetic_total / (BOLTZMANN * num_degrees_of_freedom as f64)
    }

    /// Calculate pressure from virial
    ///
    /// P = (2*KE - Virial) / (3*V)
    pub fn pressure(&self) -> f64 {
        let volume = self.box_config.volume();
        let virial_trace = self.virial_tensor.x_axis.x as f64 +
                          self.virial_tensor.y_axis.y as f64 +
                          self.virial_tensor.z_axis.z as f64;

        (2.0 * self.energies.kinetic_total - virial_trace) / (3.0 * volume)
    }
}

/// Configuration - manages current and previous states with double buffering
///
/// GROMOS uses double buffering (pointer swapping) for efficient state management.
/// This is translated to Rust using two State structs and an index.
#[derive(Debug, Clone)]
pub struct Configuration {
    state1: State,
    state2: State,
    current_idx: usize,  // 0 or 1

    /// Lattice shifts for periodic boundary tracking (FEP)
    /// Stores how many times each atom crossed periodic boundaries
    /// Dimension: [n_atoms][3] for x, y, z shifts
    pub lattice_shifts: Vec<[i32; 3]>,
}

impl Configuration {
    pub fn new(n_atoms: usize, num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            state1: State::new(n_atoms, num_temperature_groups, num_energy_groups),
            state2: State::new(n_atoms, num_temperature_groups, num_energy_groups),
            current_idx: 0,
            lattice_shifts: vec![[0, 0, 0]; n_atoms],
        }
    }

    /// Get reference to current state
    #[inline]
    pub fn current(&self) -> &State {
        if self.current_idx == 0 {
            &self.state1
        } else {
            &self.state2
        }
    }

    /// Get mutable reference to current state
    #[inline]
    pub fn current_mut(&mut self) -> &mut State {
        if self.current_idx == 0 {
            &mut self.state1
        } else {
            &mut self.state2
        }
    }

    /// Get reference to old (previous) state
    #[inline]
    pub fn old(&self) -> &State {
        if self.current_idx == 0 {
            &self.state2
        } else {
            &self.state1
        }
    }

    /// Get mutable reference to old state
    #[inline]
    pub fn old_mut(&mut self) -> &mut State {
        if self.current_idx == 0 {
            &mut self.state2
        } else {
            &mut self.state1
        }
    }

    /// Exchange current and old states (zero-cost pointer swap)
    ///
    /// This is the Rust equivalent of GROMOS C++ pointer swapping
    #[inline]
    pub fn exchange_state(&mut self) {
        self.current_idx = 1 - self.current_idx;
    }

    /// Copy current state to old state
    pub fn copy_current_to_old(&mut self) {
        if self.current_idx == 0 {
            self.state2 = self.state1.clone();
        } else {
            self.state1 = self.state2.clone();
        }
    }
}

/// Stochastic dynamics additional variables
#[derive(Debug, Clone)]
pub struct StochasticVariables {
    pub gamma: Vec<f64>,      // Friction coefficients per atom
    pub random_force: Vec<Vec3>,  // Random forces for Langevin dynamics
}

impl StochasticVariables {
    pub fn new(n_atoms: usize) -> Self {
        Self {
            gamma: vec![0.0; n_atoms],
            random_force: vec![Vec3::ZERO; n_atoms],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_box_rectangular() {
        let box_config = Box::rectangular(3.0, 4.0, 5.0);

        assert_eq!(box_config.box_type, BoxType::Rectangular);
        assert_eq!(box_config.dimensions(), Vec3::new(3.0, 4.0, 5.0));
        assert!((box_config.volume() - 60.0).abs() < 1e-6);
    }

    #[test]
    fn test_energy_total() {
        let mut energy = Energy::new(1, 1);

        energy.kinetic_total = 100.0;
        energy.bond_total = 50.0;
        energy.lj_total = 30.0;
        energy.crf_total = 20.0;

        energy.update_potential_total();

        assert_eq!(energy.potential_total, 100.0);
        assert_eq!(energy.total(), 200.0);
    }

    #[test]
    fn test_state_kinetic_energy() {
        let mut state = State::new(2, 1, 1);
        let masses = vec![1.0, 2.0];

        // Set velocities
        state.vel[0] = Vec3::new(1.0, 0.0, 0.0);
        state.vel[1] = Vec3::new(0.0, 2.0, 0.0);

        state.calculate_kinetic_energy(&masses);

        // KE = 0.5 * (1.0 * 1.0^2 + 2.0 * 2.0^2) = 0.5 * (1 + 8) = 4.5
        assert!((state.energies.kinetic_total - 4.5).abs() < 1e-6);
    }

    #[test]
    fn test_temperature_calculation() {
        let mut state = State::new(100, 1, 1);
        let masses = vec![1.0; 100];

        // Set random velocities for testing
        for vel in &mut state.vel {
            *vel = Vec3::new(0.1, 0.1, 0.1);
        }

        state.calculate_kinetic_energy(&masses);

        // 3N degrees of freedom (no constraints)
        let temp = state.temperature(300);

        assert!(temp > 0.0);
        assert!(temp < 100.0); // Sanity check
    }

    #[test]
    fn test_state_exchange() {
        let mut conf = Configuration::new(2, 1, 1);

        // Set some values in current state
        conf.current_mut().pos[0] = Vec3::new(1.0, 2.0, 3.0);

        // Exchange
        conf.exchange_state();

        // Old state should now have the values
        assert_eq!(conf.old().pos[0], Vec3::new(1.0, 2.0, 3.0));

        // Current state should be fresh
        assert_eq!(conf.current().pos[0], Vec3::ZERO);
    }

    #[test]
    fn test_clear_forces() {
        let mut state = State::new(10, 1, 1);

        // Set some forces
        state.force[0] = Vec3::new(1.0, 2.0, 3.0);
        state.virial_tensor = Mat3::from_cols(
            Vec3::new(10.0, 0.0, 0.0),
            Vec3::new(0.0, 10.0, 0.0),
            Vec3::new(0.0, 0.0, 10.0),
        );

        state.clear_forces();

        assert_eq!(state.force[0], Vec3::ZERO);
        assert_eq!(state.virial_tensor, Mat3::ZERO);
    }
}
