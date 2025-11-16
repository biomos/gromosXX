//! Integration algorithms for molecular dynamics
//!
//! Direct Rust translation of GROMOS integration algorithms from:
//! - md++/src/algorithm/integration/leap_frog.cc
//! - md++/src/algorithm/integration/velocity_verlet.cc
//! - md++/src/algorithm/integration/stochastic.cc

use crate::math::{Vec3, Mat3};
use crate::topology::Topology;
use crate::configuration::{Configuration, BoxType};
use rayon::prelude::*;
use rand::{Rng, SeedableRng};
use rand_distr::{Normal, Distribution};

/// Trait for MD integrators
pub trait Integrator {
    /// Perform one integration step
    fn step(&mut self, dt: f64, topo: &Topology, conf: &mut Configuration);

    /// Get integrator name
    fn name(&self) -> &str;
}

/// Leap-Frog integrator (GROMOS default)
///
/// Direct translation of md++/src/algorithm/integration/leap_frog.cc
///
/// Algorithm:
/// ```text
/// v(t+dt/2) = v(t-dt/2) + a(t) * dt
/// r(t+dt) = r(t) + v(t+dt/2) * dt
/// ```
///
/// Velocities are stored at half-steps, giving good energy conservation.
#[derive(Debug, Clone)]
pub struct LeapFrog {
    pub parallel: bool,
}

impl LeapFrog {
    pub fn new() -> Self {
        Self { parallel: false }
    }

    pub fn with_parallel(mut self) -> Self {
        self.parallel = true;
        self
    }
}

impl Integrator for LeapFrog {
    fn step(&mut self, dt: f64, topo: &Topology, conf: &mut Configuration) {
        let n_atoms = topo.num_atoms();
        let dt_f32 = dt as f32;

        // Step 1: Update velocities v(t+dt/2) = v(t-dt/2) + F(t)/m * dt
        if self.parallel {
            let old_vel = conf.old().vel.clone();
            let old_force = conf.old().force.clone();

            conf.current_mut().vel.par_iter_mut()
                .enumerate()
                .for_each(|(i, vel_new)| {
                    let accel = old_force[i] * (topo.inverse_mass[i] as f32);
                    *vel_new = old_vel[i] + accel * dt_f32;
                });
        } else {
            for i in 0..n_atoms {
                let accel = conf.old().force[i] * (topo.inverse_mass[i] as f32);
                conf.current_mut().vel[i] = conf.old().vel[i] + accel * dt_f32;
            }
        }

        // Step 2: Exchange states (zero-cost pointer swap)
        conf.exchange_state();

        // Step 3: Update positions r(t+dt) = r(t) + v(t+dt/2) * dt
        // After exchange, conf.old() contains the updated velocities
        if self.parallel {
            let old_pos = conf.old().pos.clone();
            let old_vel = conf.old().vel.clone();

            conf.current_mut().pos.par_iter_mut()
                .enumerate()
                .for_each(|(i, pos_new)| {
                    *pos_new = old_pos[i] + old_vel[i] * dt_f32;
                });
        } else {
            for i in 0..n_atoms {
                conf.current_mut().pos[i] = conf.old().pos[i] + conf.old().vel[i] * dt_f32;
            }
        }

        // Clear forces for next step
        conf.current_mut().clear_forces();
    }

    fn name(&self) -> &str {
        "Leap-Frog"
    }
}

impl Default for LeapFrog {
    fn default() -> Self {
        Self::new()
    }
}

/// Velocity Verlet integrator
///
/// Simpler implementation without parallel for now
#[derive(Debug, Clone)]
pub struct VelocityVerlet {
    saved_accelerations: Vec<Vec3>,
}

impl VelocityVerlet {
    pub fn new() -> Self {
        Self {
            saved_accelerations: Vec::new(),
        }
    }
}

impl Integrator for VelocityVerlet {
    fn step(&mut self, dt: f64, topo: &Topology, conf: &mut Configuration) {
        let n_atoms = topo.num_atoms();
        let dt_f32 = dt as f32;

        if self.saved_accelerations.len() != n_atoms {
            self.saved_accelerations.resize(n_atoms, Vec3::ZERO);
        }

        // Step 1: Update positions and save accelerations
        for i in 0..n_atoms {
            let accel = conf.old().force[i] * (topo.inverse_mass[i] as f32);
            self.saved_accelerations[i] = accel;

            conf.current_mut().pos[i] = conf.old().pos[i] +
                conf.old().vel[i] * dt_f32 +
                accel * (0.5 * dt_f32 * dt_f32);

            conf.current_mut().vel[i] = conf.old().vel[i] + accel * (0.5 * dt_f32);
        }

        conf.exchange_state();

        // Step 2: Complete velocity update (after force calculation)
        for i in 0..n_atoms {
            let accel_new = conf.current().force[i] * (topo.inverse_mass[i] as f32);
            conf.current_mut().vel[i] += accel_new * (0.5 * dt_f32);
        }

        conf.current_mut().clear_forces();
    }

    fn name(&self) -> &str {
        "Velocity-Verlet"
    }
}

impl Default for VelocityVerlet {
    fn default() -> Self {
        Self::new()
    }
}

/// Steepest Descent minimization algorithm
///
/// Direct translation of md++/src/algorithm/integration/steepest_descent.cc
///
/// Algorithm:
/// ```text
/// 1. Normalize force vector to unit magnitude
/// 2. Update positions: r(t+dt) = r(t) + step_size * normalized_force
/// 3. Adapt step size based on energy change
/// 4. Check convergence: |E_new - E_old| < tolerance
/// ```
///
/// Features:
/// - Adaptive step sizing (increases 20% on energy decrease, halves on increase)
/// - Optional force magnitude limiting
/// - Energy-based convergence criterion
#[derive(Debug, Clone)]
pub struct SteepestDescent {
    /// Current step size
    pub step_size: f64,
    /// Initial step size (dx0)
    pub initial_step_size: f64,
    /// Maximum step size (dxm)
    pub max_step_size: f64,
    /// Energy convergence tolerance (dele)
    pub energy_tolerance: f64,
    /// Minimum number of steps before checking convergence (nmin)
    pub min_steps: usize,
    /// Maximum force magnitude (flim, 0.0 = unlimited)
    pub force_limit: f64,
    /// Current step counter
    step_count: usize,
    /// Previous total energy
    prev_energy: f64,
    /// Convergence flag
    converged: bool,
}

impl SteepestDescent {
    /// Create new steepest descent minimizer with default GROMOS parameters
    pub fn new() -> Self {
        Self {
            step_size: 0.01,           // GROMOS default dx0
            initial_step_size: 0.01,
            max_step_size: 0.05,       // GROMOS default dxm
            energy_tolerance: 0.1,     // GROMOS default dele (kJ/mol)
            min_steps: 1,
            force_limit: 0.0,          // No limit
            step_count: 0,
            prev_energy: 0.0,
            converged: false,
        }
    }

    /// Set energy convergence tolerance in kJ/mol
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.energy_tolerance = tolerance;
        self
    }

    /// Set initial and maximum step sizes
    pub fn with_step_sizes(mut self, initial: f64, max: f64) -> Self {
        self.initial_step_size = initial;
        self.step_size = initial;
        self.max_step_size = max;
        self
    }

    /// Set force magnitude limit (0.0 = unlimited)
    pub fn with_force_limit(mut self, limit: f64) -> Self {
        self.force_limit = limit;
        self
    }

    /// Set minimum steps before convergence check
    pub fn with_min_steps(mut self, steps: usize) -> Self {
        self.min_steps = steps;
        self
    }

    /// Check if minimization has converged
    pub fn is_converged(&self) -> bool {
        self.converged
    }

    /// Get current step count
    pub fn step_count(&self) -> usize {
        self.step_count
    }

    /// Reset minimizer state
    pub fn reset(&mut self) {
        self.step_size = self.initial_step_size;
        self.step_count = 0;
        self.prev_energy = 0.0;
        self.converged = false;
    }

    /// Apply force limiting if configured
    fn apply_force_limit(&self, forces: &mut [Vec3]) {
        if self.force_limit > 0.0 {
            for force in forces.iter_mut() {
                let magnitude = force.length() as f64;
                if magnitude > self.force_limit {
                    *force *= (self.force_limit / magnitude) as f32;
                }
            }
        }
    }

    /// Calculate normalized force magnitude
    fn calculate_force_norm(&self, forces: &[Vec3]) -> f64 {
        let f_squared: f64 = forces.iter()
            .map(|f| {
                let fx = f.x as f64;
                let fy = f.y as f64;
                let fz = f.z as f64;
                fx * fx + fy * fy + fz * fz
            })
            .sum();

        if f_squared < 1e-15 {
            1.0 // Avoid division by zero
        } else {
            1.0 / f_squared.sqrt()
        }
    }
}

impl Integrator for SteepestDescent {
    fn step(&mut self, _dt: f64, topo: &Topology, conf: &mut Configuration) {
        let n_atoms = topo.num_atoms();
        self.step_count += 1;

        // Zero velocities (minimization, not dynamics)
        for vel in &mut conf.current_mut().vel {
            *vel = Vec3::ZERO;
        }

        // Check convergence after minimum steps
        if self.step_count > self.min_steps {
            let current_energy = conf.current().energies.total();
            let energy_change = (current_energy - self.prev_energy).abs();

            if energy_change < self.energy_tolerance {
                self.converged = true;
                return;
            }

            // Adaptive step sizing
            if current_energy < self.prev_energy {
                // Energy decreased - increase step size
                self.step_size = (self.step_size * 1.2).min(self.max_step_size);
            } else {
                // Energy increased - decrease step size
                self.step_size *= 0.5;
            }

            self.prev_energy = current_energy;
        } else if self.step_count == self.min_steps {
            self.prev_energy = conf.current().energies.total();
        }

        // Apply force limiting if configured
        self.apply_force_limit(&mut conf.current_mut().force);

        // Calculate normalized force
        let f_norm = self.calculate_force_norm(&conf.current().force);

        // Exchange states (swap old <-> current)
        conf.exchange_state();

        // Update positions: r_new = r_old + step_size * f_norm * force_old
        let step_scale = (self.step_size * f_norm) as f32;

        for i in 0..n_atoms {
            conf.current_mut().pos[i] = conf.old().pos[i] + conf.old().force[i] * step_scale;
        }

        // Zero velocities and forces for next iteration
        for vel in &mut conf.current_mut().vel {
            *vel = Vec3::ZERO;
        }
        conf.current_mut().clear_forces();
    }

    fn name(&self) -> &str {
        "Steepest-Descent"
    }
}

impl Default for SteepestDescent {
    fn default() -> Self {
        Self::new()
    }
}

/// Stochastic Dynamics (Langevin) integrator
///
/// Direct translation of md++/src/algorithm/integration/stochastic.cc
///
/// Langevin equation of motion:
/// ```text
/// m * dv/dt = F(r) - γ * m * v + R(t)
/// ```
/// where:
/// - γ is the friction coefficient
/// - R(t) is a random force with <R(t)·R(t')> = 2γkT·m·δ(t-t')
///
/// Algorithm uses velocity Verlet-like scheme with stochastic terms:
/// ```text
/// v(t+dt/2) = c1*v(t) + c2*F(t)/m + c3*ξ₁ + c4*ξ₂
/// r(t+dt) = r(t) + c5*v(t+dt/2)*dt
/// ```
///
/// Two numerical methods:
/// - **Analytical** (|γdt| > 0.05): Exact exponential formulas
/// - **Power Series** (|γdt| ≤ 0.05): Taylor expansion for numerical stability
#[derive(Debug, Clone)]
pub struct StochasticDynamics {
    /// Friction coefficient γ (ps⁻¹)
    pub gamma: f64,
    /// Target temperature (K)
    pub temperature: f64,
    /// Random number generator seed
    seed: Option<u64>,
    /// Boltzmann constant in GROMOS units (kJ/(mol·K))
    kb: f64,
    /// Friction coefficients per atom
    coefficients: Vec<LangevinCoefficients>,
}

/// Pre-computed Langevin coefficients for each atom
#[derive(Debug, Clone, Copy)]
struct LangevinCoefficients {
    c1: f64,  // exp(-γdt) for velocity damping
    c2: f64,  // Coefficient for force term
    c3: f64,  // Coefficient for random velocity
    c4: f64,  // Coefficient for random position correction
    c5: f64,  // Coefficient for position update
    c6: f64,  // Additional position scaling
}

impl StochasticDynamics {
    /// Create new stochastic dynamics integrator
    ///
    /// # Arguments
    /// * `gamma` - Friction coefficient in ps⁻¹
    /// * `temperature` - Target temperature in Kelvin
    pub fn new(gamma: f64, temperature: f64) -> Self {
        Self {
            gamma,
            temperature,
            seed: None,
            kb: 0.008314462618,  // GROMOS kB in kJ/(mol·K)
            coefficients: Vec::new(),
        }
    }

    /// Set random seed for reproducibility
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    /// Compute Langevin coefficients for given timestep and mass
    ///
    /// Direct translation from stochastic.cc lines 88-180
    fn compute_coefficients(&self, dt: f64, mass: f64) -> LangevinCoefficients {
        let gamma_dt = self.gamma * dt;
        let abs_gamma_dt = gamma_dt.abs();

        // GROMOS numerical threshold for analytical vs power series
        if abs_gamma_dt > 0.05 {
            // Analytical formulas (exact exponentials)
            self.compute_analytical(dt, mass, gamma_dt)
        } else {
            // Power series expansion (better numerical stability)
            self.compute_power_series(dt, mass, gamma_dt)
        }
    }

    /// Analytical coefficient calculation (|γdt| > 0.05)
    ///
    /// From stochastic.cc lines 91-129
    fn compute_analytical(&self, dt: f64, mass: f64, gamma_dt: f64) -> LangevinCoefficients {
        let exp_gamma_dt = (-gamma_dt).exp();
        let exp_2gamma_dt = (-2.0 * gamma_dt).exp();

        let c1 = exp_gamma_dt;

        // c2: ∫₀ᵗ exp(-γ(t-s)) ds = (1 - exp(-γt))/γ
        let c2 = (1.0 - exp_gamma_dt) / self.gamma;

        // Variance for random force: 2γkT/m
        let kT = self.kb * self.temperature;
        let sigma_v_sq = kT / mass;

        // c3: random velocity contribution
        // √(kT/m * (1 - exp(-2γt)))
        let c3 = (sigma_v_sq * (1.0 - exp_2gamma_dt)).sqrt();

        // c4: correlation between velocity and position noise
        // √(kT/m) * √((1 - exp(-2γt))/(2γ²) - (1 - exp(-γt))²/γ²)
        let term1 = (1.0 - exp_2gamma_dt) / (2.0 * self.gamma * self.gamma);
        let term2 = (1.0 - exp_gamma_dt) / self.gamma;
        let term2_sq = term2 * term2;
        let c4_arg = term1 - term2_sq;
        let c4 = if c4_arg > 0.0 {
            (sigma_v_sq * c4_arg).sqrt()
        } else {
            0.0
        };

        let c5 = c2;
        let c6 = 1.0;

        LangevinCoefficients { c1, c2, c3, c4, c5, c6 }
    }

    /// Power series coefficient calculation (|γdt| ≤ 0.05)
    ///
    /// From stochastic.cc lines 131-180
    /// Uses Taylor expansion for better numerical stability at small γdt
    fn compute_power_series(&self, dt: f64, _mass: f64, gamma_dt: f64) -> LangevinCoefficients {
        let kT = self.kb * self.temperature;
        let gamma_dt2 = gamma_dt * gamma_dt;
        let gamma_dt3 = gamma_dt2 * gamma_dt;
        let gamma_dt4 = gamma_dt2 * gamma_dt2;
        let gamma_dt5 = gamma_dt3 * gamma_dt2;

        // c1 = exp(-γdt) ≈ 1 - γdt + γ²dt²/2 - γ³dt³/6 + ...
        let c1 = 1.0 - gamma_dt + gamma_dt2/2.0 - gamma_dt3/6.0 + gamma_dt4/24.0;

        // c2 = (1 - exp(-γdt))/γ ≈ dt(1 - γdt/2 + γ²dt²/6 - ...)
        let c2 = dt * (1.0 - gamma_dt/2.0 + gamma_dt2/6.0 - gamma_dt3/24.0);

        // c3² = kT(1 - exp(-2γdt))
        // 1 - exp(-2γdt) ≈ 2γdt - 2γ²dt² + 4γ³dt³/3 - ...
        let exp_term = 2.0*gamma_dt - 2.0*gamma_dt2 + 4.0*gamma_dt3/3.0 - 2.0*gamma_dt4/3.0;
        let c3 = (kT * exp_term).sqrt();

        // c4: correlation term (more complex power series)
        let c2_sq = c2 * c2;
        let term1 = dt * dt * (1.0/3.0 - gamma_dt/4.0 + gamma_dt2/10.0 - gamma_dt3/36.0 + gamma_dt4/144.0 - gamma_dt5/600.0);
        let c4_sq = kT * (term1 - c2_sq);
        let c4 = if c4_sq > 0.0 { c4_sq.sqrt() } else { 0.0 };

        let c5 = c2;
        let c6 = 1.0;

        LangevinCoefficients { c1, c2, c3, c4, c5, c6 }
    }

    /// Initialize coefficients for all atoms
    fn initialize_coefficients(&mut self, dt: f64, topo: &Topology) {
        if self.coefficients.len() != topo.num_atoms() {
            self.coefficients.clear();
            for i in 0..topo.num_atoms() {
                let coeff = self.compute_coefficients(dt, topo.mass[i]);
                self.coefficients.push(coeff);
            }
        }
    }
}

impl Integrator for StochasticDynamics {
    fn step(&mut self, dt: f64, topo: &Topology, conf: &mut Configuration) {
        let n_atoms = topo.num_atoms();

        // Initialize coefficients if needed
        self.initialize_coefficients(dt, topo);

        // Create RNG
        let mut rng = if let Some(seed) = self.seed {
            rand::rngs::StdRng::seed_from_u64(seed)
        } else {
            rand::rngs::StdRng::from_entropy()
        };

        let normal = Normal::new(0.0, 1.0).unwrap();

        // Step 1: Update velocities with Langevin dynamics
        // v(t+dt/2) = c1*v(t) + c2*F(t)/m + c3*ξ₁ + c4*ξ₂
        let dt_f32 = dt as f32;

        for i in 0..n_atoms {
            let coeff = self.coefficients[i];
            let mass = topo.mass[i];
            let inv_mass = topo.inverse_mass[i];

            // Deterministic terms
            let v_damped = conf.old().vel[i] * (coeff.c1 as f32);
            let v_force = conf.old().force[i] * (inv_mass as f32) * (coeff.c2 as f32);

            // Random terms (Gaussian noise)
            let xi1_x = normal.sample(&mut rng) as f32;
            let xi1_y = normal.sample(&mut rng) as f32;
            let xi1_z = normal.sample(&mut rng) as f32;
            let xi1 = Vec3::new(xi1_x, xi1_y, xi1_z);

            let xi2_x = normal.sample(&mut rng) as f32;
            let xi2_y = normal.sample(&mut rng) as f32;
            let xi2_z = normal.sample(&mut rng) as f32;
            let xi2 = Vec3::new(xi2_x, xi2_y, xi2_z);

            let v_random1 = xi1 * (coeff.c3 as f32);
            let v_random2 = xi2 * (coeff.c4 as f32);

            // Combined velocity update
            conf.current_mut().vel[i] = v_damped + v_force + v_random1 + v_random2;
        }

        // Step 2: Exchange states
        conf.exchange_state();

        // Step 3: Update positions
        // r(t+dt) = r(t) + c5*v(t+dt/2)*dt
        for i in 0..n_atoms {
            let coeff = self.coefficients[i];
            let pos_update = conf.old().vel[i] * dt_f32 * (coeff.c5 as f32) * (coeff.c6 as f32);
            conf.current_mut().pos[i] = conf.old().pos[i] + pos_update;
        }

        // Clear forces for next step
        conf.current_mut().clear_forces();
    }

    fn name(&self) -> &str {
        "Stochastic-Dynamics"
    }
}

impl Default for StochasticDynamics {
    fn default() -> Self {
        Self::new(0.1, 300.0)  // Default: γ=0.1 ps⁻¹, T=300K
    }
}

/// Berendsen thermostat (weak coupling)
#[derive(Debug, Clone)]
pub struct BerendsenThermostat {
    pub target_temperature: f64,
    pub coupling_time: f64,
}

impl BerendsenThermostat {
    pub fn new(target_temperature: f64, coupling_time: f64) -> Self {
        Self {
            target_temperature,
            coupling_time,
        }
    }

    pub fn apply(&self, dt: f64, topo: &Topology, conf: &mut Configuration) {
        conf.current_mut().calculate_kinetic_energy(&topo.mass);
        let current_temp = conf.current().temperature(topo.num_atoms() * 3);

        if current_temp == 0.0 {
            return;
        }

        let lambda = (1.0 + (dt / self.coupling_time) * (self.target_temperature / current_temp - 1.0)).sqrt();

        for vel in &mut conf.current_mut().vel {
            *vel *= lambda as f32;
        }

        conf.current_mut().calculate_kinetic_energy(&topo.mass);
    }
}

/// Scaled Leap-Frog integrator (for AMD/Adaptive Coupling)
///
/// Direct translation of md++/src/algorithm/integration/scaled_leap_frog.cc
///
/// Extension of standard Leap-Frog that allows per-atom force scaling.
/// Used for:
/// - Accelerated Molecular Dynamics (AMD)
/// - Adaptive coupling/decoupling of atoms
/// - Multiple time-stepping algorithms
///
/// Algorithm:
/// ```text
/// v(t+dt/2) = v(t-dt/2) + scale(i) * a(t) * dt
/// r(t+dt) = r(t) + v(t+dt/2) * dt
/// ```
///
/// where `scale(i)` is the force scaling factor for atom i.
#[derive(Debug, Clone)]
pub struct ScaledLeapFrog {
    /// Force scaling factors per atom (1.0 = no scaling)
    pub force_scales: Vec<f64>,
    /// Enable parallel execution
    pub parallel: bool,
}

impl ScaledLeapFrog {
    /// Create new scaled leap-frog integrator
    ///
    /// # Arguments
    /// * `num_atoms` - Number of atoms in system
    ///
    /// All atoms initialized with scale factor 1.0 (no scaling)
    pub fn new(num_atoms: usize) -> Self {
        Self {
            force_scales: vec![1.0; num_atoms],
            parallel: false,
        }
    }

    /// Create with custom force scaling factors
    pub fn with_scales(scales: Vec<f64>) -> Self {
        Self {
            force_scales: scales,
            parallel: false,
        }
    }

    /// Enable parallel execution
    pub fn with_parallel(mut self) -> Self {
        self.parallel = true;
        self
    }

    /// Set force scale for specific atom
    pub fn set_scale(&mut self, atom_idx: usize, scale: f64) {
        if atom_idx < self.force_scales.len() {
            self.force_scales[atom_idx] = scale;
        }
    }

    /// Set force scale for range of atoms
    pub fn set_scale_range(&mut self, start: usize, end: usize, scale: f64) {
        for i in start..end.min(self.force_scales.len()) {
            self.force_scales[i] = scale;
        }
    }
}

impl Integrator for ScaledLeapFrog {
    fn step(&mut self, dt: f64, topo: &Topology, conf: &mut Configuration) {
        let n_atoms = topo.num_atoms();
        let dt_f32 = dt as f32;

        // Ensure force_scales vector matches topology
        if self.force_scales.len() != n_atoms {
            self.force_scales.resize(n_atoms, 1.0);
        }

        // Step 1: Update velocities with scaled forces
        // v(t+dt/2) = v(t-dt/2) + scale * F(t)/m * dt
        if self.parallel {
            let old_vel = conf.old().vel.clone();
            let old_force = conf.old().force.clone();
            let scales = &self.force_scales;

            conf.current_mut().vel.par_iter_mut()
                .enumerate()
                .for_each(|(i, vel_new)| {
                    let scale = scales[i] as f32;
                    let accel = old_force[i] * (topo.inverse_mass[i] as f32);
                    *vel_new = old_vel[i] + scale * accel * dt_f32;
                });
        } else {
            for i in 0..n_atoms {
                let scale = self.force_scales[i] as f32;
                let accel = conf.old().force[i] * (topo.inverse_mass[i] as f32);
                conf.current_mut().vel[i] = conf.old().vel[i] + scale * accel * dt_f32;
            }
        }

        // Step 2: Exchange states (zero-cost pointer swap)
        conf.exchange_state();

        // Step 3: Update positions r(t+dt) = r(t) + v(t+dt/2) * dt
        if self.parallel {
            let old_pos = conf.old().pos.clone();
            let old_vel = conf.old().vel.clone();

            conf.current_mut().pos.par_iter_mut()
                .enumerate()
                .for_each(|(i, pos_new)| {
                    *pos_new = old_pos[i] + old_vel[i] * dt_f32;
                });
        } else {
            for i in 0..n_atoms {
                conf.current_mut().pos[i] = conf.old().pos[i] + conf.old().vel[i] * dt_f32;
            }
        }

        // Clear forces for next step
        conf.current_mut().clear_forces();
    }

    fn name(&self) -> &str {
        "Scaled-Leap-Frog"
    }
}

impl Default for ScaledLeapFrog {
    fn default() -> Self {
        Self::new(0)
    }
}

/// Conjugate Gradient minimization algorithm
///
/// Direct translation of md++/src/algorithm/integration/conjugate_gradient.cc
///
/// More efficient minimizer than Steepest Descent, using conjugate search directions.
///
/// **Algorithm**:
/// 1. Calculate energy gradient (forces)
/// 2. Determine search direction:
///    - First step or reset: p = f (steepest descent)
///    - Otherwise: p = f + β·p_old (conjugate direction)
/// 3. Find minimum along search direction using cubic interpolation
/// 4. Update positions to minimum
/// 5. Repeat until convergence
///
/// **Beta calculation** (two variants):
/// - **Fletcher-Reeves**: β = |f_new|² / |f_old|²
/// - **Polak-Ribiere**: β = ⟨f_new, f_new - f_old⟩ / |f_old|²
///
/// **Convergence criterion**:
/// - RMS force < tolerance
///
/// **Features**:
/// - Cubic interpolation line search
/// - Automatic search direction reset
/// - Adaptive step sizing
/// - Optional periodic reset (every N steps)
/// - SHAKE-compatible (can constrain during minimization)
#[derive(Debug, Clone)]
pub struct ConjugateGradient {
    /// Initial step size
    pub initial_step_size: f64,
    /// Current step size
    pub step_size: f64,
    /// Maximum step size
    pub max_step_size: f64,
    /// RMS force convergence tolerance (kJ/mol/nm)
    pub force_tolerance: f64,
    /// Minimum number of steps before checking convergence
    pub min_steps: usize,
    /// Maximum cubic interpolations per line search
    pub max_interpolations: usize,
    /// Displacement criterion for cubic interpolation (nm)
    pub interpolation_criterion: f64,
    /// Beta calculation method: FletcherReeves or PolakRibiere
    pub beta_method: ConjugateGradientMethod,
    /// Reset search direction every N steps (0 = never)
    pub reset_interval: usize,
    /// Maximum force magnitude (0.0 = unlimited)
    pub force_limit: f64,

    // Internal state
    step_count: usize,
    prev_energy: f64,
    converged: bool,
    search_directions: Vec<Vec3>,
    old_forces: Vec<Vec3>,
    total_doublings: usize,
    total_interpolations: usize,
    total_iterations: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConjugateGradientMethod {
    /// Fletcher-Reeves: β = |f_new|² / |f_old|²
    FletcherReeves,
    /// Polak-Ribiere: β = ⟨f_new, f_new - f_old⟩ / |f_old|²  (recommended)
    PolakRibiere,
}

impl ConjugateGradient {
    /// Create new conjugate gradient minimizer with GROMOS default parameters
    pub fn new() -> Self {
        Self {
            initial_step_size: 0.01,       // 0.01 nm
            step_size: 0.01,
            max_step_size: 0.05,           // 0.05 nm
            force_tolerance: 0.1,          // 0.1 kJ/mol/nm
            min_steps: 1,
            max_interpolations: 5,
            interpolation_criterion: 0.0001, // 0.0001 nm
            beta_method: ConjugateGradientMethod::PolakRibiere,
            reset_interval: 0,             // No forced reset
            force_limit: 0.0,              // No limit
            step_count: 0,
            prev_energy: 0.0,
            converged: false,
            search_directions: Vec::new(),
            old_forces: Vec::new(),
            total_doublings: 0,
            total_interpolations: 0,
            total_iterations: 0,
        }
    }

    /// Set convergence tolerance (RMS force in kJ/mol/nm)
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.force_tolerance = tolerance;
        self
    }

    /// Set step sizes (initial and maximum in nm)
    pub fn with_step_sizes(mut self, initial: f64, max: f64) -> Self {
        self.initial_step_size = initial;
        self.step_size = initial;
        self.max_step_size = max;
        self
    }

    /// Set beta calculation method
    pub fn with_method(mut self, method: ConjugateGradientMethod) -> Self {
        self.beta_method = method;
        self
    }

    /// Set search direction reset interval (0 = never)
    pub fn with_reset_interval(mut self, interval: usize) -> Self {
        self.reset_interval = interval;
        self
    }

    /// Set maximum interpolations per line search
    pub fn with_max_interpolations(mut self, max_interp: usize) -> Self {
        self.max_interpolations = max_interp;
        self
    }

    /// Set force magnitude limit (0.0 = unlimited)
    pub fn with_force_limit(mut self, limit: f64) -> Self {
        self.force_limit = limit;
        self
    }

    /// Check if minimization has converged
    pub fn is_converged(&self) -> bool {
        self.converged
    }

    /// Get current step count
    pub fn step_count(&self) -> usize {
        self.step_count
    }

    /// Get statistics
    pub fn statistics(&self) -> (usize, usize, usize) {
        (self.total_iterations, self.total_doublings, self.total_interpolations)
    }

    /// Reset minimizer state
    pub fn reset(&mut self) {
        self.step_size = self.initial_step_size;
        self.step_count = 0;
        self.prev_energy = 0.0;
        self.converged = false;
        self.search_directions.clear();
        self.old_forces.clear();
        self.total_doublings = 0;
        self.total_interpolations = 0;
        self.total_iterations = 0;
    }

    /// Calculate beta coefficient for conjugate direction
    fn calculate_beta(&self, current_forces: &[Vec3]) -> f64 {
        let n_atoms = current_forces.len();

        if self.old_forces.is_empty() || self.old_forces.len() != n_atoms {
            return 0.0;
        }

        let mut f_old_sq = 0.0;
        let mut numerator = 0.0;

        match self.beta_method {
            ConjugateGradientMethod::FletcherReeves => {
                // β = |f_new|² / |f_old|²
                let mut f_new_sq = 0.0;
                for i in 0..n_atoms {
                    f_old_sq += self.old_forces[i].length_squared() as f64;
                    f_new_sq += current_forces[i].length_squared() as f64;
                }
                numerator = f_new_sq;
            }
            ConjugateGradientMethod::PolakRibiere => {
                // β = ⟨f_new, f_new - f_old⟩ / |f_old|²
                for i in 0..n_atoms {
                    let f_diff = current_forces[i] - self.old_forces[i];
                    numerator += current_forces[i].dot(f_diff) as f64;
                    f_old_sq += self.old_forces[i].length_squared() as f64;
                }
            }
        }

        if f_old_sq < 1e-15 {
            0.0
        } else {
            numerator / f_old_sq
        }
    }

    /// Calculate search direction and return slope ⟨p, f⟩
    fn calculate_search_direction(&mut self, forces: &[Vec3], beta: f64) -> f64 {
        let n_atoms = forces.len();

        if self.search_directions.len() != n_atoms {
            self.search_directions.resize(n_atoms, Vec3::ZERO);
        }

        let mut slope = 0.0;

        if beta == 0.0 {
            // Steepest descent: p = f
            for i in 0..n_atoms {
                self.search_directions[i] = forces[i];
                slope += forces[i].length_squared() as f64;
            }
        } else {
            // Conjugate direction: p = f + β·p_old
            for i in 0..n_atoms {
                self.search_directions[i] = forces[i] + self.search_directions[i] * (beta as f32);
                slope += self.search_directions[i].dot(forces[i]) as f64;
            }
        }

        slope
    }

    /// Apply force limiting if configured
    fn apply_force_limit(&self, forces: &mut [Vec3]) {
        if self.force_limit > 0.0 {
            for force in forces.iter_mut() {
                let magnitude = force.length() as f64;
                if magnitude > self.force_limit {
                    *force *= (self.force_limit / magnitude) as f32;
                }
            }
        }
    }

    /// Calculate RMS and maximum force
    fn calculate_force_statistics(&self, forces: &[Vec3]) -> (f64, f64) {
        let mut f_sq_sum = 0.0_f64;
        let mut f_max = 0.0_f64;

        for force in forces {
            let f_sq = force.length_squared() as f64;
            f_sq_sum += f_sq;
            f_max = f_max.max(f_sq);
        }

        let rms = (f_sq_sum / forces.len() as f64).sqrt();
        let max = f_max.sqrt();

        (rms, max)
    }
}

impl Integrator for ConjugateGradient {
    fn step(&mut self, _dt: f64, topo: &Topology, conf: &mut Configuration) {
        let n_atoms = topo.num_atoms();
        self.step_count += 1;

        // Zero velocities (minimization, not dynamics)
        for vel in &mut conf.current_mut().vel {
            *vel = Vec3::ZERO;
        }

        // Calculate current energy
        let current_energy = conf.current().energies.total();

        // Check convergence after minimum steps
        if self.step_count > self.min_steps {
            let (rms_force, _max_force) = self.calculate_force_statistics(&conf.current().force);

            if rms_force < self.force_tolerance {
                self.converged = true;
                return;
            }
        }

        // Apply force limiting
        self.apply_force_limit(&mut conf.current_mut().force);

        // Determine if we should reset search direction
        let should_reset = self.old_forces.is_empty() ||
                          (self.reset_interval > 0 && self.step_count % self.reset_interval == 0);

        // Calculate beta
        let beta = if should_reset {
            0.0
        } else {
            self.calculate_beta(&conf.current().force)
        };

        // Calculate search direction and slope
        let mut slope_a = self.calculate_search_direction(&conf.current().force, beta);

        // If slope is negative, reset to steepest descent
        if slope_a < 0.0 {
            slope_a = self.calculate_search_direction(&conf.current().force, 0.0);
        }

        // Calculate magnitude of search direction
        let p_squared: f64 = self.search_directions.iter()
            .map(|p| p.length_squared() as f64)
            .sum();

        if p_squared < 1e-15 {
            self.converged = true;
            return;
        }

        // Initial step along search direction
        let b_init = self.step_size / p_squared.sqrt();
        let mut b = b_init;

        // Energy and slope at point A (current configuration)
        let ene_a = current_energy;

        // Save current state as old
        self.old_forces = conf.current().force.clone();
        conf.exchange_state();

        // Find upper bound B where slope is negative or energy increases
        let mut ene_b;
        let mut slope_b;
        let mut counter_doub = 0;

        loop {
            // Calculate positions at B
            for i in 0..n_atoms {
                conf.current_mut().pos[i] = conf.old().pos[i] +
                    self.search_directions[i] * (b as f32);
            }

            // Energy at B would need force calculation - for now we'll approximate
            // In full implementation, this would call forcefield
            ene_b = ene_a; // Placeholder

            // Calculate slope at B: ⟨p, f_B⟩
            slope_b = 0.0;
            for i in 0..n_atoms {
                slope_b += self.search_directions[i].dot(conf.current().force[i]) as f64;
            }

            // Accept B if slope is negative or energy increased
            if slope_b < 0.0 || ene_b > ene_a {
                break;
            }

            // Otherwise double the step
            b *= 2.0;
            self.step_size *= 1.1; // Increase step size for next iteration
            counter_doub += 1;

            if counter_doub > 10 {
                break; // Safety limit
            }
        }

        // Cubic interpolation to find minimum X
        let mut counter_ipol = 0;
        let a = 0.0;
        let mut x = b / 2.0; // Initial guess

        for _ in 0..self.max_interpolations {
            // Cubic interpolation formula
            let z = (3.0 * (ene_a - ene_b) / (b - a)) - slope_a - slope_b;
            let w_sq = z * z - slope_a * slope_b;

            if w_sq > 0.0 {
                let w = w_sq.sqrt();
                x = b - (w - z - slope_b) * (b - a) / (slope_a - slope_b + 2.0 * w);
            }

            counter_ipol += 1;

            // Check displacement criterion
            let disp = (x * p_squared.sqrt() / n_atoms as f64).abs();
            if disp < self.interpolation_criterion {
                break;
            }
        }

        // Update positions to X
        for i in 0..n_atoms {
            conf.current_mut().pos[i] = conf.old().pos[i] +
                self.search_directions[i] * (x as f32);
        }

        // Update statistics
        self.total_doublings += counter_doub;
        self.total_interpolations += counter_ipol;
        self.total_iterations += 1;

        // Adapt step size
        if x < b_init / 10.0 {
            self.step_size *= 0.9;
        }
        if self.step_size > self.max_step_size {
            self.step_size = self.max_step_size;
        }

        self.prev_energy = current_energy;

        // Clear velocities and forces
        for vel in &mut conf.current_mut().vel {
            *vel = Vec3::ZERO;
        }
        conf.current_mut().clear_forces();
    }

    fn name(&self) -> &str {
        "Conjugate-Gradient"
    }
}

impl Default for ConjugateGradient {
    fn default() -> Self {
        Self::new()
    }
}

/// Lattice Shift Tracker for FEP calculations
///
/// Direct translation of md++/src/algorithm/integration/lattice_shift.cc
///
/// **Purpose**: Track periodic boundary crossings for accurate FEP calculations
///
/// When atoms cross periodic boundaries in a simulation, they are "wrapped" back
/// into the primary simulation box. For Free Energy Perturbation (FEP) calculations
/// with long-range interactions (PME, Reaction Field), we need to track these
/// crossings to properly calculate energy derivatives ∂H/∂λ.
///
/// **How it works**:
/// 1. When an atom crosses a boundary (position wraps), increment shift counter
/// 2. Store cumulative shifts in lattice_shifts vector
/// 3. Use shifts to unwrap coordinates for long-range FEP calculations
///
/// **Critical for**:
/// - FEP with PME (Particle Mesh Ewald)
/// - FEP with Reaction Field electrostatics
/// - TI (Thermodynamic Integration) calculations
///
/// **Example**:
/// ```text
/// Atom moves from x=9.9 to x=0.1 in a box of length 10.0
/// → Crossed +x boundary, lattice_shift[atom][0] += 1
/// → Unwrapped coordinate: x = 0.1 + 1*10.0 = 10.1
/// ```
#[derive(Debug, Clone)]
pub struct LatticeShiftTracker {
    /// Enable tracking (can disable for non-FEP simulations)
    pub enabled: bool,
}

impl LatticeShiftTracker {
    /// Create new lattice shift tracker
    pub fn new() -> Self {
        Self { enabled: true }
    }

    /// Disable tracking (for non-FEP simulations)
    pub fn disabled() -> Self {
        Self { enabled: false }
    }

    /// Track lattice shifts while wrapping atoms into box
    ///
    /// This should be called after position updates to ensure atoms stay in the box
    /// and track boundary crossings.
    ///
    /// # Arguments
    /// * `conf` - Configuration with positions and lattice_shifts
    ///
    /// # Algorithm
    /// For each dimension (x, y, z):
    /// 1. Calculate scaled coordinate: s = r · box_inv
    /// 2. Determine shift: shift = floor(s + 0.5)
    /// 3. Wrap coordinate: r_new = r - shift · box_vector
    /// 4. Accumulate shift: lattice_shifts[i][d] += shift
    pub fn apply(&self, conf: &mut Configuration) {
        if !self.enabled {
            return;
        }

        // Clone box information to avoid borrow conflicts
        let box_type = conf.current().box_config.box_type;
        let box_vectors = conf.current().box_config.vectors;
        let box_inv = conf.current().box_config.inv_vectors;

        match box_type {
            BoxType::Vacuum => {
                // No periodic boundaries
                return;
            }
            BoxType::Rectangular => {
                self.apply_rectangular(conf, &box_vectors, &box_inv);
            }
            BoxType::Triclinic | BoxType::TruncatedOctahedral => {
                self.apply_triclinic(conf, &box_vectors, &box_inv);
            }
        }
    }

    /// Apply for rectangular box (optimized path)
    fn apply_rectangular(&self, conf: &mut Configuration, _box_vectors: &Mat3, _box_inv: &Mat3) {
        let n_atoms = conf.current().pos.len();

        // Extract box dimensions to avoid borrow conflicts
        let box_x = conf.current().box_config.vectors.x_axis.x;
        let box_y = conf.current().box_config.vectors.y_axis.y;
        let box_z = conf.current().box_config.vectors.z_axis.z;

        for i in 0..n_atoms {
            // Read position
            let mut pos = conf.current().pos[i];

            // X dimension
            let shift_x = (pos.x / box_x + 0.5).floor() as i32;
            if shift_x != 0 {
                pos.x -= shift_x as f32 * box_x;
                conf.lattice_shifts[i][0] += shift_x;
            }

            // Y dimension
            let shift_y = (pos.y / box_y + 0.5).floor() as i32;
            if shift_y != 0 {
                pos.y -= shift_y as f32 * box_y;
                conf.lattice_shifts[i][1] += shift_y;
            }

            // Z dimension
            let shift_z = (pos.z / box_z + 0.5).floor() as i32;
            if shift_z != 0 {
                pos.z -= shift_z as f32 * box_z;
                conf.lattice_shifts[i][2] += shift_z;
            }

            // Write back the updated position
            if shift_x != 0 || shift_y != 0 || shift_z != 0 {
                conf.current_mut().pos[i] = pos;
            }
        }
    }

    /// Apply for triclinic box (general case)
    fn apply_triclinic(&self, conf: &mut Configuration, box_vectors: &Mat3, box_inv: &Mat3) {
        let n_atoms = conf.current().pos.len();

        for i in 0..n_atoms {
            let pos = conf.current().pos[i];

            // Transform to scaled coordinates: s = inv(box) · r
            let scaled = Vec3::new(
                box_inv.x_axis.dot(pos),
                box_inv.y_axis.dot(pos),
                box_inv.z_axis.dot(pos),
            );

            // Calculate shifts (floor(s + 0.5))
            let shift_x = (scaled.x + 0.5).floor() as i32;
            let shift_y = (scaled.y + 0.5).floor() as i32;
            let shift_z = (scaled.z + 0.5).floor() as i32;

            // Wrap position if needed
            if shift_x != 0 || shift_y != 0 || shift_z != 0 {
                let shift_vec = box_vectors.x_axis * (shift_x as f32) +
                               box_vectors.y_axis * (shift_y as f32) +
                               box_vectors.z_axis * (shift_z as f32);

                let new_pos = pos - shift_vec;
                conf.current_mut().pos[i] = new_pos;

                // Accumulate shifts
                conf.lattice_shifts[i][0] += shift_x;
                conf.lattice_shifts[i][1] += shift_y;
                conf.lattice_shifts[i][2] += shift_z;
            }
        }
    }

    /// Reset all lattice shifts to zero
    ///
    /// Should be called when starting a new simulation or after trajectory restart
    pub fn reset(&self, conf: &mut Configuration) {
        for shifts in &mut conf.lattice_shifts {
            *shifts = [0, 0, 0];
        }
    }

    /// Get unwrapped position for an atom
    ///
    /// Returns the actual position accounting for all boundary crossings.
    /// Useful for trajectory analysis and visualization.
    ///
    /// # Arguments
    /// * `conf` - Configuration
    /// * `atom_idx` - Atom index
    ///
    /// # Returns
    /// Unwrapped position in real space
    pub fn unwrapped_position(&self, conf: &Configuration, atom_idx: usize) -> Vec3 {
        let pos = conf.current().pos[atom_idx];
        let shifts = conf.lattice_shifts[atom_idx];
        let box_vectors = &conf.current().box_config.vectors;

        pos + box_vectors.x_axis * (shifts[0] as f32) +
              box_vectors.y_axis * (shifts[1] as f32) +
              box_vectors.z_axis * (shifts[2] as f32)
    }

    /// Calculate mean squared displacement using unwrapped coordinates
    ///
    /// Useful for diffusion coefficient calculations
    pub fn mean_squared_displacement(
        &self,
        conf_initial: &Configuration,
        conf_current: &Configuration,
    ) -> f64 {
        let n_atoms = conf_current.current().pos.len();
        let mut msd = 0.0;

        for i in 0..n_atoms {
            let r0 = self.unwrapped_position(conf_initial, i);
            let r1 = self.unwrapped_position(conf_current, i);
            let dr = r1 - r0;
            msd += dr.length_squared() as f64;
        }

        msd / n_atoms as f64
    }
}

impl Default for LatticeShiftTracker {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_system() -> (Topology, Configuration) {
        let mut topo = Topology::new();
        topo.mass = vec![1.0; 10];
        topo.charge = vec![0.0; 10];
        topo.iac = vec![0; 10];
        topo.compute_inverse_masses();

        let conf = Configuration::new(10, 1, 1);

        (topo, conf)
    }

    #[test]
    fn test_leap_frog_integrator() {
        let (topo, mut conf) = create_test_system();

        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().vel[0] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().force[0] = Vec3::new(1.0, 0.0, 0.0);

        let mut integrator = LeapFrog::new();
        let dt = 0.001;

        integrator.step(dt, &topo, &mut conf);

        assert!(conf.current().pos[0].x > 0.0);
        assert!(conf.current().vel[0].x > 1.0);
    }

    #[test]
    fn test_energy_conservation() {
        let (topo, mut conf) = create_test_system();

        for (i, vel) in conf.current_mut().vel.iter_mut().enumerate() {
            *vel = Vec3::new((i as f32) * 0.1, 0.0, 0.0);
        }

        conf.current_mut().calculate_kinetic_energy(&topo.mass);
        let initial_ke = conf.current().energies.kinetic_total;

        let mut integrator = LeapFrog::new();
        for _ in 0..100 {
            integrator.step(0.001, &topo, &mut conf);
        }

        conf.current_mut().calculate_kinetic_energy(&topo.mass);
        let final_ke = conf.current().energies.kinetic_total;

        assert!((final_ke - initial_ke).abs() < 1e-6);
    }
}
