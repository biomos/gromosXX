//! Integration algorithms for molecular dynamics
//!
//! Direct Rust translation of GROMOS integration algorithms from:
//! - md++/src/algorithm/integration/leap_frog.cc
//! - md++/src/algorithm/integration/velocity_verlet.cc
//! - md++/src/algorithm/integration/stochastic.cc

use crate::math::Vec3;
use crate::topology::Topology;
use crate::configuration::Configuration;
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
