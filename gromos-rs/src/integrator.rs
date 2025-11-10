//! Integration algorithms for molecular dynamics
//!
//! Direct Rust translation of GROMOS integration algorithms from:
//! - md++/src/algorithm/integration/leap_frog.cc
//! - md++/src/algorithm/integration/velocity_verlet.cc

use crate::math::Vec3;
use crate::topology::Topology;
use crate::configuration::Configuration;
use rayon::prelude::*;

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
        if self.parallel {
            let old_pos = conf.old().pos.clone();
            let current_vel = conf.current().vel.clone();

            conf.current_mut().pos.par_iter_mut()
                .enumerate()
                .for_each(|(i, pos_new)| {
                    *pos_new = old_pos[i] + current_vel[i] * dt_f32;
                });
        } else {
            for i in 0..n_atoms {
                conf.current_mut().pos[i] = conf.old().pos[i] + conf.current().vel[i] * dt_f32;
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
