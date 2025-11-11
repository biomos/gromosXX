//! Restraint interactions for MD simulations
//!
//! Direct translation of GROMOS restraint interactions from:
//! - md++/src/interaction/special/position_restraint_interaction.cc
//! - md++/src/interaction/special/distance_restraint_interaction.cc

use crate::math::{Vec3, Periodicity};
use crate::configuration::Configuration;

/// Position restraint - restrains atom to reference position
///
/// Direct translation from position_restraint_interaction.cc
///
/// Energy: E = 0.5 * k * |r - r_ref|²
/// Force: F = -k * (r - r_ref)
#[derive(Debug, Clone)]
pub struct PositionRestraint {
    /// Atom index
    pub atom: usize,
    /// Reference position
    pub reference_pos: Vec3,
    /// Force constant (kJ mol⁻¹ nm⁻²)
    pub force_constant: f64,
}

impl PositionRestraint {
    pub fn new(atom: usize, reference_pos: Vec3, force_constant: f64) -> Self {
        Self {
            atom,
            reference_pos,
            force_constant,
        }
    }

    /// Calculate restraint energy and forces
    pub fn calculate(
        &self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
    ) -> f64 {
        // Get displacement vector using periodic boundary conditions
        // nearest_image returns (reference - current), but we need (current - reference)
        let displacement = periodicity.nearest_image(
            self.reference_pos,
            conf.current().pos[self.atom]
        );

        let dist_sq = displacement.length_squared() as f64;
        let dist = (dist_sq).sqrt();

        // Energy: 0.5 * k * r²
        let energy = 0.5 * self.force_constant * dist_sq;

        // Force: F = -k * (r - r0) = -k * displacement
        let force_magnitude = -self.force_constant;
        let force = displacement * (force_magnitude as f32);

        conf.current_mut().force[self.atom] += force;

        energy
    }
}

/// Collection of position restraints
#[derive(Debug, Clone, Default)]
pub struct PositionRestraints {
    pub restraints: Vec<PositionRestraint>,
}

impl PositionRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: PositionRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(
        &self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
    ) -> f64 {
        self.restraints.iter()
            .map(|r| r.calculate(conf, periodicity))
            .sum()
    }
}

/// Distance restraint - restrains distance between two atoms
///
/// Direct translation from distance_restraint_interaction.cc
///
/// Supports:
/// - Harmonic/linear force function
/// - Time-averaging for NOE restraints
/// - Repulsive/attractive/harmonic modes
#[derive(Debug, Clone)]
pub struct DistanceRestraint {
    /// First atom index
    pub atom_i: usize,
    /// Second atom index
    pub atom_j: usize,
    /// Target distance (nm)
    pub r0: f64,
    /// Force constant (kJ mol⁻¹ nm⁻²)
    pub force_constant: f64,
    /// Linear region width (nm)
    pub r_linear: f64,
    /// Restraint type: -1 (repulsive), 0 (harmonic), +1 (attractive)
    pub restraint_type: i32,
    /// Time-averaged distance (for NOE restraints, r⁻³ averaging)
    pub time_averaged_distance: Option<f64>,
    /// Time constant for averaging (ps)
    pub tau: Option<f64>,
}

impl DistanceRestraint {
    pub fn new(
        atom_i: usize,
        atom_j: usize,
        r0: f64,
        force_constant: f64,
    ) -> Self {
        Self {
            atom_i,
            atom_j,
            r0,
            force_constant,
            r_linear: 0.0,  // Default: pure harmonic
            restraint_type: 0,  // Harmonic by default
            time_averaged_distance: None,
            tau: None,
        }
    }

    /// Set linear region width
    pub fn with_linear_region(mut self, r_linear: f64) -> Self {
        self.r_linear = r_linear;
        self
    }

    /// Set restraint type: -1 (repulsive), 0 (harmonic), +1 (attractive)
    pub fn with_type(mut self, restraint_type: i32) -> Self {
        self.restraint_type = restraint_type;
        self
    }

    /// Enable time averaging with exponential decay
    pub fn with_time_averaging(mut self, tau: f64) -> Self {
        self.tau = Some(tau);
        self.time_averaged_distance = Some(self.r0);
        self
    }

    /// Calculate restraint energy and forces
    ///
    /// Direct translation from distance_restraint_interaction.cc
    pub fn calculate(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        // Get distance vector using periodic boundary conditions
        let r_vec = periodicity.nearest_image(
            conf.current().pos[self.atom_j],
            conf.current().pos[self.atom_i]
        );

        let mut distance = r_vec.length() as f64;

        // Time averaging (NOE r⁻³ averaging)
        if let (Some(tau), Some(avg_dist)) = (self.tau, self.time_averaged_distance) {
            let exp_term = (-dt / tau).exp();
            let r_inv_3 = 1.0 / (distance * distance * distance);
            let avg_r_inv_3 = 1.0 / (avg_dist * avg_dist * avg_dist);

            // Update average: <r⁻³> = (1-e^(-dt/τ)) * r⁻³ + e^(-dt/τ) * <r⁻³>_old
            let new_avg_r_inv_3 = (1.0 - exp_term) * r_inv_3 + exp_term * avg_r_inv_3;
            let new_avg_dist = new_avg_r_inv_3.powf(-1.0/3.0);

            self.time_averaged_distance = Some(new_avg_dist);
            distance = new_avg_dist;
        }

        let delta = distance - self.r0;

        // Check if restraint is violated based on type
        let violated = match self.restraint_type {
            -1 => distance < self.r0,   // Repulsive: violated if too close
            0  => true,                  // Harmonic: always active
            1  => distance > self.r0,   // Attractive: violated if too far
            _ => true,
        };

        if !violated {
            return 0.0;  // No restraint force
        }

        // Calculate energy and force
        let (energy, force_magnitude) = if delta.abs() < self.r_linear || self.r_linear == 0.0 {
            // Harmonic region: E = 0.5 * k * (r - r0)²
            let energy = 0.5 * self.force_constant * delta * delta;
            let force_mag = -self.force_constant * delta;
            (energy, force_mag)
        } else {
            // Linear region (prevents unrealistic forces at large distances)
            let sign = if delta > 0.0 { 1.0 } else { -1.0 };
            let energy = self.force_constant * self.r_linear *
                        (delta.abs() - 0.5 * self.r_linear);
            let force_mag = -sign * self.force_constant * self.r_linear;
            (energy, force_mag)
        };

        // Apply force along bond direction
        if distance > 1e-10 {
            let force = r_vec * ((force_magnitude / distance) as f32);
            conf.current_mut().force[self.atom_i] -= force;
            conf.current_mut().force[self.atom_j] += force;
        }

        energy
    }
}

/// Collection of distance restraints
#[derive(Debug, Clone, Default)]
pub struct DistanceRestraints {
    pub restraints: Vec<DistanceRestraint>,
}

impl DistanceRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: DistanceRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        self.restraints.iter_mut()
            .map(|r| r.calculate(conf, periodicity, dt))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_position_restraint() {
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(1.0, 0.0, 0.0);

        let restraint = PositionRestraint::new(
            0,
            Vec3::new(0.0, 0.0, 0.0),
            100.0
        );

        let periodicity = Periodicity::Rectangular(
            crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0))
        );

        let energy = restraint.calculate(&mut conf, &periodicity);

        // E = 0.5 * 100 * 1² = 50 kJ/mol
        assert!((energy - 50.0).abs() < 1e-6);

        // Force should point toward origin: F = -k * r = -100 * 1 = -100
        // (pointing in negative x direction)
        assert!((conf.current().force[0].x + 100.0).abs() < 1e-3);
    }

    #[test]
    fn test_distance_restraint_harmonic() {
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(2.0, 0.0, 0.0);

        let mut restraint = DistanceRestraint::new(0, 1, 1.0, 100.0);

        let periodicity = Periodicity::Rectangular(
            crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0))
        );

        let energy = restraint.calculate(&mut conf, &periodicity, 0.001);

        // Distance = 2.0, r0 = 1.0, delta = 1.0
        // E = 0.5 * 100 * 1² = 50 kJ/mol
        assert!((energy - 50.0).abs() < 1e-6);

        // Forces should be equal and opposite
        let force_i = conf.current().force[0];
        let force_j = conf.current().force[1];
        assert!((force_i.x + force_j.x).abs() < 1e-3);
    }

    #[test]
    fn test_distance_restraint_linear() {
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(5.0, 0.0, 0.0);

        let mut restraint = DistanceRestraint::new(0, 1, 1.0, 100.0)
            .with_linear_region(0.1);

        let periodicity = Periodicity::Rectangular(
            crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0))
        );

        let energy = restraint.calculate(&mut conf, &periodicity, 0.001);

        // Distance = 5.0, r0 = 1.0, delta = 4.0
        // Linear region since |delta| > 0.1
        // E = k * r_linear * (|delta| - 0.5*r_linear)
        // E = 100 * 0.1 * (4.0 - 0.05) = 10 * 3.95 = 39.5
        assert!((energy - 39.5).abs() < 1e-6);
    }

    #[test]
    fn test_restraint_types() {
        let mut conf = Configuration::new(2, 1, 1);
        let periodicity = Periodicity::Rectangular(
            crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0))
        );

        // Repulsive restraint: only active when distance < r0
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.5, 0.0, 0.0);
        let mut restraint = DistanceRestraint::new(0, 1, 1.0, 100.0)
            .with_type(-1);  // Repulsive
        let energy = restraint.calculate(&mut conf, &periodicity, 0.001);
        assert!(energy > 0.0);  // Should have energy

        // Same restraint but distance > r0 - no energy
        conf.current_mut().pos[1] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().clear_forces();
        let energy = restraint.calculate(&mut conf, &periodicity, 0.001);
        assert_eq!(energy, 0.0);  // Should have no energy
    }
}
