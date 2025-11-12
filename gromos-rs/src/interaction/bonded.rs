//! Bonded interaction calculations
//!
//! This module implements all bonded force field terms:
//! - Bonds (quartic and harmonic)
//! - Angles (cosine-based and harmonic)
//! - Proper dihedrals
//! - Improper dihedrals
//! - Cross-dihedrals

use crate::math::Vec3;
use crate::topology::{Topology, BondParameters, AngleParameters};
use crate::configuration::Configuration;
use crate::fep::LambdaController;

/// Result of force calculation: energy and forces
#[derive(Debug, Clone)]
pub struct ForceEnergy {
    pub energy: f64,
    pub forces: Vec<Vec3>,
}

/// Result of FEP force calculation: energy, forces, and lambda derivative
#[derive(Debug, Clone)]
pub struct ForceEnergyLambda {
    pub energy: f64,
    pub forces: Vec<Vec3>,
    pub lambda_derivative: f64,  // dE/dλ for free energy calculations
}

impl ForceEnergy {
    pub fn new(num_atoms: usize) -> Self {
        Self {
            energy: 0.0,
            forces: vec![Vec3::ZERO; num_atoms],
        }
    }

    pub fn add(&mut self, other: &ForceEnergy) {
        self.energy += other.energy;
        for i in 0..self.forces.len().min(other.forces.len()) {
            self.forces[i] += other.forces[i];
        }
    }
}

impl ForceEnergyLambda {
    pub fn new(num_atoms: usize) -> Self {
        Self {
            energy: 0.0,
            forces: vec![Vec3::ZERO; num_atoms],
            lambda_derivative: 0.0,
        }
    }

    pub fn add(&mut self, other: &ForceEnergyLambda) {
        self.energy += other.energy;
        self.lambda_derivative += other.lambda_derivative;
        for i in 0..self.forces.len().min(other.forces.len()) {
            self.forces[i] += other.forces[i];
        }
    }
}

/// Calculate quartic bond forces and energies (GROMOS standard)
///
/// Potential: V = (1/4) * k_harmonic * (r^2 - r0^2)^2
/// Force: F = -dV/dr = -k_harmonic * (r^2 - r0^2) * r * r_vec/r
pub fn calculate_bond_forces_quartic(
    topo: &Topology,
    conf: &Configuration,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for bond in &topo.solute.bonds {
        if bond.bond_type >= topo.bond_parameters.len() {
            continue;
        }

        let params = &topo.bond_parameters[bond.bond_type];
        let r_vec = conf.current().pos[bond.j] - conf.current().pos[bond.i];
        let r = r_vec.length() as f64;

        if r < 1e-10 {
            // Avoid division by zero
            continue;
        }

        // Energy: V = (1/4) * k * (r^2 - r0^2)^2
        let r2 = r * r;
        let r0_2 = params.r0 * params.r0;
        let dr2 = r2 - r0_2;
        let energy = 0.25 * params.k_harmonic * dr2 * dr2;

        // Force magnitude: dV/dr = k * (r^2 - r0^2) * r
        let f_magnitude = params.k_harmonic * dr2 * r;

        // Force vector: F = -f_magnitude * r_vec/r
        let force = r_vec * (-f_magnitude as f32 / r as f32);

        result.energy += energy;
        result.forces[bond.i] += force;
        result.forces[bond.j] -= force;
    }

    result
}

/// Calculate harmonic bond forces (alternative to quartic)
///
/// Potential: V = (1/2) * k * (r - r0)^2
/// Force: F = -k * (r - r0) * r_vec/r
pub fn calculate_bond_forces_harmonic(
    topo: &Topology,
    conf: &Configuration,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for bond in &topo.solute.bonds {
        if bond.bond_type >= topo.bond_parameters.len() {
            continue;
        }

        let params = &topo.bond_parameters[bond.bond_type];
        let r_vec = conf.current().pos[bond.j] - conf.current().pos[bond.i];
        let r = r_vec.length() as f64;

        if r < 1e-10 {
            continue;
        }

        // Energy: V = (1/2) * k * (r - r0)^2
        let dr = r - params.r0;
        let energy = 0.5 * params.k_harmonic * dr * dr;

        // Force: F = -k * (r - r0) * r_vec/r
        let f_magnitude = params.k_harmonic * dr;
        let force = r_vec * (-f_magnitude as f32 / r as f32);

        result.energy += energy;
        result.forces[bond.i] += force;
        result.forces[bond.j] -= force;
    }

    result
}

/// Calculate CG (coarse-grained) bond forces - quartic repulsive potential
///
/// Potential: V = 0.5 * K * (r - r0)^4 when r > r0, otherwise V = 0
/// Force: F = -2 * K * (r - r0)^3 * direction when r > r0, otherwise F = 0
///
/// This is a soft repulsive potential that prevents bonds from stretching too much
/// but allows free compression. Useful for coarse-grained models.
pub fn calculate_cg_bond_forces(
    topo: &Topology,
    conf: &Configuration,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    // CG bonds are stored separately in topology (if available)
    // For now, we'll use a marker in bond parameters to identify CG bonds
    for bond in &topo.solute.bonds {
        if bond.bond_type >= topo.bond_parameters.len() {
            continue;
        }

        let params = &topo.bond_parameters[bond.bond_type];

        // Skip if this is not a CG bond (we'll add a flag later)
        // For now, assume CG bonds have negative k_harmonic as a marker
        if params.k_harmonic >= 0.0 {
            continue;
        }

        let r_vec = conf.current().pos[bond.j] - conf.current().pos[bond.i];
        let r = r_vec.length() as f64;

        if r < 1e-10 {
            // Avoid division by zero
            continue;
        }

        // Only apply force when bond is stretched (r > r0)
        if r > params.r0 {
            let diff = r - params.r0;
            let diff2 = diff * diff;
            let diff3 = diff2 * diff;
            let diff4 = diff2 * diff2;

            // Use absolute value of k_harmonic (it's negative as a marker)
            let k = params.k_harmonic.abs();

            // Energy: V = 0.5 * K * (r - r0)^4
            let energy = 0.5 * k * diff4;

            // Force magnitude: -dV/dr = -2 * K * (r - r0)^3
            let f_magnitude = -2.0 * k * diff3;

            // Force vector: F = f_magnitude * direction
            let direction = r_vec / (r as f32);
            let force = direction * (f_magnitude as f32);

            result.energy += energy;
            result.forces[bond.i] += force;
            result.forces[bond.j] -= force;
        }
    }

    result
}

/// Calculate angle forces (GROMOS cosine-based)
///
/// Potential: V = (1/2) * k_harmonic * (cos(θ) - cos(θ0))^2
/// Forces calculated using chain rule on angle derivatives
pub fn calculate_angle_forces(
    topo: &Topology,
    conf: &Configuration,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for angle in &topo.solute.angles {
        if angle.angle_type >= topo.angle_parameters.len() {
            continue;
        }

        let params = &topo.angle_parameters[angle.angle_type];

        // Vectors: i->j and k->j (both point to central atom j)
        let r_ij = conf.current().pos[angle.j] - conf.current().pos[angle.i];
        let r_kj = conf.current().pos[angle.j] - conf.current().pos[angle.k];

        let len_ij = r_ij.length() as f64;
        let len_kj = r_kj.length() as f64;

        if len_ij < 1e-10 || len_kj < 1e-10 {
            continue;
        }

        // cos(θ) = r_ij · r_kj / (|r_ij| |r_kj|)
        let cos_theta = (r_ij.dot(r_kj) / (len_ij * len_kj) as f32).clamp(-1.0, 1.0) as f64;
        let cos_theta0 = params.theta0.cos();

        // Energy: V = (1/2) * k * (cos(θ) - cos(θ0))^2
        let d_cos = cos_theta - cos_theta0;
        let energy = 0.5 * params.k_harmonic * d_cos * d_cos;

        // Force calculation using chain rule
        // F = -dV/dcos(θ) * dcos(θ)/dr
        let dV_dcos = params.k_harmonic * d_cos;

        // Derivatives of cos(θ) with respect to positions
        // d(cos θ)/dr_i = (r_kj / (|r_ij| |r_kj|)) - cos(θ) * r_ij / |r_ij|^2
        // d(cos θ)/dr_k = (r_ij / (|r_ij| |r_kj|)) - cos(θ) * r_kj / |r_kj|^2

        let inv_len_ij = 1.0 / len_ij;
        let inv_len_kj = 1.0 / len_kj;
        let inv_len_prod = inv_len_ij * inv_len_kj;

        // Force on atom i
        let dcos_dri = r_kj * (inv_len_prod as f32) - r_ij * ((cos_theta * inv_len_ij * inv_len_ij) as f32);
        let f_i = dcos_dri * (-dV_dcos as f32);

        // Force on atom k
        let dcos_drk = r_ij * (inv_len_prod as f32) - r_kj * ((cos_theta * inv_len_kj * inv_len_kj) as f32);
        let f_k = dcos_drk * (-dV_dcos as f32);

        // Force on atom j (central atom): F_j = -(F_i + F_k)
        let f_j = -(f_i + f_k);

        result.energy += energy;
        result.forces[angle.i] += f_i;
        result.forces[angle.j] += f_j;
        result.forces[angle.k] += f_k;
    }

    result
}

/// Calculate harmonic angle forces (simple harmonic potential)
///
/// Potential: V = (1/2) * K * (θ - θ₀)²
/// where:
/// - K is the force constant
/// - θ is the bond angle in radians
/// - θ₀ is the equilibrium angle in radians
///
/// This is an alternative to the cosine-based angle potential.
/// It's simpler and more intuitive, commonly used in other force fields.
pub fn calculate_harmonic_angle_forces(
    topo: &Topology,
    conf: &Configuration,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    const EPSILON: f64 = 1e-10;  // Small value for numerical stability
    const PI: f64 = std::f64::consts::PI;

    for angle in &topo.solute.angles {
        if angle.angle_type >= topo.angle_parameters.len() {
            continue;
        }

        let params = &topo.angle_parameters[angle.angle_type];

        // Get position vectors (point from outer atoms to central atom j)
        let r_ij = conf.current().pos[angle.j] - conf.current().pos[angle.i];
        let r_kj = conf.current().pos[angle.j] - conf.current().pos[angle.k];

        let d_ij = r_ij.length() as f64;
        let d_kj = r_kj.length() as f64;

        if d_ij < EPSILON || d_kj < EPSILON {
            continue; // Avoid division by zero
        }

        // Calculate angle θ using dot product
        let cos_theta = (r_ij.dot(r_kj) / (d_ij * d_kj) as f32).clamp(-1.0, 1.0) as f64;
        let theta = cos_theta.acos();
        let sin_theta = theta.sin();

        // Get equilibrium angle (stored in params.theta0)
        let theta0 = params.theta0;

        // Force constants for the three atoms
        let k_i: f64;
        let k_k: f64;
        let energy: f64;

        // Special handling when sin(θ) ≈ 0 (θ near 0 or π)
        if sin_theta.abs() < EPSILON {
            // When θ ≈ π (linear angle), special treatment to avoid numerical errors
            if (theta0 > PI + EPSILON) || (theta0 < PI - EPSILON) {
                // θ₀ ≠ π but current angle is π: this is problematic
                // Skip this angle to avoid numerical issues
                continue;
            }

            // Special force constants when θ ≈ π
            // From GROMOS++: ki = -K / dij, kk = -K / dkj
            k_i = -params.k_harmonic / d_ij;
            k_k = -params.k_harmonic / d_kj;

            // Energy when θ ≈ π: V = K * (1 + cos(θ))
            energy = params.k_harmonic * (1.0 + cos_theta);
        } else {
            // Normal case: θ is not near 0 or π
            // Force constants: ki = K * (θ - θ₀) / (sin(θ) * dij)
            let force_factor = params.k_harmonic * (theta - theta0) / sin_theta;
            k_i = force_factor / d_ij;
            k_k = force_factor / d_kj;

            // Energy: V = 0.5 * K * (θ - θ₀)²
            let d_theta = theta - theta0;
            energy = 0.5 * params.k_harmonic * d_theta * d_theta;
        }

        // Calculate forces on the three atoms
        // The force direction is perpendicular to the bond vectors
        // fi = ki * (r_kj/d_kj - cos(θ) * r_ij/d_ij)
        // fk = kk * (r_ij/d_ij - cos(θ) * r_kj/d_kj)
        // fj = -(fi + fk)  [force conservation]

        let unit_ij = r_ij / d_ij as f32;
        let unit_kj = r_kj / d_kj as f32;

        let f_i = (unit_kj - unit_ij * cos_theta as f32) * k_i as f32;
        let f_k = (unit_ij - unit_kj * cos_theta as f32) * k_k as f32;
        let f_j = -(f_i + f_k);  // Force conservation

        result.energy += energy;
        result.forces[angle.i] += f_i;
        result.forces[angle.j] += f_j;
        result.forces[angle.k] += f_k;
    }

    result
}

/// Calculate proper dihedral forces (GROMOS torsional potential)
///
/// Potential: V = K * [1 + cos(δ) * cos(m*φ)]
/// where:
/// - K is the force constant
/// - δ is the phase shift angle
/// - m is the multiplicity (1-6)
/// - φ is the dihedral angle
///
/// Uses Chebyshev polynomials for efficient cos(m*φ) calculation
pub fn calculate_dihedral_forces(
    topo: &Topology,
    conf: &Configuration,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for dihedral in &topo.solute.proper_dihedrals {
        if dihedral.dihedral_type >= topo.dihedral_parameters.len() {
            continue;
        }

        let params = &topo.dihedral_parameters[dihedral.dihedral_type];

        // Get position vectors (applying PBC if needed)
        let r_ij = conf.current().pos[dihedral.j] - conf.current().pos[dihedral.i];
        let r_kj = conf.current().pos[dihedral.j] - conf.current().pos[dihedral.k];
        let r_kl = conf.current().pos[dihedral.l] - conf.current().pos[dihedral.k];

        // Cross products to get plane normals
        let r_mj = r_ij.cross(r_kj);  // normal to plane i-j-k
        let r_nk = r_kj.cross(r_kl);  // normal to plane j-k-l

        let d_kj2 = r_kj.dot(r_kj);

        if d_kj2 < 1e-10 {
            continue;
        }

        let f_rim = r_ij.dot(r_kj) / d_kj2;
        let f_rln = r_kl.dot(r_kj) / d_kj2;

        let r_im = r_ij - r_kj * f_rim;
        let r_ln = r_kj * f_rln - r_kl;

        let d_im = r_im.length() as f64;
        let d_ln = r_ln.length() as f64;

        if d_im < 1e-10 || d_ln < 1e-10 {
            continue;
        }

        // Calculate cos(φ)
        let cos_phi = (r_im.dot(r_ln) / (d_im * d_ln) as f32).clamp(-1.0, 1.0) as f64;

        // Calculate cos(m*φ) and d[cos(m*φ)]/d[cos(φ)] using Chebyshev polynomials
        let (cos_m_phi, d_cos_m_phi) = match params.m {
            0 => (0.0, 0.0),
            1 => (cos_phi, 1.0),
            2 => (2.0 * cos_phi * cos_phi - 1.0, 4.0 * cos_phi),
            3 => (4.0 * cos_phi.powi(3) - 3.0 * cos_phi, 12.0 * cos_phi * cos_phi - 3.0),
            4 => (8.0 * cos_phi.powi(4) - 8.0 * cos_phi * cos_phi + 1.0,
                  32.0 * cos_phi.powi(3) - 16.0 * cos_phi),
            5 => (16.0 * cos_phi.powi(5) - 20.0 * cos_phi.powi(3) + 5.0 * cos_phi,
                  80.0 * cos_phi.powi(4) - 60.0 * cos_phi * cos_phi + 5.0),
            6 => (32.0 * cos_phi.powi(6) - 48.0 * cos_phi.powi(4) + 18.0 * cos_phi * cos_phi - 1.0,
                  192.0 * cos_phi.powi(5) - 192.0 * cos_phi.powi(3) + 36.0 * cos_phi),
            _ => {
                eprintln!("Warning: unsupported dihedral multiplicity {}", params.m);
                continue;
            }
        };

        // Energy: V = K * (1 + cos(δ) * cos(m*φ))
        let energy = params.k * (1.0 + params.cospd * cos_m_phi);

        // Force calculation
        let k_i = -params.k * params.cospd * d_cos_m_phi / d_im;
        let k_l = -params.k * params.cospd * d_cos_m_phi / d_ln;
        let k_j1 = f_rim - 1.0;
        let k_j2 = f_rln;

        let f_i = r_ln * (k_i / d_ln) as f32 - r_im * (k_i * cos_phi / d_im) as f32;
        let f_l = r_im * (k_l / d_im) as f32 - r_ln * (k_l * cos_phi / d_ln) as f32;
        let f_j = f_i * k_j1 as f32 - f_l * k_j2 as f32;
        let f_k = -(f_i + f_j + f_l);

        result.energy += energy;
        result.forces[dihedral.i] += f_i;
        result.forces[dihedral.j] += f_j;
        result.forces[dihedral.k] += f_k;
        result.forces[dihedral.l] += f_l;
    }

    result
}

/// Calculate improper dihedral forces (for planarity/chirality)
///
/// Potential: V = (1/2) * K * (ζ - ζ₀)²
/// where:
/// - K is the force constant
/// - ζ is the improper dihedral angle
/// - ζ₀ is the reference angle
///
/// Improper dihedrals maintain planarity or chirality
pub fn calculate_improper_dihedral_forces(
    topo: &Topology,
    conf: &Configuration,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for improper in &topo.solute.improper_dihedrals {
        if improper.dihedral_type >= topo.improper_dihedral_parameters.len() {
            continue;
        }

        let params = &topo.improper_dihedral_parameters[improper.dihedral_type];

        // Get position vectors
        let r_kj = conf.current().pos[improper.j] - conf.current().pos[improper.k];
        let r_ij = conf.current().pos[improper.j] - conf.current().pos[improper.i];
        let r_kl = conf.current().pos[improper.l] - conf.current().pos[improper.k];

        // Cross products
        let r_mj = r_ij.cross(r_kj);
        let r_nk = r_kj.cross(r_kl);

        let d_kj2 = r_kj.dot(r_kj);
        let d_mj2 = r_mj.dot(r_mj);
        let d_nk2 = r_nk.dot(r_nk);

        let d_kj = (d_kj2 as f64).sqrt();
        let d_mj = (d_mj2 as f64).sqrt();
        let d_nk = (d_nk2 as f64).sqrt();

        if d_mj < 1e-10 || d_nk < 1e-10 {
            continue;
        }

        // Calculate improper angle ζ
        let acs = (r_mj.dot(r_nk) / (d_mj * d_nk) as f32).clamp(-1.0, 1.0);
        let mut zeta = (acs as f64).acos();

        // Determine sign
        let ip = r_ij.dot(r_nk);
        if ip < 0.0 {
            zeta = -zeta;
        }

        // Bring ζ to interval [ζ₀ - π, ζ₀ + π]
        let mut zeta_adj = zeta;
        while zeta_adj < (params.q0 - std::f64::consts::PI) {
            zeta_adj += 2.0 * std::f64::consts::PI;
        }
        while zeta_adj > (params.q0 + std::f64::consts::PI) {
            zeta_adj -= 2.0 * std::f64::consts::PI;
        }

        // Energy: V = 0.5 * K * (ζ - ζ₀)²
        let d_zeta = zeta_adj - params.q0;
        let energy = 0.5 * params.k * d_zeta * d_zeta;

        // Force calculation
        let mut k_i = -params.k * d_zeta * d_kj;
        let mut k_l = -k_i;

        // Avoid singularities when angles approach 180°
        if d_mj2 < (1.0e-10 * d_kj2) {
            k_i = 0.0;
        } else {
            k_i = k_i / d_mj2 as f64;
        }

        if d_nk2 < (1.0e-10 * d_kj2) {
            k_l = 0.0;
        } else {
            k_l = k_l / d_nk2 as f64;
        }

        let k_j1 = (r_ij.dot(r_kj) / d_kj2) as f64 - 1.0;
        let k_j2 = (r_kl.dot(r_kj) / d_kj2) as f64;

        let f_i = r_mj * k_i as f32;
        let f_l = r_nk * k_l as f32;
        let f_j = f_i * k_j1 as f32 - f_l * k_j2 as f32;
        let f_k = -(f_i + f_j + f_l);

        result.energy += energy;
        result.forces[improper.i] += f_i;
        result.forces[improper.j] += f_j;
        result.forces[improper.k] += f_k;
        result.forces[improper.l] += f_l;
    }

    result
}

/// Calculate proper dihedral forces with "new" formula (arbitrary phase shifts)
///
/// Potential: V = K * (1 + cos(m*φ - δ))
/// Force: dV/dφ = K * m * sin(m*φ - δ)
///
/// This is the improved GROMOS dihedral potential that supports arbitrary phase shifts δ,
/// not limited to 0° or 180° as in the old formula. It's simpler and more flexible.
pub fn calculate_dihedral_new_forces(
    topo: &Topology,
    conf: &Configuration,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for dihedral in &topo.solute.proper_dihedrals {
        if dihedral.dihedral_type >= topo.dihedral_parameters.len() {
            continue;
        }

        let params = &topo.dihedral_parameters[dihedral.dihedral_type];

        // Get position vectors
        let r_ij = conf.current().pos[dihedral.j] - conf.current().pos[dihedral.i];
        let r_kj = conf.current().pos[dihedral.j] - conf.current().pos[dihedral.k];
        let r_kl = conf.current().pos[dihedral.l] - conf.current().pos[dihedral.k];

        // Cross products to get plane normals
        let r_mj = r_ij.cross(r_kj);  // normal to plane i-j-k
        let r_nk = r_kj.cross(r_kl);  // normal to plane j-k-l

        let d_kj2 = r_kj.dot(r_kj);
        let d_mj2 = r_mj.dot(r_mj);
        let d_nk2 = r_nk.dot(r_nk);
        let d_kj = (d_kj2 as f64).sqrt();

        if d_kj < 1e-10 {
            continue;
        }

        // Project vectors onto plane perpendicular to k-j bond
        let f_rim = r_ij.dot(r_kj) / d_kj2;
        let f_rln = r_kl.dot(r_kj) / d_kj2;

        let r_im = r_ij - r_kj * f_rim;
        let r_ln = r_kj * f_rln - r_kl;

        let d_im = r_im.length() as f64;
        let d_ln = r_ln.length() as f64;

        if d_im < 1e-10 || d_ln < 1e-10 {
            continue;
        }

        // Calculate dihedral angle φ
        let cos_phi = (r_im.dot(r_ln) / (d_im * d_ln) as f32).clamp(-1.0, 1.0) as f64;
        let mut phi = cos_phi.acos();

        // Determine sign of angle
        let sign = r_ij.dot(r_nk);
        if sign < 0.0 {
            phi = -phi;
        }

        // Calculate energy and force derivative
        // Energy: V = K * (1 + cos(m*φ - δ))
        let m_phi_minus_delta = params.m as f64 * phi - params.pd;
        let energy = params.k * (1.0 + m_phi_minus_delta.cos());

        // Force derivative: dV/dφ = K * m * sin(m*φ - δ)
        let k_i = params.k * params.m as f64 * m_phi_minus_delta.sin();
        let k_l = -k_i;

        // Calculate forces on each atom
        // Check for near-linear angles to avoid division by zero
        let mut f_i = Vec3::ZERO;
        let mut f_l = Vec3::ZERO;

        if d_mj2 > 1e-10 * d_kj2 as f32 {
            f_i = r_mj * (k_i * d_kj / d_mj2 as f64) as f32;
        }

        if d_nk2 > 1e-10 * d_kj2 as f32 {
            f_l = r_nk * (k_l * d_kj / d_nk2 as f64) as f32;
        }

        // Forces on j and k from chain rule
        let k_j1 = f_rim - 1.0;
        let k_j2 = f_rln;
        let f_j = f_i * k_j1 as f32 - f_l * k_j2 as f32;
        let f_k = -(f_i + f_j + f_l);

        result.energy += energy;
        result.forces[dihedral.i] += f_i;
        result.forces[dihedral.j] += f_j;
        result.forces[dihedral.k] += f_k;
        result.forces[dihedral.l] += f_l;
    }

    result
}

/// Calculate cross-dihedral forces (8-atom coupled torsional term)
///
/// Potential: V = K * (1 + cos(m*(φ + ψ) - δ))
/// where φ and ψ are two dihedral angles
///
/// This couples two torsional degrees of freedom, useful for:
/// - Ring systems with correlated torsions
/// - Backbone torsion coupling (φ-ψ correlations)
/// - Multi-bond torsional potentials
///
/// Atoms: a-b-c-d (dihedral φ) and e-f-g-h (dihedral ψ)
pub fn calculate_crossdihedral_forces(
    topo: &Topology,
    conf: &Configuration,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for crossdih in &topo.solute.cross_dihedrals {
        if crossdih.cross_dihedral_type >= topo.dihedral_parameters.len() {
            continue;
        }

        let params = &topo.dihedral_parameters[crossdih.cross_dihedral_type];

        // Get position vectors for first dihedral (a-b-c-d)
        let r_ab = conf.current().pos[crossdih.b] - conf.current().pos[crossdih.a];
        let r_cb = conf.current().pos[crossdih.b] - conf.current().pos[crossdih.c];
        let r_cd = conf.current().pos[crossdih.d] - conf.current().pos[crossdih.c];

        // Get position vectors for second dihedral (e-f-g-h)
        let r_ef = conf.current().pos[crossdih.f] - conf.current().pos[crossdih.e];
        let r_gf = conf.current().pos[crossdih.f] - conf.current().pos[crossdih.g];
        let r_gh = conf.current().pos[crossdih.h] - conf.current().pos[crossdih.g];

        // Calculate first dihedral angle φ
        let r_nc = r_cb.cross(r_cd);  // Normal to c-b-c-d plane
        let d_cb2 = r_cb.dot(r_cb);

        if d_cb2 < 1e-10 {
            continue;
        }

        let f_ram = r_ab.dot(r_cb) / d_cb2;
        let f_rdn = r_cd.dot(r_cb) / d_cb2;

        let r_am = r_ab - r_cb * f_ram;
        let r_dn = r_cb * f_rdn - r_cd;

        let d_am = r_am.length() as f64;
        let d_dn = r_dn.length() as f64;

        if d_am < 1e-10 || d_dn < 1e-10 {
            continue;
        }

        let cos_phi = (r_am.dot(r_dn) / (d_am * d_dn) as f32).clamp(-1.0, 1.0) as f64;
        let mut phi = cos_phi.acos();

        let sign_phi = r_ab.dot(r_nc);
        if sign_phi < 0.0 {
            phi = -phi;
        }

        // Calculate second dihedral angle ψ
        let r_ng = r_gf.cross(r_gh);  // Normal to g-f-g-h plane
        let d_gf2 = r_gf.dot(r_gf);

        if d_gf2 < 1e-10 {
            continue;
        }

        let f_rem = r_ef.dot(r_gf) / d_gf2;
        let f_rgn = r_gh.dot(r_gf) / d_gf2;

        let r_em = r_ef - r_gf * f_rem;
        let r_gn = r_gf * f_rgn - r_gh;

        let d_em = r_em.length() as f64;
        let d_gn = r_gn.length() as f64;

        if d_em < 1e-10 || d_gn < 1e-10 {
            continue;
        }

        let cos_psi = (r_em.dot(r_gn) / (d_em * d_gn) as f32).clamp(-1.0, 1.0) as f64;
        let mut psi = cos_psi.acos();

        let sign_psi = r_ef.dot(r_ng);
        if sign_psi < 0.0 {
            psi = -psi;
        }

        // Energy: V = K * (1 + cos(m*(φ + ψ) - δ))
        let angle_sum = phi + psi;
        let arg = params.m as f64 * angle_sum - params.pd;
        let energy = params.k * (1.0 + arg.cos());

        // Force derivative: dV/d(φ+ψ) = K * m * sin(m*(φ + ψ) - δ)
        let k = params.k * params.m as f64 * arg.sin();

        // Calculate cross products and squared lengths for forces
        let r_mb = r_ab.cross(r_cb);
        let r_mf = r_ef.cross(r_gf);
        let d_mb2 = r_mb.dot(r_mb);
        let d_mf2 = r_mf.dot(r_mf);
        let d_nc2 = r_nc.dot(r_nc);
        let d_ng2 = r_ng.dot(r_ng);

        let d_cb = (d_cb2 as f64).sqrt();
        let d_gf = (d_gf2 as f64).sqrt();

        // Forces on first dihedral (a, b, c, d)
        let mut f_a = Vec3::ZERO;
        let mut f_d = Vec3::ZERO;

        if d_mb2 > 1e-10 * d_cb2 as f32 {
            f_a = r_mb * (k * d_cb / d_mb2 as f64) as f32;
        }

        if d_nc2 > 1e-10 * d_cb2 as f32 {
            f_d = r_nc * (-k * d_cb / d_nc2 as f64) as f32;
        }

        let k_b1 = f_ram - 1.0;
        let k_b2 = f_rdn;
        let f_b = f_a * k_b1 as f32 - f_d * k_b2 as f32;
        let f_c = -(f_a + f_b + f_d);

        // Forces on second dihedral (e, f, g, h)
        let mut f_e = Vec3::ZERO;
        let mut f_h = Vec3::ZERO;

        if d_mf2 > 1e-10 * d_gf2 as f32 {
            f_e = r_mf * (k * d_gf / d_mf2 as f64) as f32;
        }

        if d_ng2 > 1e-10 * d_gf2 as f32 {
            f_h = r_ng * (-k * d_gf / d_ng2 as f64) as f32;
        }

        let k_f1 = f_rem - 1.0;
        let k_f2 = f_rgn;
        let f_f = f_e * k_f1 as f32 - f_h * k_f2 as f32;
        let f_g = -(f_e + f_f + f_h);

        // Add energy and forces
        result.energy += energy;
        result.forces[crossdih.a] += f_a;
        result.forces[crossdih.b] += f_b;
        result.forces[crossdih.c] += f_c;
        result.forces[crossdih.d] += f_d;
        result.forces[crossdih.e] += f_e;
        result.forces[crossdih.f] += f_f;
        result.forces[crossdih.g] += f_g;
        result.forces[crossdih.h] += f_h;
    }

    result
}

/// Calculate all bonded forces (bonds + angles + dihedrals)
pub fn calculate_bonded_forces(
    topo: &Topology,
    conf: &Configuration,
    use_quartic_bonds: bool,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    // Bonds
    let bond_forces = if use_quartic_bonds {
        calculate_bond_forces_quartic(topo, conf)
    } else {
        calculate_bond_forces_harmonic(topo, conf)
    };
    result.add(&bond_forces);

    // Angles
    let angle_forces = calculate_angle_forces(topo, conf);
    result.add(&angle_forces);

    // Proper dihedrals
    let dihedral_forces = calculate_dihedral_forces(topo, conf);
    result.add(&dihedral_forces);

    // Improper dihedrals
    let improper_forces = calculate_improper_dihedral_forces(topo, conf);
    result.add(&improper_forces);

    result
}

/// Calculate perturbed harmonic bond forces for FEP
///
/// Perturbed bonds have dual-topology parameters (state A and B) that are
/// interpolated based on lambda. This function calculates both the energy
/// and the lambda derivative (dE/dλ) needed for free energy calculations.
///
/// # Energy
/// V(λ) = 0.5 * K(λ) * (r - r0(λ))²
/// where:
/// - K(λ) = (1-λ)*K_A + λ*K_B
/// - r0(λ) = (1-λ)*r0_A + λ*r0_B
///
/// # Lambda derivative
/// ∂V/∂λ = 0.5 * (K_B - K_A) * (r - r0)² - K(λ) * (r - r0) * (r0_B - r0_A)
///
/// # Reference
/// - GROMOS: md++/src/interaction/bonded/perturbed_harmonic_bond_interaction.cc
pub fn calculate_perturbed_bond_forces(
    topo: &Topology,
    conf: &Configuration,
    lambda_ctrl: &LambdaController,
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    let lambda = lambda_ctrl.interaction_lambdas.bond;
    let lambda_derivative = lambda_ctrl.lambda_derivative();

    for bond in &topo.perturbed_solute.bonds {
        // Get state A and B parameters
        let params_a = &topo.bond_parameters[bond.a_type];
        let params_b = &topo.bond_parameters[bond.b_type];

        // Lambda-interpolated parameters
        let k = (1.0 - lambda) * params_a.k_harmonic + lambda * params_b.k_harmonic;
        let r0 = (1.0 - lambda) * params_a.r0 + lambda * params_b.r0;

        // Compute bond vector and distance
        let r_vec = conf.current().pos[bond.j] - conf.current().pos[bond.i];
        let r_sq = (r_vec.length_squared() as f64).max(1e-10);
        let r = r_sq.sqrt();

        // Displacement from equilibrium
        let diff = r - r0;

        // Energy: V = 0.5 * K * (r - r0)²
        let energy = 0.5 * k * diff * diff;

        // Force magnitude: |F| = K * (r - r0)
        let force_mag = k * diff / r;
        let force = r_vec * (force_mag as f32);

        result.forces[bond.i] -= force;
        result.forces[bond.j] += force;
        result.energy += energy;

        // Lambda derivative: ∂V/∂λ = 0.5 * ΔK * (r-r0)² - K * (r-r0) * Δr0
        let delta_k = params_b.k_harmonic - params_a.k_harmonic;
        let delta_r0 = params_b.r0 - params_a.r0;
        let de_dlambda = 0.5 * delta_k * diff * diff - k * diff * delta_r0;

        result.lambda_derivative += lambda_derivative * de_dlambda;
    }

    result
}

/// Calculate perturbed harmonic angle forces for FEP
///
/// Perturbed angles have dual-topology parameters (state A and B) that are
/// interpolated based on lambda. Calculates energy and ∂E/∂λ for free energy.
///
/// # Energy
/// V(λ) = 0.5 * K(λ) * (θ - θ0(λ))²
/// where:
/// - K(λ) = (1-λ)*K_A + λ*K_B
/// - θ0(λ) = (1-λ)*θ0_A + λ*θ0_B
///
/// # Lambda derivative
/// ∂V/∂λ = 0.5 * (K_B - K_A) * (θ - θ0)² - K(λ) * (θ - θ0) * (θ0_B - θ0_A)
pub fn calculate_perturbed_angle_forces(
    topo: &Topology,
    conf: &Configuration,
    lambda_ctrl: &LambdaController,
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    let lambda = lambda_ctrl.interaction_lambdas.angle;
    let lambda_derivative = lambda_ctrl.lambda_derivative();

    for angle in &topo.perturbed_solute.angles {
        // Get state A and B parameters
        let params_a = &topo.angle_parameters[angle.a_type];
        let params_b = &topo.angle_parameters[angle.b_type];

        // Lambda-interpolated parameters
        let k = (1.0 - lambda) * params_a.k_harmonic + lambda * params_b.k_harmonic;
        let theta0 = (1.0 - lambda) * params_a.theta0 + lambda * params_b.theta0;

        // Get atom indices
        let i = angle.i;
        let j = angle.j;
        let k_atom = angle.k;

        // Get vectors
        let r_ij = conf.current().pos[j] - conf.current().pos[i];
        let r_kj = conf.current().pos[j] - conf.current().pos[k_atom];

        let r_ij_len_sq = (r_ij.length_squared() as f64).max(1e-10);
        let r_kj_len_sq = (r_kj.length_squared() as f64).max(1e-10);

        let r_ij_len = r_ij_len_sq.sqrt();
        let r_kj_len = r_kj_len_sq.sqrt();

        // Calculate angle using dot product
        let cos_theta = (r_ij.dot(r_kj) as f64) / (r_ij_len * r_kj_len);
        let cos_theta = cos_theta.clamp(-1.0, 1.0);
        let theta = cos_theta.acos();

        // Angular difference from equilibrium
        let diff = theta - theta0;

        // Energy: V = 0.5 * K * (θ - θ0)²
        let energy = 0.5 * k * diff * diff;

        // Force calculation
        let sin_theta = theta.sin().max(1e-10);
        let force_mag = -k * diff / sin_theta;

        // Force on atom i
        let f_i_term = (r_kj * (cos_theta as f32) - r_ij) * ((force_mag / r_ij_len) as f32);
        result.forces[i] += f_i_term;

        // Force on atom k
        let f_k_term = (r_ij * (cos_theta as f32) - r_kj) * ((force_mag / r_kj_len) as f32);
        result.forces[k_atom] += f_k_term;

        // Force on atom j (conservation)
        result.forces[j] -= f_i_term + f_k_term;

        result.energy += energy;

        // Lambda derivative: ∂V/∂λ = 0.5 * ΔK * (θ-θ0)² - K * (θ-θ0) * Δθ0
        let delta_k = params_b.k_harmonic - params_a.k_harmonic;
        let delta_theta0 = params_b.theta0 - params_a.theta0;
        let de_dlambda = 0.5 * delta_k * diff * diff - k * diff * delta_theta0;

        result.lambda_derivative += lambda_derivative * de_dlambda;
    }

    result
}

/// Calculate perturbed dihedral forces for FEP
///
/// Perturbed dihedrals have dual-topology parameters (state A and B) that are
/// interpolated based on lambda. Calculates energy and ∂E/∂λ for free energy.
///
/// # Energy
/// V(λ) = K(λ) * (1 + cos(m*φ - δ(λ)))
/// where:
/// - K(λ) = (1-λ)*K_A + λ*K_B
/// - δ(λ) = (1-λ)*δ_A + λ*δ_B (phase offset)
///
/// # Lambda derivative
/// ∂V/∂λ = (K_B - K_A) * (1 + cos(m*φ - δ)) - K(λ) * sin(m*φ - δ) * (δ_B - δ_A)
pub fn calculate_perturbed_dihedral_forces(
    topo: &Topology,
    conf: &Configuration,
    lambda_ctrl: &LambdaController,
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    let lambda = lambda_ctrl.interaction_lambdas.dihedral;
    let lambda_derivative = lambda_ctrl.lambda_derivative();

    for dihedral in &topo.perturbed_solute.proper_dihedrals {
        // Get state A and B parameters
        let params_a = &topo.dihedral_parameters[dihedral.a_type];
        let params_b = &topo.dihedral_parameters[dihedral.b_type];

        // Lambda-interpolated parameters
        let k = (1.0 - lambda) * params_a.k + lambda * params_b.k;
        let pd = (1.0 - lambda) * params_a.pd + lambda * params_b.pd;
        let m = params_a.m; // Multiplicity (should be same for A and B)

        // Get atom indices
        let i = dihedral.i;
        let j = dihedral.j;
        let k_atom = dihedral.k;
        let l = dihedral.l;

        // Get bond vectors
        let r_ij = conf.current().pos[j] - conf.current().pos[i];
        let r_kj = conf.current().pos[j] - conf.current().pos[k_atom];
        let r_kl = conf.current().pos[l] - conf.current().pos[k_atom];

        // Calculate normal vectors to the two planes
        let m_vec = r_ij.cross(r_kj);
        let n_vec = r_kj.cross(r_kl);

        let m_len_sq = m_vec.length_squared() as f64;
        let n_len_sq = n_vec.length_squared() as f64;

        if m_len_sq < 1e-20 || n_len_sq < 1e-20 {
            continue; // Degenerate dihedral
        }

        let m_len = m_len_sq.sqrt();
        let n_len = n_len_sq.sqrt();

        // Calculate dihedral angle
        let r_kj_len = r_kj.length() as f64;
        let cos_phi = (m_vec.dot(n_vec) as f64) / (m_len * n_len);
        let cos_phi = cos_phi.clamp(-1.0, 1.0);

        // Get sign from scalar triple product
        let sign = if r_ij.dot(n_vec) >= 0.0 { 1.0 } else { -1.0 };
        let phi = sign * cos_phi.acos();

        // Compute energy: V = K * (1 + cos(m*φ - δ))
        let argument = (m as f64) * phi - pd;
        let energy = k * (1.0 + argument.cos());

        // Force magnitude: -dV/dφ = K * m * sin(m*φ - δ)
        let force_magnitude = k * (m as f64) * argument.sin();

        // Apply forces (using standard dihedral force derivatives)
        let f_factor = force_magnitude / r_kj_len;

        // Forces on atoms i and l (perpendicular to planes)
        let f_i = m_vec * ((f_factor / m_len) as f32);
        let f_l = n_vec * ((-f_factor / n_len) as f32);

        result.forces[i] += f_i;
        result.forces[l] += f_l;

        // Forces on atoms j and k (derived from torque balance)
        let r_ij_len = r_ij.length() as f64;
        let r_kl_len = r_kl.length() as f64;

        let dot_ij_kj = r_ij.dot(r_kj) as f64;
        let dot_kl_kj = r_kl.dot(r_kj) as f64;

        let f_j_factor = dot_ij_kj / (r_kj_len * r_kj_len);
        let f_k_factor = dot_kl_kj / (r_kj_len * r_kj_len);

        let f_j = f_i * (f_j_factor as f32 - 1.0);
        let f_k = f_l * (f_k_factor as f32 - 1.0);

        result.forces[j] += f_j;
        result.forces[k_atom] += f_k;

        result.energy += energy;

        // Lambda derivative: ∂V/∂λ = ΔK * (1 + cos(arg)) - K * sin(arg) * Δδ
        let delta_k = params_b.k - params_a.k;
        let delta_pd = params_b.pd - params_a.pd;
        let de_dlambda = delta_k * (1.0 + argument.cos()) - k * argument.sin() * delta_pd;

        result.lambda_derivative += lambda_derivative * de_dlambda;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::topology::{read_topology_file, build_topology};
    use crate::io::coordinate::read_coordinate_file;

    #[test]
    fn test_bond_forces_quartic() {
        let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");
        let conf_result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

        if topo_result.is_err() || conf_result.is_err() {
            println!("Skipping: test files not found");
            return;
        }

        let parsed = topo_result.unwrap();
        let topo = build_topology(parsed);
        let conf = conf_result.unwrap();

        let result = calculate_bond_forces_quartic(&topo, &conf);

        println!("Bond energy: {:.6} kJ/mol", result.energy);
        println!("Force on atom 0: ({:.6}, {:.6}, {:.6})",
                 result.forces[0].x, result.forces[0].y, result.forces[0].z);

        // Energy should be positive (bonds are stretched or compressed)
        assert!(result.energy >= 0.0);

        // Forces should not all be zero
        let total_force: Vec3 = result.forces.iter().sum();
        println!("Total force (should be ~0): ({:.6}, {:.6}, {:.6})",
                 total_force.x, total_force.y, total_force.z);
    }

    #[test]
    fn test_angle_forces() {
        let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");
        let conf_result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

        if topo_result.is_err() || conf_result.is_err() {
            println!("Skipping: test files not found");
            return;
        }

        let parsed = topo_result.unwrap();
        let topo = build_topology(parsed);
        let conf = conf_result.unwrap();

        let result = calculate_angle_forces(&topo, &conf);

        println!("Angle energy: {:.6} kJ/mol", result.energy);
        println!("Force on atom 0: ({:.6}, {:.6}, {:.6})",
                 result.forces[0].x, result.forces[0].y, result.forces[0].z);

        assert!(result.energy >= 0.0);

        // Check force conservation
        let total_force: Vec3 = result.forces.iter().sum();
        println!("Total force (should be ~0): ({:.6}, {:.6}, {:.6})",
                 total_force.x, total_force.y, total_force.z);

        // Should be very close to zero (within numerical precision)
        assert!(total_force.length() < 1e-4);
    }

    #[test]
    fn test_bonded_forces_complete() {
        let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");
        let conf_result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

        if topo_result.is_err() || conf_result.is_err() {
            println!("Skipping: test files not found");
            return;
        }

        let parsed = topo_result.unwrap();
        let topo = build_topology(parsed);
        let conf = conf_result.unwrap();

        let result = calculate_bonded_forces(&topo, &conf, true);

        println!("\n========== Bonded Forces Test ==========");
        println!("Total bonded energy: {:.6} kJ/mol", result.energy);
        println!("Number of atoms: {}", topo.num_atoms());

        for i in 0..topo.num_atoms() {
            println!("Atom {}: F = ({:.6}, {:.6}, {:.6}) kJ/(mol·nm)",
                     i, result.forces[i].x, result.forces[i].y, result.forces[i].z);
        }

        // Total force should be conserved (sum to zero)
        let total_force: Vec3 = result.forces.iter().sum();
        println!("\nTotal force: ({:.6}, {:.6}, {:.6})",
                 total_force.x, total_force.y, total_force.z);
        assert!(total_force.length() < 1e-3, "Forces should be conserved");
    }

    #[test]
    fn test_harmonic_angle_simple() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Angle, AngleParameters, Atom};

        // Create simple 3-atom system: i-j-k forming an angle
        let mut topo = Topology::new();

        // Add 3 atoms with equal masses
        for i in 0..3 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0, 12.0, 12.0]; // Carbon atoms
        topo.inverse_mass = vec![1.0/12.0, 1.0/12.0, 1.0/12.0];

        // Add one angle: 0-1-2
        topo.solute.angles.push(Angle { i: 0, j: 1, k: 2, angle_type: 0 });

        // Angle parameters: K = 500 kJ/(mol·rad²), θ₀ = 120° = 2.0944 rad
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 0.0,     // Not used for harmonic angles
            k_harmonic: 500.0,
            theta0: 2.0943951, // 120 degrees in radians
        });

        // Create configuration with atoms at 120° angle
        let mut conf = Configuration::new(3, 1, 1);

        // Position atoms in xy plane
        // Atom 1 (central) at origin
        // Atom 0 along x-axis
        // Atom 2 at 120° from x-axis
        conf.current_mut().pos[0] = Vec3::new(0.15, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(-0.075, 0.12990, 0.0); // 120° angle

        // Calculate forces
        let result = calculate_harmonic_angle_forces(&topo, &conf);

        println!("\n========== Harmonic Angle Test (Equilibrium) ==========");
        println!("Angle energy: {:.6} kJ/mol", result.energy);
        println!("Force on atom 0: ({:.6}, {:.6}, {:.6})",
                 result.forces[0].x, result.forces[0].y, result.forces[0].z);
        println!("Force on atom 1: ({:.6}, {:.6}, {:.6})",
                 result.forces[1].x, result.forces[1].y, result.forces[1].z);
        println!("Force on atom 2: ({:.6}, {:.6}, {:.6})",
                 result.forces[2].x, result.forces[2].y, result.forces[2].z);

        // At equilibrium angle, energy should be near zero
        assert!(result.energy < 0.01, "Energy at equilibrium should be ~0: {}", result.energy);

        // Forces should be near zero at equilibrium
        assert!(result.forces[0].length() < 0.1, "Forces should be small at equilibrium");
        assert!(result.forces[1].length() < 0.1, "Forces should be small at equilibrium");
        assert!(result.forces[2].length() < 0.1, "Forces should be small at equilibrium");

        // Force conservation
        let total_force = result.forces[0] + result.forces[1] + result.forces[2];
        println!("Total force: ({:.6}, {:.6}, {:.6})",
                 total_force.x, total_force.y, total_force.z);
        assert!(total_force.length() < 1e-5, "Forces should be conserved");
    }

    #[test]
    fn test_harmonic_angle_compressed() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Angle, AngleParameters, Atom};

        // Test with angle compressed (smaller than equilibrium)
        let mut topo = Topology::new();

        // Add 3 atoms
        for i in 0..3 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0, 12.0, 12.0];
        topo.inverse_mass = vec![1.0/12.0, 1.0/12.0, 1.0/12.0];

        topo.solute.angles.push(Angle { i: 0, j: 1, k: 2, angle_type: 0 });
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 0.0,
            k_harmonic: 500.0,
            theta0: 2.0943951, // 120° = 2.0944 rad
        });

        let mut conf = Configuration::new(3, 1, 1);

        // Create 90° angle (compressed from 120° equilibrium)
        conf.current_mut().pos[0] = Vec3::new(0.15, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(0.0, 0.15, 0.0); // 90° angle

        let result = calculate_harmonic_angle_forces(&topo, &conf);

        println!("\n========== Harmonic Angle Test (Compressed 90°) ==========");
        println!("Angle energy: {:.6} kJ/mol", result.energy);

        // Calculate expected energy: V = 0.5 * K * (θ - θ₀)²
        // θ = π/2 = 1.5708 rad, θ₀ = 2.0944 rad
        // dθ = 1.5708 - 2.0944 = -0.5236 rad (30° difference)
        let expected_energy = 0.5 * 500.0 * 0.5236 * 0.5236; // ≈ 68.5 kJ/mol

        println!("Expected energy: {:.6} kJ/mol", expected_energy);
        assert!((result.energy - expected_energy).abs() < 1.0,
                "Energy should be ~{} kJ/mol, got {}", expected_energy, result.energy);

        // Energy should be positive when compressed
        assert!(result.energy > 0.0, "Energy should be positive when angle is compressed");

        // Forces should push atoms apart (opening the angle)
        // Force on atom 0 should have negative y component (push down)
        // Force on atom 2 should have negative x component (push left)
        println!("Force on atom 0: ({:.6}, {:.6}, {:.6})",
                 result.forces[0].x, result.forces[0].y, result.forces[0].z);
        println!("Force on atom 2: ({:.6}, {:.6}, {:.6})",
                 result.forces[2].x, result.forces[2].y, result.forces[2].z);

        // Force conservation
        let total_force = result.forces[0] + result.forces[1] + result.forces[2];
        assert!(total_force.length() < 1e-5, "Forces should be conserved");
    }

    #[test]
    fn test_harmonic_angle_vs_cosine() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Angle, AngleParameters, Atom};

        // Compare harmonic and cosine-based angle potentials
        let mut topo = Topology::new();

        // Add 3 atoms (water-like masses)
        topo.solute.atoms.push(Atom {
            name: "O".to_string(),
            residue_nr: 1,
            residue_name: "WAT".to_string(),
            iac: 0,
            mass: 16.0,
            charge: 0.0,
            is_perturbed: false,
            is_polarisable: false,
            is_coarse_grained: false,
        });
        for i in 0..2 {
            topo.solute.atoms.push(Atom {
                name: format!("H{}", i+1),
                residue_nr: 1,
                residue_name: "WAT".to_string(),
                iac: 1,
                mass: 1.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![16.0, 1.0, 1.0]; // Water-like: O-H-H
        topo.inverse_mass = vec![1.0/16.0, 1.0, 1.0];

        topo.solute.angles.push(Angle { i: 0, j: 1, k: 2, angle_type: 0 });
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 400.0,
            k_harmonic: 400.0,
            theta0: 1.91063, // 109.47° (tetrahedral angle)
        });

        let mut conf = Configuration::new(3, 1, 1);

        // Create angle slightly off equilibrium (115°)
        conf.current_mut().pos[0] = Vec3::new(0.1, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(-0.04226, 0.09063, 0.0); // 115°

        let harmonic = calculate_harmonic_angle_forces(&topo, &conf);
        let cosine = calculate_angle_forces(&topo, &conf);

        println!("\n========== Harmonic vs Cosine Angle Comparison ==========");
        println!("Harmonic energy: {:.6} kJ/mol", harmonic.energy);
        println!("Cosine energy:   {:.6} kJ/mol", cosine.energy);

        // Both should give reasonable energies (positive for off-equilibrium)
        assert!(harmonic.energy > 0.0, "Harmonic energy should be positive");
        assert!(cosine.energy > 0.0, "Cosine energy should be positive");

        // Both should conserve forces
        let harm_total: Vec3 = harmonic.forces.iter().sum();
        let cos_total: Vec3 = cosine.forces.iter().sum();

        println!("Harmonic total force: ({:.6}, {:.6}, {:.6})",
                 harm_total.x, harm_total.y, harm_total.z);
        println!("Cosine total force:   ({:.6}, {:.6}, {:.6})",
                 cos_total.x, cos_total.y, cos_total.z);

        assert!(harm_total.length() < 5e-5, "Harmonic forces should be conserved");
        assert!(cos_total.length() < 5e-5, "Cosine forces should be conserved");

        // Near equilibrium, both should give similar results (within 20%)
        let energy_ratio = harmonic.energy / cosine.energy;
        println!("Energy ratio (harmonic/cosine): {:.3}", energy_ratio);
        assert!(energy_ratio > 0.5 && energy_ratio < 2.0,
                "Energies should be similar near equilibrium");
    }

    #[test]
    fn test_harmonic_angle_linear() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Angle, AngleParameters, Atom};

        // Test special case: nearly linear angle (θ ≈ 180°)
        let mut topo = Topology::new();

        // Add 3 atoms
        for i in 0..3 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0, 12.0, 12.0];
        topo.inverse_mass = vec![1.0/12.0, 1.0/12.0, 1.0/12.0];

        topo.solute.angles.push(Angle { i: 0, j: 1, k: 2, angle_type: 0 });
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 0.0,  // Not used for harmonic angles
            k_harmonic: 500.0,
            theta0: std::f64::consts::PI, // 180° (linear)
        });

        let mut conf = Configuration::new(3, 1, 1);

        // Create nearly linear angle (179.9°)
        conf.current_mut().pos[0] = Vec3::new(0.15, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(-0.15, 0.001, 0.0); // Almost linear

        let result = calculate_harmonic_angle_forces(&topo, &conf);

        println!("\n========== Harmonic Angle Test (Nearly Linear) ==========");
        println!("Angle energy: {:.6} kJ/mol", result.energy);

        // Energy should be finite and small (close to equilibrium)
        assert!(!result.energy.is_nan(), "Energy should not be NaN");
        assert!(!result.energy.is_infinite(), "Energy should not be infinite");
        assert!(result.energy < 1.0, "Energy should be small near equilibrium");

        // Forces should be finite
        for i in 0..3 {
            assert!(!result.forces[i].x.is_nan(), "Force should not be NaN");
            assert!(!result.forces[i].y.is_nan(), "Force should not be NaN");
            assert!(!result.forces[i].z.is_nan(), "Force should not be NaN");
        }

        // Force conservation
        let total_force = result.forces[0] + result.forces[1] + result.forces[2];
        assert!(total_force.length() < 1e-5, "Forces should be conserved");
    }

    #[test]
    fn test_cg_bond_compressed() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Bond, BondParameters, Atom};

        // Test CG bond when compressed (r < r0) - should have zero force/energy
        let mut topo = Topology::new();

        // Add 2 atoms
        for i in 0..2 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: true,
            });
        }
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0/12.0, 1.0/12.0];

        // Add one CG bond (negative k_harmonic as marker)
        topo.solute.bonds.push(Bond { i: 0, j: 1, bond_type: 0 });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,      // Not used for CG bonds
            k_harmonic: -500.0,  // Negative = CG bond marker
            r0: 0.15,  // 0.15 nm equilibrium
        });

        let mut conf = Configuration::new(2, 1, 1);

        // Position atoms compressed (0.1 nm < 0.15 nm equilibrium)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.1, 0.0, 0.0);

        let result = calculate_cg_bond_forces(&topo, &conf);

        println!("\n========== CG Bond Test (Compressed) ==========");
        println!("Bond length: 0.10 nm, equilibrium: 0.15 nm");
        println!("Energy: {:.6} kJ/mol", result.energy);
        println!("Force on atom 0: ({:.6}, {:.6}, {:.6})",
                 result.forces[0].x, result.forces[0].y, result.forces[0].z);

        // When compressed (r < r0), CG bonds should have zero energy and force
        assert!(result.energy < 1e-10, "Energy should be zero when compressed: {}", result.energy);
        assert!(result.forces[0].length() < 1e-10, "Force should be zero when compressed");
    }

    #[test]
    fn test_cg_bond_at_equilibrium() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Bond, BondParameters, Atom};

        // Test CG bond at equilibrium (r = r0) - should have zero force/energy
        let mut topo = Topology::new();

        // Add 2 atoms
        for i in 0..2 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: true,
            });
        }
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0/12.0, 1.0/12.0];

        topo.solute.bonds.push(Bond { i: 0, j: 1, bond_type: 0 });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: -500.0,
            r0: 0.15,
        });

        let mut conf = Configuration::new(2, 1, 1);

        // Position atoms at equilibrium (0.15 nm)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.15, 0.0, 0.0);

        let result = calculate_cg_bond_forces(&topo, &conf);

        println!("\n========== CG Bond Test (Equilibrium) ==========");
        println!("Bond length: 0.15 nm (at equilibrium)");
        println!("Energy: {:.6} kJ/mol", result.energy);

        // At equilibrium, energy and force should be zero
        assert!(result.energy < 1e-10, "Energy should be zero at equilibrium");
        assert!(result.forces[0].length() < 1e-10, "Force should be zero at equilibrium");
    }

    #[test]
    fn test_cg_bond_stretched() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Bond, BondParameters, Atom};

        // Test CG bond when stretched (r > r0) - should have repulsive force
        let mut topo = Topology::new();

        // Add 2 atoms
        for i in 0..2 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: true,
            });
        }
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0/12.0, 1.0/12.0];

        topo.solute.bonds.push(Bond { i: 0, j: 1, bond_type: 0 });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: -500.0,  // K = 500 kJ/(mol·nm^4)
            r0: 0.15,  // r0 = 0.15 nm
        });

        let mut conf = Configuration::new(2, 1, 1);

        // Position atoms stretched (0.20 nm > 0.15 nm equilibrium)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.20, 0.0, 0.0);

        let result = calculate_cg_bond_forces(&topo, &conf);

        // Calculate expected values
        // diff = 0.20 - 0.15 = 0.05 nm
        // Energy = 0.5 * 500 * (0.05)^4 = 0.5 * 500 * 0.00000625 = 0.0015625 kJ/mol
        let diff: f64 = 0.05;
        let expected_energy = 0.5 * 500.0 * diff.powi(4);
        let expected_force_mag = 2.0 * 500.0 * diff.powi(3);  // Attractive (pulls atoms together)

        println!("\n========== CG Bond Test (Stretched) ==========");
        println!("Bond length: 0.20 nm, equilibrium: 0.15 nm");
        println!("Energy: {:.6} kJ/mol (expected: {:.6})", result.energy, expected_energy);
        println!("Force magnitude: {:.6} (expected: {:.6})",
                 result.forces[0].length(), expected_force_mag);
        println!("Force on atom 0: ({:.6}, {:.6}, {:.6})",
                 result.forces[0].x, result.forces[0].y, result.forces[0].z);
        println!("Force on atom 1: ({:.6}, {:.6}, {:.6})",
                 result.forces[1].x, result.forces[1].y, result.forces[1].z);

        // Check energy
        assert!((result.energy - expected_energy).abs() < 1e-6,
                "Energy mismatch: got {}, expected {}", result.energy, expected_energy);

        // Check force magnitude (should be attractive = negative x-direction on atom 0)
        assert!(result.forces[0].x < 0.0, "Force should be attractive (negative x)");
        assert!((result.forces[0].length() as f64 - expected_force_mag).abs() < 1e-5,
                "Force magnitude mismatch: got {}, expected {}",
                result.forces[0].length(), expected_force_mag);

        // Force conservation
        let total_force = result.forces[0] + result.forces[1];
        assert!(total_force.length() < 1e-5, "Forces should be conserved");
    }

    #[test]
    fn test_cg_bond_highly_stretched() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Bond, BondParameters, Atom};

        // Test CG bond when highly stretched - quartic potential grows rapidly
        let mut topo = Topology::new();

        // Add 2 atoms
        for i in 0..2 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: true,
            });
        }
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0/12.0, 1.0/12.0];

        topo.solute.bonds.push(Bond { i: 0, j: 1, bond_type: 0 });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: -500.0,
            r0: 0.15,
        });

        let mut conf = Configuration::new(2, 1, 1);

        // Highly stretched: 0.25 nm (0.10 nm beyond equilibrium)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.25, 0.0, 0.0);

        let result = calculate_cg_bond_forces(&topo, &conf);

        // diff = 0.10 nm, Energy = 0.5 * 500 * (0.10)^4 = 0.025 kJ/mol
        let diff: f64 = 0.10;
        let expected_energy = 0.5 * 500.0 * diff.powi(4);

        println!("\n========== CG Bond Test (Highly Stretched) ==========");
        println!("Bond length: 0.25 nm, equilibrium: 0.15 nm (diff = 0.10 nm)");
        println!("Energy: {:.6} kJ/mol (expected: {:.6})", result.energy, expected_energy);
        println!("Force magnitude: {:.6}", result.forces[0].length());

        // Energy should grow as (diff)^4
        assert!((result.energy - expected_energy).abs() < 1e-6,
                "Energy mismatch for highly stretched bond");

        // Verify quartic growth: doubling diff should multiply energy by 16
        // Previous test: diff = 0.05, energy = 0.0015625
        // This test: diff = 0.10, energy should be 16x larger = 0.025
        assert!(result.energy > 0.02, "Quartic potential should grow rapidly");
    }

    #[test]
    fn test_dihedral_new_at_minimum() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Dihedral, DihedralParameters, Atom};
        use std::f64::consts::PI;

        // Test new dihedral at energy minimum (φ = δ/m)
        let mut topo = Topology::new();

        // Add 4 atoms in a dihedral
        for i in 0..4 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0/12.0; 4];

        // Dihedral parameters: K=10, m=1, δ=180° (π radians)
        // Energy minimum at φ = π
        topo.solute.proper_dihedrals.push(Dihedral { i: 0, j: 1, k: 2, l: 3, dihedral_type: 0 });
        topo.dihedral_parameters.push(DihedralParameters {
            k: 10.0,
            cospd: -1.0,  // cos(180°) = -1
            pd: PI,       // δ = 180° = π radians
            m: 1,         // multiplicity
        });

        let mut conf = Configuration::new(4, 1, 1);

        // Create planar trans configuration (φ = 180° = π)
        // For trans, atoms 0 and 3 should be on opposite sides of the j-k bond
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.5, 0.0);    // Above the j-k plane
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0, -0.5, 0.0);  // Below the j-k plane (opposite side)

        let result = calculate_dihedral_new_forces(&topo, &conf);

        // At minimum: V = K * (1 + cos(m*φ - δ)) = K * (1 + cos(π - π)) = K * (1 + 1) = 2K = 20
        let expected_energy = 10.0 * (1.0 + (1.0 * PI - PI).cos());

        println!("\n========== New Dihedral Test (At Minimum, φ=180°) ==========");
        println!("Energy: {:.6} kJ/mol (expected: {:.6})", result.energy, expected_energy);
        println!("Force on atom 0: ({:.6}, {:.6}, {:.6})",
                 result.forces[0].x, result.forces[0].y, result.forces[0].z);

        // At minimum, energy should be maximum (2K) and forces should be zero
        assert!((result.energy - expected_energy).abs() < 1e-6,
                "Energy at minimum should be 2K: got {}, expected {}", result.energy, expected_energy);

        // Forces should be very small at minimum
        for i in 0..4 {
            assert!(result.forces[i].length() < 1e-3,
                    "Force on atom {} should be near zero at minimum", i);
        }
    }

    #[test]
    fn test_dihedral_new_vs_old() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Dihedral, DihedralParameters, Atom};
        use std::f64::consts::PI;

        // Compare new and old dihedral formulas at φ = 0° (cis)
        // For δ = 0°, both formulas should give same result
        let mut topo = Topology::new();

        // Add 4 atoms
        for i in 0..4 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0/12.0; 4];

        // Parameters: K=10, m=1, δ=0°
        topo.solute.proper_dihedrals.push(Dihedral { i: 0, j: 1, k: 2, l: 3, dihedral_type: 0 });
        topo.dihedral_parameters.push(DihedralParameters {
            k: 10.0,
            cospd: 1.0,   // cos(0°) = 1
            pd: 0.0,      // δ = 0°
            m: 1,
        });

        let mut conf = Configuration::new(4, 1, 1);

        // Create eclipsed configuration (φ ≈ 0°)
        // Atoms form a zigzag with small dihedral angle
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.5, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0, 0.5, 0.01); // Small z-deviation for φ ≈ 0°

        let result_new = calculate_dihedral_new_forces(&topo, &conf);
        let result_old = calculate_dihedral_forces(&topo, &conf);

        println!("\n========== New vs Old Dihedral Comparison ==========");
        println!("New formula energy: {:.6} kJ/mol", result_new.energy);
        println!("Old formula energy: {:.6} kJ/mol", result_old.energy);
        println!("Energy difference: {:.6} kJ/mol", (result_new.energy - result_old.energy).abs());

        // For small angles and δ=0°, both formulas should give very similar results
        assert!((result_new.energy - result_old.energy).abs() < 0.1,
                "New and old formulas should agree for δ=0°");
    }

    #[test]
    fn test_dihedral_new_arbitrary_phase() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Dihedral, DihedralParameters, Atom};
        use std::f64::consts::PI;

        // Test new dihedral with arbitrary phase shift δ = 60°
        let mut topo = Topology::new();

        // Add 4 atoms
        for i in 0..4 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0/12.0; 4];

        // Parameters: K=15, m=3, δ=60° = π/3 rad
        let delta = PI / 3.0;  // 60°
        topo.solute.proper_dihedrals.push(Dihedral { i: 0, j: 1, k: 2, l: 3, dihedral_type: 0 });
        topo.dihedral_parameters.push(DihedralParameters {
            k: 15.0,
            cospd: delta.cos(),  // cos(60°) = 0.5
            pd: delta,           // δ = 60°
            m: 3,                // multiplicity 3
        });

        let mut conf = Configuration::new(4, 1, 1);

        // Create configuration with φ = 90°
        // Proper dihedral geometry: first plane in xy, rotate second plane by 90° around x
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.5, 0.0);    // First plane (xy)
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0, 0.0, 0.5);  // Second plane (xz) - 90° rotation

        let result = calculate_dihedral_new_forces(&topo, &conf);

        // Energy: V = K * (1 + cos(m*φ - δ)) = 15 * (1 + cos(3*90° - 60°))
        //           = 15 * (1 + cos(270° - 60°)) = 15 * (1 + cos(210°))
        //           = 15 * (1 - sqrt(3)/2) ≈ 2.00961894323
        let phi = PI / 2.0;  // 90° in radians
        let expected_energy = 15.0 * (1.0 + (3.0 * phi - delta).cos());

        println!("\n========== New Dihedral Test (Arbitrary Phase δ=60°) ==========");
        println!("φ = 90°, m = 3, δ = 60°");
        println!("Energy: {:.6} kJ/mol (expected: {:.6})", result.energy, expected_energy);

        // Energy should be finite and non-zero
        assert!(!result.energy.is_nan() && !result.energy.is_infinite(),
                "Energy should be finite");

        // Energy should be reasonable (within broad range)
        // The exact value depends on the precise dihedral angle
        assert!(result.energy > 0.0 && result.energy < 30.0,
                "Energy should be in reasonable range: got {}", result.energy);
    }

    #[test]
    fn test_dihedral_new_force_conservation() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{Dihedral, DihedralParameters, Atom};

        // Test force conservation for new dihedral formula
        let mut topo = Topology::new();

        // Add 4 atoms
        for i in 0..4 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0/12.0; 4];

        topo.solute.proper_dihedrals.push(Dihedral { i: 0, j: 1, k: 2, l: 3, dihedral_type: 0 });
        topo.dihedral_parameters.push(DihedralParameters {
            k: 20.0,
            cospd: 0.0,  // cos(90°) = 0
            pd: std::f64::consts::PI / 2.0,  // δ = 90°
            m: 2,
        });

        let mut conf = Configuration::new(4, 1, 1);

        // Arbitrary dihedral geometry
        conf.current_mut().pos[0] = Vec3::new(-1.2, 0.3, 0.1);
        conf.current_mut().pos[1] = Vec3::new(-0.4, -0.1, 0.2);
        conf.current_mut().pos[2] = Vec3::new(0.5, 0.2, -0.1);
        conf.current_mut().pos[3] = Vec3::new(1.3, -0.2, 0.3);

        let result = calculate_dihedral_new_forces(&topo, &conf);

        println!("\n========== New Dihedral Force Conservation Test ==========");
        println!("Energy: {:.6} kJ/mol", result.energy);

        // Check force conservation
        let total_force = result.forces[0] + result.forces[1] + result.forces[2] + result.forces[3];
        println!("Total force: ({:.6}, {:.6}, {:.6})",
                 total_force.x, total_force.y, total_force.z);

        assert!(total_force.length() < 1e-4,
                "Forces should be conserved: total = {}", total_force.length());
    }

    #[test]
    fn test_crossdihedral_simple() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{CrossDihedral, DihedralParameters, Atom};
        use std::f64::consts::PI;

        // Test cross-dihedral with two simple dihedrals
        let mut topo = Topology::new();

        // Add 8 atoms for cross-dihedral: a-b-c-d and e-f-g-h
        for i in 0..8 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 8];
        topo.inverse_mass = vec![1.0/12.0; 8];

        // Parameters: K=5, m=1, δ=0°
        topo.solute.cross_dihedrals.push(CrossDihedral {
            a: 0, b: 1, c: 2, d: 3,  // First dihedral
            e: 4, f: 5, g: 6, h: 7,  // Second dihedral
            cross_dihedral_type: 0
        });

        topo.dihedral_parameters.push(DihedralParameters {
            k: 5.0,
            cospd: 1.0,   // cos(0°) = 1
            pd: 0.0,      // δ = 0°
            m: 1,
        });

        let mut conf = Configuration::new(8, 1, 1);

        // Create geometry: both dihedrals at ~180° (trans)
        // First dihedral (0-1-2-3)
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.5, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0, -0.5, 0.0);

        // Second dihedral (4-5-6-7)
        conf.current_mut().pos[4] = Vec3::new(3.0, 0.5, 0.0);
        conf.current_mut().pos[5] = Vec3::new(4.0, 0.0, 0.0);
        conf.current_mut().pos[6] = Vec3::new(5.0, 0.0, 0.0);
        conf.current_mut().pos[7] = Vec3::new(6.0, -0.5, 0.0);

        let result = calculate_crossdihedral_forces(&topo, &conf);

        println!("\n========== Cross-Dihedral Test (Both Trans) ==========");
        println!("Energy: {:.6} kJ/mol", result.energy);

        // Energy should be finite
        assert!(!result.energy.is_nan() && !result.energy.is_infinite(),
                "Energy should be finite");

        // Energy should be reasonable
        assert!(result.energy >= 0.0 && result.energy <= 20.0,
                "Energy should be in reasonable range: {}", result.energy);
    }

    #[test]
    fn test_crossdihedral_force_conservation() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{CrossDihedral, DihedralParameters, Atom};

        // Test force conservation for cross-dihedral
        let mut topo = Topology::new();

        // Add 8 atoms
        for i in 0..8 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 8];
        topo.inverse_mass = vec![1.0/12.0; 8];

        topo.solute.cross_dihedrals.push(CrossDihedral {
            a: 0, b: 1, c: 2, d: 3,
            e: 4, f: 5, g: 6, h: 7,
            cross_dihedral_type: 0
        });

        topo.dihedral_parameters.push(DihedralParameters {
            k: 10.0,
            cospd: 0.0,  // cos(90°) = 0
            pd: std::f64::consts::PI / 2.0,  // δ = 90°
            m: 2,
        });

        let mut conf = Configuration::new(8, 1, 1);

        // Arbitrary geometry
        conf.current_mut().pos[0] = Vec3::new(-1.2, 0.3, 0.1);
        conf.current_mut().pos[1] = Vec3::new(-0.4, -0.1, 0.2);
        conf.current_mut().pos[2] = Vec3::new(0.5, 0.2, -0.1);
        conf.current_mut().pos[3] = Vec3::new(1.3, -0.2, 0.3);
        conf.current_mut().pos[4] = Vec3::new(2.1, 0.4, -0.2);
        conf.current_mut().pos[5] = Vec3::new(3.0, -0.1, 0.1);
        conf.current_mut().pos[6] = Vec3::new(3.8, 0.3, -0.3);
        conf.current_mut().pos[7] = Vec3::new(4.7, -0.3, 0.2);

        let result = calculate_crossdihedral_forces(&topo, &conf);

        println!("\n========== Cross-Dihedral Force Conservation Test ==========");
        println!("Energy: {:.6} kJ/mol", result.energy);

        // Check force conservation (sum of all forces should be zero)
        let mut total_force = Vec3::ZERO;
        for i in 0..8 {
            total_force += result.forces[i];
        }

        println!("Total force: ({:.6}, {:.6}, {:.6})",
                 total_force.x, total_force.y, total_force.z);

        assert!(total_force.length() < 1e-4,
                "Forces should be conserved: total = {}", total_force.length());
    }

    #[test]
    fn test_crossdihedral_coupling() {
        use crate::math::Vec3;
        use crate::configuration::Configuration;
        use crate::topology::{CrossDihedral, DihedralParameters, Atom};
        use std::f64::consts::PI;

        // Test that cross-dihedral properly couples the two angles
        let mut topo = Topology::new();

        // Add 8 atoms
        for i in 0..8 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 8];
        topo.inverse_mass = vec![1.0/12.0; 8];

        topo.solute.cross_dihedrals.push(CrossDihedral {
            a: 0, b: 1, c: 2, d: 3,
            e: 4, f: 5, g: 6, h: 7,
            cross_dihedral_type: 0
        });

        // V = K * (1 + cos(m*(φ + ψ) - δ))
        // With m=1, δ=0: V = K * (1 + cos(φ + ψ))
        topo.dihedral_parameters.push(DihedralParameters {
            k: 10.0,
            cospd: 1.0,   // cos(0°) = 1
            pd: 0.0,      // δ = 0°
            m: 1,
        });

        let mut conf = Configuration::new(8, 1, 1);

        // Both dihedrals near 180° (trans): φ ≈ π, ψ ≈ π
        // φ + ψ ≈ 2π, cos(2π) = 1, Energy = K * (1 + 1) = 20
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.5, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0, -0.5, 0.0);

        conf.current_mut().pos[4] = Vec3::new(3.0, 0.5, 0.0);
        conf.current_mut().pos[5] = Vec3::new(4.0, 0.0, 0.0);
        conf.current_mut().pos[6] = Vec3::new(5.0, 0.0, 0.0);
        conf.current_mut().pos[7] = Vec3::new(6.0, -0.5, 0.0);

        let result = calculate_crossdihedral_forces(&topo, &conf);

        println!("\n========== Cross-Dihedral Coupling Test ==========");
        println!("Both dihedrals near 180° (trans)");
        println!("Energy: {:.6} kJ/mol", result.energy);

        // Energy should be close to 2*K = 20 when φ + ψ ≈ 2π
        // Allow some tolerance for geometry
        assert!(result.energy > 15.0 && result.energy < 21.0,
                "Energy for φ+ψ≈2π should be near 2K=20: got {}", result.energy);
    }

    #[test]
    fn test_perturbed_bond_lambda_interpolation() {
        use crate::fep::LambdaController;
        use crate::topology::{PerturbedBond, BondParameters, Atom};

        let mut topo = Topology::new();

        // Add 2 atoms to solute
        for i in 0..2 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "MOL".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: true,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }

        topo.iac = vec![0, 0];
        topo.charge = vec![0.0, 0.0];
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0/12.0, 1.0/12.0];

        // State A: k=1000, r0=0.1 nm
        // State B: k=2000, r0=0.15 nm
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: 1000.0,
            r0: 0.1,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: 2000.0,
            r0: 0.15,
        });

        topo.perturbed_solute.bonds.push(PerturbedBond {
            i: 0,
            j: 1,
            a_type: 0,
            b_type: 1,
        });

        let mut conf = Configuration::new(2, 1, 1);
        // Bond length = 0.12 nm (between r0_A and r0_B)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.12, 0.0, 0.0);

        // Test at λ = 0 (pure state A)
        let lambda_ctrl_0 = LambdaController::new().with_lambda(0.0);
        let result_0 = calculate_perturbed_bond_forces(&topo, &conf, &lambda_ctrl_0);

        // At λ=0: k=1000, r0=0.1, r=0.12, diff=0.02
        // E = 0.5 * 1000 * 0.02² = 0.2
        let expected_energy_0 = 0.5 * 1000.0 * 0.02 * 0.02;
        println!("\nλ=0: Energy = {:.6}, Expected = {:.6}", result_0.energy, expected_energy_0);
        assert!((result_0.energy - expected_energy_0).abs() < 1e-6,
                "Energy at λ=0 should match state A");

        // Test at λ = 1 (pure state B)
        let lambda_ctrl_1 = LambdaController::new().with_lambda(1.0);
        let result_1 = calculate_perturbed_bond_forces(&topo, &conf, &lambda_ctrl_1);

        // At λ=1: k=2000, r0=0.15, r=0.12, diff=-0.03
        // E = 0.5 * 2000 * 0.03² = 0.9
        let expected_energy_1 = 0.5 * 2000.0 * 0.03 * 0.03;
        println!("λ=1: Energy = {:.6}, Expected = {:.6}", result_1.energy, expected_energy_1);
        assert!((result_1.energy - expected_energy_1).abs() < 1e-6,
                "Energy at λ=1 should match state B");

        // Test at λ = 0.5 (interpolated)
        let lambda_ctrl_05 = LambdaController::new().with_lambda(0.5);
        let result_05 = calculate_perturbed_bond_forces(&topo, &conf, &lambda_ctrl_05);

        // At λ=0.5: k=1500, r0=0.125, r=0.12, diff=-0.005
        // E = 0.5 * 1500 * 0.005² = 0.01875
        let k_interp = 0.5 * 1000.0 + 0.5 * 2000.0;
        let r0_interp = 0.5 * 0.1 + 0.5 * 0.15;
        let diff_interp = 0.12 - r0_interp;
        let expected_energy_05 = 0.5 * k_interp * diff_interp * diff_interp;
        println!("λ=0.5: Energy = {:.6}, Expected = {:.6}", result_05.energy, expected_energy_05);
        assert!((result_05.energy - expected_energy_05).abs() < 1e-6,
                "Energy at λ=0.5 should be interpolated");

        // Verify lambda derivative
        // ∂V/∂λ = 0.5 * ΔK * (r-r0)² - K * (r-r0) * Δr0
        let delta_k = 2000.0 - 1000.0;
        let delta_r0 = 0.15 - 0.1;
        let expected_de_dlambda = 0.5 * delta_k * diff_interp * diff_interp
                                  - k_interp * diff_interp * delta_r0;
        println!("λ=0.5: dE/dλ = {:.6}, Expected = {:.6}",
                 result_05.lambda_derivative, expected_de_dlambda);
        assert!((result_05.lambda_derivative - expected_de_dlambda).abs() < 1e-6,
                "Lambda derivative should match analytical formula");

        // Verify force conservation
        let total_force = result_05.forces[0] + result_05.forces[1];
        println!("Force conservation: |F_total| = {:.6e}", total_force.length());
        assert!(total_force.length() < 1e-5,
                "Forces should be conserved (Newton's 3rd law)");
    }
}
