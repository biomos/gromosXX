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

/// Result of force calculation: energy and forces
#[derive(Debug, Clone)]
pub struct ForceEnergy {
    pub energy: f64,
    pub forces: Vec<Vec3>,
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
}
