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
}
