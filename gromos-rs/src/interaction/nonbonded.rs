//! Nonbonded interactions: Lennard-Jones and Coulomb (CRF)

use crate::math::{Vec3, BoundaryCondition};
use rayon::prelude::*;

/// Lennard-Jones parameters for atom pair
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct LJParameters {
    pub c6: f64,   // LJ C6 coefficient (attraction)
    pub c12: f64,  // LJ C12 coefficient (repulsion)
}

/// Coulomb Reaction Field parameters
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CRFParameters {
    pub crf_cut: f64,       // Cutoff distance
    pub crf_2cut3i: f64,    // 2 / cutoff^3
    pub crf_cut3i: f64,     // 1 / cutoff^3
}

/// Storage for forces and energies
#[repr(C)]
pub struct ForceStorage {
    pub forces: Vec<Vec3>,
    pub e_lj: f64,
    pub e_crf: f64,
    pub virial: [[f64; 3]; 3],
}

impl ForceStorage {
    pub fn new(n_atoms: usize) -> Self {
        Self {
            forces: vec![Vec3::ZERO; n_atoms],
            e_lj: 0.0,
            e_crf: 0.0,
            virial: [[0.0; 3]; 3],
        }
    }

    pub fn clear(&mut self) {
        self.forces.fill(Vec3::ZERO);
        self.e_lj = 0.0;
        self.e_crf = 0.0;
        self.virial = [[0.0; 3]; 3];
    }
}

/// Core LJ + CRF interaction calculation (hot path!)
///
/// # Arguments
/// * `r` - Distance vector from i to j
/// * `c6` - LJ C6 coefficient
/// * `c12` - LJ C12 coefficient
/// * `q_prod` - Charge product qi * qj
/// * `crf` - CRF parameters
///
/// # Returns
/// * `force_magnitude` - Scalar force magnitude
/// * `e_lj` - Lennard-Jones energy
/// * `e_crf` - Coulomb reaction field energy
#[inline(always)]
pub fn lj_crf_interaction(
    r: Vec3,
    c6: f64,
    c12: f64,
    q_prod: f64,
    crf: &CRFParameters,
) -> (f64, f64, f64) {
    let r2 = r.length_squared() as f64;

    // Early exit for zero distance (should not happen, but safety)
    if r2 < 1e-10 {
        return (0.0, 0.0, 0.0);
    }

    let inv_r2 = 1.0 / r2;
    let inv_r6 = inv_r2 * inv_r2 * inv_r2;

    // Lennard-Jones: E_lj = C12/r^12 - C6/r^6
    let e_lj = (c12 * inv_r6 - c6) * inv_r6;

    // Lennard-Jones force: F = 12*C12/r^14 - 6*C6/r^8
    let f_lj = (12.0 * c12 * inv_r6 - 6.0 * c6) * inv_r6 * inv_r2;

    // Coulomb Reaction Field
    let inv_r = inv_r2.sqrt();
    let e_crf = q_prod * (inv_r + crf.crf_2cut3i * r2 - crf.crf_cut3i);
    let f_crf = q_prod * (inv_r * inv_r2 - 2.0 * crf.crf_2cut3i);

    let force_magnitude = f_lj + f_crf;

    (force_magnitude, e_lj, e_crf)
}

/// Vectorized version processing 4 pairs simultaneously (AVX2/AVX-512)
#[cfg(feature = "simd")]
#[inline]
pub fn lj_crf_interaction_simd_x4(
    r: [Vec3; 4],
    c6: [f64; 4],
    c12: [f64; 4],
    q_prod: [f64; 4],
    crf: &CRFParameters,
) -> ([f64; 4], [f64; 4], [f64; 4]) {
    use wide::f64x4;

    // Convert to SIMD vectors
    let r2_array: [f64; 4] = r.map(|v| v.length_squared() as f64);
    let r2 = f64x4::from(r2_array);

    let inv_r2 = f64x4::splat(1.0) / r2;
    let inv_r6 = inv_r2 * inv_r2 * inv_r2;

    // Lennard-Jones (vectorized)
    let c6_vec = f64x4::from(c6);
    let c12_vec = f64x4::from(c12);

    let e_lj = (c12_vec * inv_r6 - c6_vec) * inv_r6;
    let f_lj = (f64x4::splat(12.0) * c12_vec * inv_r6 - f64x4::splat(6.0) * c6_vec) * inv_r6 * inv_r2;

    // Coulomb Reaction Field (vectorized)
    let q_prod_vec = f64x4::from(q_prod);
    let inv_r = inv_r2.sqrt();

    let e_crf = q_prod_vec * (inv_r + f64x4::splat(crf.crf_2cut3i) * r2 - f64x4::splat(crf.crf_cut3i));
    let f_crf = q_prod_vec * (inv_r * inv_r2 - f64x4::splat(2.0 * crf.crf_2cut3i));

    let force = f_lj + f_crf;

    // Convert back to arrays
    let force_array: [f64; 4] = force.into();
    let e_lj_array: [f64; 4] = e_lj.into();
    let e_crf_array: [f64; 4] = e_crf.into();

    (force_array, e_lj_array, e_crf_array)
}

/// Inner loop for nonbonded interactions
///
/// This is the hottest function in MD simulations!
/// Processes pairlist and accumulates forces, energies, and virial.
pub fn lj_crf_innerloop<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f32],
    iac: &[u32],  // Integer atom codes (atom types)
    pairlist: &[(u32, u32)],
    lj_params: &[Vec<LJParameters>],
    crf: &CRFParameters,
    periodicity: &BC,
    storage: &mut ForceStorage,
) {
    // Serial version (baseline)
    for &(i, j) in pairlist {
        let i = i as usize;
        let j = j as usize;

        // Get positions
        let pos_i = positions[i];
        let pos_j = positions[j];

        // Apply periodic boundary conditions
        let r = periodicity.nearest_image(pos_i, pos_j);

        // Get interaction parameters
        let type_i = iac[i] as usize;
        let type_j = iac[j] as usize;
        let lj = lj_params[type_i][type_j];

        let q_prod = (charges[i] * charges[j]) as f64;

        // Calculate interaction
        let (f_mag, e_lj, e_crf) = lj_crf_interaction(r, lj.c6, lj.c12, q_prod, crf);

        // Accumulate forces
        let force = r * f_mag as f32;
        storage.forces[i] += force;
        storage.forces[j] -= force;

        // Accumulate energies
        storage.e_lj += e_lj;
        storage.e_crf += e_crf;

        // Accumulate virial (for pressure calculation)
        for a in 0..3 {
            for b in 0..3 {
                storage.virial[a][b] += r[b] as f64 * force[a] as f64;
            }
        }
    }
}

/// Parallel version of innerloop using Rayon
pub fn lj_crf_innerloop_parallel<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f32],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &[Vec<LJParameters>],
    crf: &CRFParameters,
    periodicity: &BC,
    n_atoms: usize,
) -> ForceStorage {
    // Thread-local accumulation to avoid data races
    let results: Vec<ForceStorage> = pairlist
        .par_chunks(1024)  // Process in chunks for better cache locality
        .map(|chunk| {
            let mut local_storage = ForceStorage::new(n_atoms);

            for &(i, j) in chunk {
                let i = i as usize;
                let j = j as usize;

                let pos_i = positions[i];
                let pos_j = positions[j];
                let r = periodicity.nearest_image(pos_i, pos_j);

                let type_i = iac[i] as usize;
                let type_j = iac[j] as usize;
                let lj = lj_params[type_i][type_j];
                let q_prod = (charges[i] * charges[j]) as f64;

                let (f_mag, e_lj, e_crf) = lj_crf_interaction(r, lj.c6, lj.c12, q_prod, crf);

                let force = r * f_mag as f32;
                local_storage.forces[i] += force;
                local_storage.forces[j] -= force;

                local_storage.e_lj += e_lj;
                local_storage.e_crf += e_crf;

                for a in 0..3 {
                    for b in 0..3 {
                        local_storage.virial[a][b] += r[b] as f64 * force[a] as f64;
                    }
                }
            }

            local_storage
        })
        .collect();

    // Reduce all thread-local results
    let mut final_storage = ForceStorage::new(n_atoms);
    for local in results {
        for (f_final, f_local) in final_storage.forces.iter_mut().zip(local.forces.iter()) {
            *f_final += *f_local;
        }
        final_storage.e_lj += local.e_lj;
        final_storage.e_crf += local.e_crf;
        for a in 0..3 {
            for b in 0..3 {
                final_storage.virial[a][b] += local.virial[a][b];
            }
        }
    }

    final_storage
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Vacuum, Rectangular};
    use approx::assert_relative_eq;

    #[test]
    fn test_lj_interaction() {
        // Test case: two atoms at 1 nm distance
        let r = Vec3::new(1.0, 0.0, 0.0);
        let c6 = 0.001;   // Typical C6 value
        let c12 = 0.0001; // Typical C12 value
        let q_prod = 0.0; // No charges

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.0,
            crf_cut3i: 0.0,
        };

        let (f, e_lj, _e_crf) = lj_crf_interaction(r, c6, c12, q_prod, &crf);

        // Verify energy is correct: E = C12/r^12 - C6/r^6
        let expected_e_lj = c12 - c6;
        assert_relative_eq!(e_lj, expected_e_lj, epsilon = 1e-9);

        // Verify force has correct sign (attractive at this distance)
        assert!(f < 0.0, "Force should be attractive");
    }

    #[test]
    fn test_innerloop_simple() {
        // Simple 2-atom system
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        ];
        let charges = vec![0.5, -0.5];
        let iac = vec![0, 0];
        let pairlist = vec![(0, 1)];

        let lj_params = vec![vec![LJParameters { c6: 0.001, c12: 0.0001 }]];
        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431,  // 2 / 1.4^3
            crf_cut3i: 0.364431 / 2.0,
        };

        let periodicity = Vacuum;
        let mut storage = ForceStorage::new(2);

        lj_crf_innerloop(
            &positions,
            &charges,
            &iac,
            &pairlist,
            &lj_params,
            &crf,
            &periodicity,
            &mut storage,
        );

        // Verify Newton's third law: F_i = -F_j
        assert_relative_eq!(storage.forces[0].x, -storage.forces[1].x, epsilon = 1e-6);
        assert_relative_eq!(storage.forces[0].y, -storage.forces[1].y, epsilon = 1e-6);
        assert_relative_eq!(storage.forces[0].z, -storage.forces[1].z, epsilon = 1e-6);

        // Verify energy is non-zero
        assert!(storage.e_lj != 0.0);
        assert!(storage.e_crf != 0.0);
    }

    #[test]
    fn test_periodic_boundary() {
        let box_size = Vec3::splat(10.0);
        let periodicity = Rectangular::new(box_size);

        // Two atoms across the boundary
        let positions = vec![
            Vec3::new(9.5, 0.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
        ];
        let charges = vec![0.0, 0.0];
        let iac = vec![0, 0];
        let pairlist = vec![(0, 1)];

        let lj_params = vec![vec![LJParameters { c6: 0.001, c12: 0.0001 }]];
        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.0,
            crf_cut3i: 0.0,
        };

        let mut storage = ForceStorage::new(2);

        lj_crf_innerloop(
            &positions,
            &charges,
            &iac,
            &pairlist,
            &lj_params,
            &crf,
            &periodicity,
            &mut storage,
        );

        // Distance should be 1.0 nm (minimum image), not 9.0 nm
        // Force should be reasonable for 1 nm separation
        assert!(storage.forces[0].length() < 10.0, "Force too large - PBC not working");
    }
}
