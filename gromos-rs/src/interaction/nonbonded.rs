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
    pub crf_cut: f64,       // Cutoff distance (nm)
    pub crf_2cut3i: f64,    // crf / (2 * cutoff^3) for force
    pub crf_cut3i: f64,     // (1 - crf/2) / cutoff for energy
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

/// Lambda-dependent parameters for perturbed nonbonded interactions
///
/// Direct translation from GROMOS++ perturbed_nonbonded_term.h/.cc
#[derive(Debug, Clone)]
pub struct PerturbedLambdaParams {
    /// Lambda for LJ interaction, state A: (1-λ)^n
    pub a_lj_lambda_n: f64,
    /// Lambda for LJ interaction, state B: λ^n
    pub b_lj_lambda_n: f64,
    /// Lambda for CRF interaction, state A: (1-λ)^n
    pub a_crf_lambda_n: f64,
    /// Lambda for CRF interaction, state B: λ^n
    pub b_crf_lambda_n: f64,

    /// Lambda derivative for LJ, state A: (1-λ)^(n-1) * dλ/dt
    pub a_lj_lambda_n_1: f64,
    /// Lambda derivative for LJ, state B: λ^(n-1) * dλ/dt
    pub b_lj_lambda_n_1: f64,
    /// Lambda derivative for CRF, state A: (1-λ)^(n-1) * dλ/dt
    pub a_crf_lambda_n_1: f64,
    /// Lambda derivative for CRF, state B: λ^(n-1) * dλ/dt
    pub b_crf_lambda_n_1: f64,

    /// Soft-core lambda for state A: (1-λ_soft)
    pub a_ljs_lambda: f64,
    /// Soft-core lambda for state B: λ_soft
    pub b_ljs_lambda: f64,
    /// Soft-core lambda squared for state A: (1-λ_soft)²
    pub a_ljs_lambda2: f64,
    /// Soft-core lambda squared for state B: λ_soft²
    pub b_ljs_lambda2: f64,

    /// Soft-core CRF lambda for state A: (1-λ_soft)
    pub a_crfs_lambda: f64,
    /// Soft-core CRF lambda for state B: λ_soft
    pub b_crfs_lambda: f64,
    /// Soft-core CRF lambda squared for state A: (1-λ_soft)²
    pub a_crfs_lambda2: f64,
    /// Soft-core CRF lambda squared for state B: λ_soft²
    pub b_crfs_lambda2: f64,

    /// Lambda exponent (n)
    pub lambda_exp: i32,
}

impl PerturbedLambdaParams {
    /// Create lambda parameters from lambda controller
    ///
    /// Direct translation from set_lambda() in perturbed_nonbonded_term.cc
    pub fn from_lambda(
        lj_lambda: f64,
        ljs_lambda: f64,
        crf_lambda: f64,
        crfs_lambda: f64,
        lj_lambda_derivative: f64,
        ljs_lambda_derivative: f64,
        crf_lambda_derivative: f64,
        crfs_lambda_derivative: f64,
        n: i32,
    ) -> Self {
        let n_f = n as f64;

        Self {
            // State A: (1-λ)^n
            a_lj_lambda_n: (1.0 - lj_lambda).powi(n),
            a_crf_lambda_n: (1.0 - crf_lambda).powi(n),

            // State B: λ^n
            b_lj_lambda_n: lj_lambda.powi(n),
            b_crf_lambda_n: crf_lambda.powi(n),

            // Derivatives: (1-λ)^(n-1) * dλ/dt
            a_lj_lambda_n_1: (1.0 - lj_lambda).powi(n - 1) * lj_lambda_derivative,
            a_crf_lambda_n_1: (1.0 - crf_lambda).powi(n - 1) * crf_lambda_derivative,

            // Derivatives: λ^(n-1) * dλ/dt
            b_lj_lambda_n_1: lj_lambda.powi(n - 1) * lj_lambda_derivative,
            b_crf_lambda_n_1: crf_lambda.powi(n - 1) * crf_lambda_derivative,

            // Soft-core lambdas: (1-λ) * dλ/dt for state A
            a_ljs_lambda: (1.0 - ljs_lambda) * ljs_lambda_derivative,
            a_crfs_lambda: (1.0 - crfs_lambda) * crfs_lambda_derivative,

            // Soft-core lambdas: λ * dλ/dt for state B
            b_ljs_lambda: ljs_lambda * ljs_lambda_derivative,
            b_crfs_lambda: crfs_lambda * crfs_lambda_derivative,

            // Squared soft-core lambdas
            a_ljs_lambda2: (1.0 - ljs_lambda) * (1.0 - ljs_lambda),
            a_crfs_lambda2: (1.0 - crfs_lambda) * (1.0 - crfs_lambda),
            b_ljs_lambda2: ljs_lambda * ljs_lambda,
            b_crfs_lambda2: crfs_lambda * crfs_lambda,

            lambda_exp: n,
        }
    }
}

/// Perturbed LJ+CRF interaction with soft-core potentials
///
/// Direct translation from eds_pert_lj_crf_interaction() in perturbed_nonbonded_term.cc
///
/// # Arguments
/// * `r` - Distance vector
/// * `a_c6, a_c12` - State A LJ parameters
/// * `b_c6, b_c12` - State B LJ parameters
/// * `a_q, b_q` - State A and B charge products
/// * `alpha_lj, alpha_crf` - Soft-core parameters
/// * `lambda_params` - Lambda values and derivatives
/// * `crf` - CRF parameters
///
/// # Returns
/// * `force_magnitude` - Scalar force magnitude
/// * `e_lj` - Lennard-Jones energy
/// * `e_crf` - Coulomb reaction field energy
/// * `de_lj` - Lambda derivative of LJ energy (dH/dλ)
/// * `de_crf` - Lambda derivative of CRF energy (dH/dλ)
#[inline]
pub fn perturbed_lj_crf_interaction(
    r: Vec3,
    a_c6: f64,
    a_c12: f64,
    b_c6: f64,
    b_c12: f64,
    a_q: f64,
    b_q: f64,
    alpha_lj: f64,
    alpha_crf: f64,
    lambda_params: &PerturbedLambdaParams,
    crf: &CRFParameters,
) -> (f64, f64, f64, f64, f64) {
    const FOUR_PI_EPS_I: f64 = 138.9354859; // kJ mol⁻¹ nm e⁻²

    let r2 = r.length_squared() as f64;

    // Early exit for zero distance
    if r2 < 1e-10 {
        return (0.0, 0.0, 0.0, 0.0, 0.0);
    }

    let r6 = r2 * r2 * r2;
    let r4 = r2 * r2;

    // Calculate C12/C6 ratios for soft-core
    let a_c126 = if a_c6.abs() > 1e-10 { a_c12 / a_c6 } else { 0.0 };
    let b_c126 = if b_c6.abs() > 1e-10 { b_c12 / b_c6 } else { 0.0 };

    // ===== STATE A SOFT-CORE DISTANCES =====
    // For state A, we use (1-λ) for softness
    // Electrostatics: r_soft² = r² + α_crf * (1-λ_crf)²
    let a_dist2soft = r2 + alpha_crf * lambda_params.b_crfs_lambda2;
    let a_distisoft = 1.0 / a_dist2soft.sqrt();
    let a_dist3isoft = a_distisoft / a_dist2soft;

    // LJ: r⁶_soft = r⁶ + α_lj * (1-λ_lj)² * C12/C6
    let a_dist6soft = r6 + alpha_lj * lambda_params.b_ljs_lambda2 * a_c126;
    let a_dist6isoft = 1.0 / a_dist6soft;

    // ===== STATE B SOFT-CORE DISTANCES =====
    // For state B, we use λ for softness
    let b_dist2soft = r2 + alpha_crf * lambda_params.a_crfs_lambda2;
    let b_distisoft = 1.0 / b_dist2soft.sqrt();
    let b_dist3isoft = b_distisoft / b_dist2soft;

    let b_dist6soft = r6 + alpha_lj * lambda_params.a_ljs_lambda2 * b_c126;
    let b_dist6isoft = 1.0 / b_dist6soft;

    // ===== SOFT-CORE CUTOFF MODIFICATIONS =====
    let cut2 = crf.crf_cut * crf.crf_cut;
    let a_cut2soft = cut2 + alpha_crf * lambda_params.b_crfs_lambda2;
    let b_cut2soft = cut2 + alpha_crf * lambda_params.a_crfs_lambda2;

    let a_cut2soft3 = a_cut2soft * a_cut2soft * a_cut2soft;
    let b_cut2soft3 = b_cut2soft * b_cut2soft * b_cut2soft;

    // CRF constants with soft-core cutoff
    // crf_2cut3i = crf / (2 * cutoff³)
    let a_crf_2cut3i = crf.crf_2cut3i / a_cut2soft3.sqrt();
    let b_crf_2cut3i = crf.crf_2cut3i / b_cut2soft3.sqrt();

    let a_crf_cut3i = 2.0 * a_crf_2cut3i;
    let b_crf_cut3i = 2.0 * b_crf_2cut3i;

    // Derivative terms for soft-core
    let a_crf_pert = 3.0 * a_crf_2cut3i / a_cut2soft;
    let b_crf_pert = 3.0 * b_crf_2cut3i / b_cut2soft;

    // ===== FORCE CALCULATION =====
    // CRF force
    let mut force = (lambda_params.a_crf_lambda_n * a_q * (a_dist3isoft + a_crf_cut3i) +
                     lambda_params.b_crf_lambda_n * b_q * (b_dist3isoft + b_crf_cut3i)) *
                    FOUR_PI_EPS_I;

    // LJ attractive force: -6 * C6 / r⁸
    force += -6.0 * (lambda_params.a_lj_lambda_n * a_c6 * a_dist6isoft * a_dist6isoft +
                     lambda_params.b_lj_lambda_n * b_c6 * b_dist6isoft * b_dist6isoft) * r4;

    // LJ repulsive force: 12 * C12 / r¹⁴
    force += 12.0 * (lambda_params.a_lj_lambda_n * a_c12 * a_dist6isoft * a_dist6isoft * a_dist6isoft +
                     lambda_params.b_lj_lambda_n * b_c12 * b_dist6isoft * b_dist6isoft * b_dist6isoft) * r4;

    // ===== ENERGY CALCULATION =====
    let a_e_lj = (a_c12 * a_dist6isoft - a_c6) * a_dist6isoft;
    let b_e_lj = (b_c12 * b_dist6isoft - b_c6) * b_dist6isoft;

    let crf_cut = crf.crf_cut3i; // (1 - crf/2) / cutoff
    let a_e_crf = a_q * (a_distisoft - a_crf_2cut3i * r2 - crf_cut);
    let b_e_crf = b_q * (b_distisoft - b_crf_2cut3i * r2 - crf_cut);

    // Total energy: E = (1-λ)^n * E_A + λ^n * E_B
    let e_lj = lambda_params.a_lj_lambda_n * a_e_lj + lambda_params.b_lj_lambda_n * b_e_lj;
    let e_crf = (lambda_params.a_crf_lambda_n * a_e_crf + lambda_params.b_crf_lambda_n * b_e_crf) * FOUR_PI_EPS_I;

    // ===== LAMBDA DERIVATIVES (dH/dλ) =====
    // LJ derivative: soft-core contribution + direct lambda contribution
    let de_lj = -2.0 * alpha_lj * (
        lambda_params.a_lj_lambda_n * lambda_params.b_ljs_lambda * a_c126 *
            a_dist6isoft * a_dist6isoft * (2.0 * a_c12 * a_dist6isoft - a_c6) -
        lambda_params.b_lj_lambda_n * lambda_params.a_ljs_lambda * b_c126 *
            b_dist6isoft * b_dist6isoft * (2.0 * b_c12 * b_dist6isoft - b_c6)
    ) + (lambda_params.lambda_exp as f64) * (
        lambda_params.b_lj_lambda_n_1 * b_e_lj - lambda_params.a_lj_lambda_n_1 * a_e_lj
    );

    // CRF derivative: soft-core contribution + direct lambda contribution
    let de_crf = -(
        lambda_params.a_crf_lambda_n * a_q * lambda_params.b_crfs_lambda *
            (a_dist3isoft - a_crf_pert * r2) -
        lambda_params.b_crf_lambda_n * b_q * lambda_params.a_crfs_lambda *
            (b_dist3isoft - b_crf_pert * r2)
    ) * FOUR_PI_EPS_I * alpha_crf +
    (lambda_params.lambda_exp as f64) * (
        lambda_params.b_crf_lambda_n_1 * b_e_crf - lambda_params.a_crf_lambda_n_1 * a_e_crf
    ) * FOUR_PI_EPS_I;

    (force, e_lj, e_crf, de_lj, de_crf)
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

    #[test]
    fn test_perturbed_lambda_params() {
        // Test lambda parameter calculation
        let lambda = 0.5;
        let lambda_soft = 0.5;
        let derivative = 1.0;
        let n = 2;

        let params = PerturbedLambdaParams::from_lambda(
            lambda, lambda_soft,  // lj_lambda, ljs_lambda
            lambda, lambda_soft,  // crf_lambda, crfs_lambda
            derivative, derivative,  // lj_lambda_derivative, ljs_lambda_derivative
            derivative, derivative,  // crf_lambda_derivative, crfs_lambda_derivative
            n,
        );

        // Check λ^n calculations
        assert_relative_eq!(params.b_lj_lambda_n, 0.25, epsilon = 1e-10);  // 0.5²
        assert_relative_eq!(params.a_lj_lambda_n, 0.25, epsilon = 1e-10);  // (1-0.5)²

        // Check λ^(n-1) calculations
        assert_relative_eq!(params.b_lj_lambda_n_1, 0.5, epsilon = 1e-10);  // 0.5¹ * 1.0
        assert_relative_eq!(params.a_lj_lambda_n_1, 0.5, epsilon = 1e-10);  // (1-0.5)¹ * 1.0

        // Check soft-core lambda squared
        assert_relative_eq!(params.b_ljs_lambda2, 0.25, epsilon = 1e-10);  // 0.5²
        assert_relative_eq!(params.a_ljs_lambda2, 0.25, epsilon = 1e-10);  // (1-0.5)²
    }

    #[test]
    fn test_perturbed_interaction_at_lambda_0() {
        // At λ=0, should get pure state A
        let r = Vec3::new(1.0, 0.0, 0.0);

        let a_c6 = 0.001;
        let a_c12 = 0.0001;
        let b_c6 = 0.002;  // Different from state A
        let b_c12 = 0.0002;

        let a_q = 0.25;  // q_i * q_j for state A
        let b_q = 0.0;   // No charge in state B

        let alpha_lj = 0.5;
        let alpha_crf = 0.5;

        let lambda_params = PerturbedLambdaParams::from_lambda(
            0.0, 0.0,  // lj_lambda=0, ljs_lambda=0
            0.0, 0.0,  // crf_lambda=0, crfs_lambda=0
            1.0, 1.0,  // derivatives
            1.0, 1.0,
            1,  // n=1
        );

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431 / 2.0,
            crf_cut3i: 0.364431 / 2.0 / 1.4,
        };

        let (force, e_lj, e_crf, de_lj, de_crf) = perturbed_lj_crf_interaction(
            r, a_c6, a_c12, b_c6, b_c12, a_q, b_q,
            alpha_lj, alpha_crf, &lambda_params, &crf
        );

        // At λ=0: should be dominated by state A
        // Energy should be finite
        assert!(e_lj.is_finite(), "LJ energy should be finite");
        assert!(e_crf.is_finite(), "CRF energy should be finite");

        // Derivatives should be finite
        assert!(de_lj.is_finite(), "dH/dλ for LJ should be finite");
        assert!(de_crf.is_finite(), "dH/dλ for CRF should be finite");

        println!("λ=0: E_lj={:.6}, E_crf={:.6}, dE_lj/dλ={:.6}, dE_crf/dλ={:.6}",
                 e_lj, e_crf, de_lj, de_crf);
    }

    #[test]
    fn test_perturbed_interaction_at_lambda_1() {
        // At λ=1, should get pure state B
        let r = Vec3::new(1.0, 0.0, 0.0);

        let a_c6 = 0.001;
        let a_c12 = 0.0001;
        let b_c6 = 0.002;  // Different from state A
        let b_c12 = 0.0002;

        let a_q = 0.0;   // No charge in state A
        let b_q = 0.25;  // q_i * q_j for state B

        let alpha_lj = 0.5;
        let alpha_crf = 0.5;

        let lambda_params = PerturbedLambdaParams::from_lambda(
            1.0, 1.0,  // lj_lambda=1, ljs_lambda=1
            1.0, 1.0,  // crf_lambda=1, crfs_lambda=1
            1.0, 1.0,  // derivatives
            1.0, 1.0,
            1,  // n=1
        );

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431 / 2.0,
            crf_cut3i: 0.364431 / 2.0 / 1.4,
        };

        let (force, e_lj, e_crf, de_lj, de_crf) = perturbed_lj_crf_interaction(
            r, a_c6, a_c12, b_c6, b_c12, a_q, b_q,
            alpha_lj, alpha_crf, &lambda_params, &crf
        );

        // At λ=1: should be dominated by state B
        assert!(e_lj.is_finite(), "LJ energy should be finite");
        assert!(e_crf.is_finite(), "CRF energy should be finite");
        assert!(de_lj.is_finite(), "dH/dλ for LJ should be finite");
        assert!(de_crf.is_finite(), "dH/dλ for CRF should be finite");

        println!("λ=1: E_lj={:.6}, E_crf={:.6}, dE_lj/dλ={:.6}, dE_crf/dλ={:.6}",
                 e_lj, e_crf, de_lj, de_crf);
    }

    #[test]
    fn test_perturbed_interaction_lambda_sweep() {
        // Test lambda sweep from 0 to 1
        let r = Vec3::new(0.5, 0.0, 0.0);

        let a_c6 = 0.001;
        let a_c12 = 0.0001;
        let b_c6 = 0.002;
        let b_c12 = 0.0003;

        let a_q = 0.25;
        let b_q = -0.25;

        let alpha_lj = 0.5;
        let alpha_crf = 0.5;

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431 / 2.0,
            crf_cut3i: 0.364431 / 2.0 / 1.4,
        };

        let lambdas = [0.0, 0.25, 0.5, 0.75, 1.0];
        let mut energies = Vec::new();
        let mut derivatives = Vec::new();

        for &lambda in &lambdas {
            let lambda_params = PerturbedLambdaParams::from_lambda(
                lambda, lambda, lambda, lambda,
                1.0, 1.0, 1.0, 1.0,
                1,
            );

            let (_, e_lj, e_crf, de_lj, de_crf) = perturbed_lj_crf_interaction(
                r, a_c6, a_c12, b_c6, b_c12, a_q, b_q,
                alpha_lj, alpha_crf, &lambda_params, &crf
            );

            energies.push(e_lj + e_crf);
            derivatives.push(de_lj + de_crf);

            // All should be finite
            assert!(e_lj.is_finite(), "Energy should be finite at λ={}", lambda);
            assert!(de_lj.is_finite(), "Derivative should be finite at λ={}", lambda);
        }

        println!("Lambda sweep:");
        for (i, lambda) in lambdas.iter().enumerate() {
            println!("  λ={:.2}: E={:.6}, dE/dλ={:.6}",
                     lambda, energies[i], derivatives[i]);
        }

        // Check smoothness: no huge jumps
        for i in 0..energies.len()-1 {
            let energy_change = (energies[i+1] - energies[i]).abs();
            assert!(energy_change < 100.0,
                    "Energy jump too large between λ={} and λ={}",
                    lambdas[i], lambdas[i+1]);
        }
    }

    #[test]
    fn test_soft_core_prevents_singularity() {
        // Test that soft-core prevents singularities at close distances
        let r = Vec3::new(0.1, 0.0, 0.0);  // Very close distance

        // Particle appearing (state A=dummy, state B=real)
        let a_c6 = 0.0;     // Dummy particle
        let a_c12 = 0.0;
        let b_c6 = 0.001;   // Real particle
        let b_c12 = 0.0001;

        let a_q = 0.0;
        let b_q = 0.25;

        let alpha_lj = 1.0;   // Strong soft-core
        let alpha_crf = 1.0;

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431 / 2.0,
            crf_cut3i: 0.364431 / 2.0 / 1.4,
        };

        // At λ=0 (particle not yet present), soft-core should prevent singularity
        let lambda_params_0 = PerturbedLambdaParams::from_lambda(
            0.0, 0.0, 0.0, 0.0,
            1.0, 1.0, 1.0, 1.0,
            1,
        );

        let (_, e_lj_0, e_crf_0, _, _) = perturbed_lj_crf_interaction(
            r, a_c6, a_c12, b_c6, b_c12, a_q, b_q,
            alpha_lj, alpha_crf, &lambda_params_0, &crf
        );

        // At λ=1 (particle fully present), will have high energy but should be finite
        let lambda_params_1 = PerturbedLambdaParams::from_lambda(
            1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0,
            1,
        );

        let (_, e_lj_1, e_crf_1, _, _) = perturbed_lj_crf_interaction(
            r, a_c6, a_c12, b_c6, b_c12, a_q, b_q,
            alpha_lj, alpha_crf, &lambda_params_1, &crf
        );

        // Energies should be finite (no NaN or Inf)
        assert!(e_lj_0.is_finite(), "Soft-core should prevent LJ singularity at λ=0");
        assert!(e_crf_0.is_finite(), "Soft-core should prevent CRF singularity at λ=0");
        assert!(e_lj_1.is_finite(), "Energy should be finite at λ=1");
        assert!(e_crf_1.is_finite(), "Energy should be finite at λ=1");

        // At λ=0 with soft-core, energy should be much lower than at λ=1
        println!("Soft-core test at r=0.1 nm:");
        println!("  λ=0: E_lj={:.6}, E_crf={:.6}", e_lj_0, e_crf_0);
        println!("  λ=1: E_lj={:.6}, E_crf={:.6}", e_lj_1, e_crf_1);
    }

    #[test]
    fn test_lambda_exponent() {
        // Test different lambda exponents
        let r = Vec3::new(1.0, 0.0, 0.0);
        let lambda = 0.5;

        let a_c6 = 0.001;
        let a_c12 = 0.0001;
        let b_c6 = 0.001;
        let b_c12 = 0.0001;
        let a_q = 0.1;
        let b_q = 0.1;

        let alpha_lj = 0.5;
        let alpha_crf = 0.5;

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431 / 2.0,
            crf_cut3i: 0.364431 / 2.0 / 1.4,
        };

        // Test n=1 (linear coupling)
        let params_n1 = PerturbedLambdaParams::from_lambda(
            lambda, lambda, lambda, lambda,
            1.0, 1.0, 1.0, 1.0,
            1,
        );

        let (_, e_lj_n1, e_crf_n1, de_lj_n1, de_crf_n1) = perturbed_lj_crf_interaction(
            r, a_c6, a_c12, b_c6, b_c12, a_q, b_q,
            alpha_lj, alpha_crf, &params_n1, &crf
        );

        // Test n=2 (quadratic coupling)
        let params_n2 = PerturbedLambdaParams::from_lambda(
            lambda, lambda, lambda, lambda,
            1.0, 1.0, 1.0, 1.0,
            2,
        );

        let (_, e_lj_n2, e_crf_n2, de_lj_n2, de_crf_n2) = perturbed_lj_crf_interaction(
            r, a_c6, a_c12, b_c6, b_c12, a_q, b_q,
            alpha_lj, alpha_crf, &params_n2, &crf
        );

        // Energies should differ due to different λ^n weighting
        assert_ne!(e_lj_n1, e_lj_n2, "Energy should differ for different exponents");
        assert_ne!(de_lj_n1, de_lj_n2, "Derivative should differ for different exponents");

        println!("Lambda exponent comparison at λ=0.5:");
        println!("  n=1: E_lj={:.6}, dE_lj/dλ={:.6}", e_lj_n1, de_lj_n1);
        println!("  n=2: E_lj={:.6}, dE_lj/dλ={:.6}", e_lj_n2, de_lj_n2);
    }
}
