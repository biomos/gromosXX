//! Long-range electrostatics methods
//!
//! This module implements various methods for handling long-range electrostatic interactions:
//! - Reaction Field (RF): Continuum dielectric approximation
//! - Particle Mesh Ewald (PME): FFT-based periodic electrostatics
//! - P3M: Particle-Particle Particle-Mesh method

use crate::math::{Vec3, Mat3};
use crate::topology::Topology;
use crate::configuration::Configuration;

/// 1/√π constant (stable alternative to unstable FRAC_1_SQRT_PI)
const INV_SQRT_PI: f64 = 0.5641895835477562869480794515607725858;

/// Electrostatics method selection
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElectrostaticsMethod {
    /// No long-range electrostatics (cutoff only)
    Cutoff,
    /// Reaction Field with continuum dielectric
    ReactionField,
    /// Particle Mesh Ewald (FFT-based)
    PME,
    /// Particle-Particle Particle-Mesh
    P3M,
}

/// Reaction Field parameters
///
/// The Reaction Field method treats the region beyond the cutoff as a
/// continuum dielectric with permittivity epsilon_rf.
///
/// # Theory
/// The electrostatic potential includes a reaction field correction:
///
/// V(r) = q_i * q_j * [1/r + C_rf * r² - C_cut]
///
/// where:
/// - C_rf = crf / (2 * cutoff³)
/// - C_cut = (1 + crf/2) / cutoff
/// - crf depends on epsilon, epsilon_rf, and kappa (ionic strength)
///
/// # Special Cases
/// - epsilon_rf = ∞: crf = -1 (conducting boundary)
/// - epsilon_rf = 1: crf = 0 (vacuum, reduces to cutoff)
/// - epsilon_rf = epsilon: crf = 0 (homogeneous medium)
///
/// # References
/// - Tironi et al., J. Chem. Phys. 1995, 102, 5451
/// - Barker & Watts, Mol. Phys. 1973, 26, 789
#[derive(Debug, Clone)]
pub struct ReactionFieldParameters {
    /// Cutoff distance (nm)
    pub cutoff: f64,

    /// Dielectric constant of the medium (typically 1.0 for vacuum)
    pub epsilon: f64,

    /// Dielectric constant outside the cutoff sphere
    /// - Use 0.0 for infinity (conducting boundary, GROMOS default)
    /// - Use 78.4 for water
    /// - Use 1.0 for vacuum
    pub epsilon_rf: f64,

    /// Ionic strength parameter kappa (nm⁻¹)
    /// For non-ionic systems, set to 0.0
    pub kappa: f64,

    // Derived parameters (computed from above)
    /// CRF coefficient
    pub crf: f64,

    /// 1 / cutoff³
    pub cut3i: f64,

    /// crf / (2 * cutoff³)
    pub crf_2cut3i: f64,

    /// crf / cutoff³
    pub crf_cut3i: f64,

    /// (1 + crf/2) / cutoff (for energy)
    pub crf_cut: f64,
}

impl ReactionFieldParameters {
    /// Create Reaction Field parameters with proper dielectric treatment
    ///
    /// # Arguments
    /// - `cutoff`: Cutoff distance in nm
    /// - `epsilon`: Dielectric constant of the medium (typically 1.0)
    /// - `epsilon_rf`: Dielectric constant beyond cutoff (0.0 = infinity)
    /// - `kappa`: Ionic strength parameter in nm⁻¹ (0.0 for non-ionic)
    ///
    /// # Examples
    /// ```
    /// // GROMOS default: conducting boundary (epsilon_rf = infinity)
    /// let rf = ReactionFieldParameters::new(1.4, 1.0, 0.0, 0.0);
    ///
    /// // Water dielectric
    /// let rf = ReactionFieldParameters::new(1.4, 1.0, 78.4, 0.0);
    /// ```
    pub fn new(cutoff: f64, epsilon: f64, epsilon_rf: f64, kappa: f64) -> Self {
        let cut3i = 1.0 / (cutoff * cutoff * cutoff);

        // Calculate CRF coefficient based on dielectric constants
        let crf = if epsilon_rf == 0.0 {
            // Infinite dielectric (conducting boundary) - GROMOS default
            // This is the simplest and most common case
            -1.0
        } else {
            // General Reaction Field formula
            // Tironi et al., J. Chem. Phys. 1995, 102, 5451
            //
            // crf = [2*(ε - ε_rf)*(1 + κR_c) - ε_rf*(κR_c)²] /
            //       [(ε + 2ε_rf)*(1 + κR_c) + ε_rf*(κR_c)²]

            let kappa_cutoff = kappa * cutoff;
            let kappa_cutoff_sq = kappa_cutoff * kappa_cutoff;

            let numerator = 2.0 * (epsilon - epsilon_rf) * (1.0 + kappa_cutoff)
                - epsilon_rf * kappa_cutoff_sq;

            let denominator = (epsilon + 2.0 * epsilon_rf) * (1.0 + kappa_cutoff)
                + epsilon_rf * kappa_cutoff_sq;

            numerator / denominator
        };

        let crf_cut3i = crf * cut3i;
        let crf_2cut3i = crf * cut3i / 2.0;
        let crf_cut = (1.0 - crf / 2.0) / cutoff;

        Self {
            cutoff,
            epsilon,
            epsilon_rf,
            kappa,
            crf,
            cut3i,
            crf_2cut3i,
            crf_cut3i,
            crf_cut,
        }
    }

    /// Create default GROMOS RF parameters (conducting boundary)
    pub fn gromos_default(cutoff: f64) -> Self {
        Self::new(cutoff, 1.0, 0.0, 0.0)
    }

    /// Create RF parameters for aqueous solution
    pub fn aqueous(cutoff: f64) -> Self {
        Self::new(cutoff, 1.0, 78.4, 0.0)
    }
}

impl Default for ReactionFieldParameters {
    fn default() -> Self {
        // GROMOS default: 1.4 nm cutoff, conducting boundary
        Self::gromos_default(1.4)
    }
}

/// Calculate Reaction Field energy and force for atom pair
///
/// # Arguments
/// - `r_vec`: Distance vector from i to j
/// - `q_i`: Charge on atom i (in e)
/// - `q_j`: Charge on atom j (in e)
/// - `rf`: Reaction Field parameters
///
/// # Returns
/// - `(force_vec, energy)`: Force vector on atom i and electrostatic energy
///
/// # Formulas
/// Energy: E = q_i * q_j * [1/r + C_rf * r² - C_cut]
/// Force:  F = -q_i * q_j * [1/r³ - 2*C_rf] * r_vec
///
/// where C_rf = crf/(2*cutoff³) and C_cut = (1 + crf/2)/cutoff
#[inline(always)]
pub fn reaction_field_interaction(
    r_vec: Vec3,
    q_i: f64,
    q_j: f64,
    rf: &ReactionFieldParameters,
) -> (Vec3, f64) {
    let r2 = r_vec.length_squared() as f64;

    if r2 < 1e-10 {
        return (Vec3::ZERO, 0.0);
    }

    let r = r2.sqrt();
    let inv_r = 1.0 / r;
    let inv_r3 = inv_r / r2;

    let q_prod = q_i * q_j;

    // Energy: E = q_i * q_j * [1/r + crf/(2*cutoff³) * r² - (1 + crf/2)/cutoff]
    let energy = q_prod * (inv_r + rf.crf_2cut3i * r2 - rf.crf_cut);

    // Force magnitude: |F| = q_i * q_j * [1/r³ - 2*crf/(2*cutoff³)]
    let force_mag = q_prod * (inv_r3 - 2.0 * rf.crf_2cut3i);

    // Force vector: F = force_mag * r_vec (pointing from i to j)
    let force = r_vec * (force_mag as f32);

    (force, energy)
}

/// Particle Mesh Ewald (PME) parameters
///
/// PME splits the electrostatic interaction into:
/// 1. Real-space sum (short-range, Gaussian screened)
/// 2. Reciprocal-space sum (long-range, computed via FFT on a grid)
/// 3. Self-energy correction
///
/// # Theory
/// The Ewald sum decomposes 1/r into:
/// - erfc(α*r)/r  (real space, short-range)
/// - erf(α*r)/r   (reciprocal space, long-range via FFT)
///
/// where α (alpha) is the Ewald splitting parameter.
///
/// # References
/// - Darden et al., J. Chem. Phys. 1993, 98, 10089 (PME)
/// - Essmann et al., J. Chem. Phys. 1995, 103, 8577 (smooth PME)
#[derive(Debug, Clone)]
pub struct PMEParameters {
    /// Real-space cutoff (nm)
    pub cutoff: f64,

    /// Ewald splitting parameter alpha (nm⁻¹)
    /// Typically alpha ≈ 3.0 / cutoff
    pub alpha: f64,

    /// FFT grid size in x direction
    pub grid_x: usize,

    /// FFT grid size in y direction
    pub grid_y: usize,

    /// FFT grid size in z direction
    pub grid_z: usize,

    /// Interpolation order (typically 4 for cubic splines)
    pub spline_order: usize,

    /// Tolerance for energy error (typically 1e-5)
    pub tolerance: f64,
}

impl PMEParameters {
    /// Create PME parameters with automatic alpha selection
    ///
    /// # Arguments
    /// - `cutoff`: Real-space cutoff in nm
    /// - `grid_size`: Approximate number of grid points per direction
    /// - `spline_order`: Interpolation order (typically 4)
    /// - `tolerance`: Energy error tolerance (typically 1e-5)
    pub fn new(cutoff: f64, grid_size: usize, spline_order: usize, tolerance: f64) -> Self {
        // Optimal alpha for given cutoff (Essmann et al.)
        let alpha = 3.0 / cutoff;

        Self {
            cutoff,
            alpha,
            grid_x: grid_size,
            grid_y: grid_size,
            grid_z: grid_size,
            spline_order,
            tolerance,
        }
    }

    /// Create PME parameters with specified alpha
    pub fn new_with_alpha(
        cutoff: f64,
        alpha: f64,
        grid_x: usize,
        grid_y: usize,
        grid_z: usize,
        spline_order: usize,
        tolerance: f64,
    ) -> Self {
        Self {
            cutoff,
            alpha,
            grid_x,
            grid_y,
            grid_z,
            spline_order,
            tolerance,
        }
    }
}

impl Default for PMEParameters {
    fn default() -> Self {
        // Typical PME settings for ~1.0 nm cutoff
        Self::new(1.0, 64, 4, 1e-5)
    }
}

/// PME real-space interaction (Gaussian-screened Coulomb)
///
/// Computes the short-range part of PME using the complementary error function.
///
/// V(r) = q_i * q_j * erfc(α*r) / r
/// F(r) = q_i * q_j * [erfc(α*r)/r³ + 2α/√π * exp(-α²r²)/r²] * r_vec
///
/// where erfc is the complementary error function: erfc(x) = 1 - erf(x)
#[inline(always)]
pub fn pme_real_space_interaction(
    r_vec: Vec3,
    q_i: f64,
    q_j: f64,
    pme: &PMEParameters,
) -> (Vec3, f64) {
    let r2 = r_vec.length_squared() as f64;

    if r2 < 1e-10 || r2 > pme.cutoff * pme.cutoff {
        return (Vec3::ZERO, 0.0);
    }

    let r = r2.sqrt();
    let inv_r = 1.0 / r;

    let alpha_r = pme.alpha * r;
    let alpha_r_sq = alpha_r * alpha_r;

    // erfc(α*r) using approximation or libm
    let erfc_alpha_r = erfc(alpha_r);

    let q_prod = q_i * q_j;

    // Energy: E = q_i * q_j * erfc(α*r) / r
    let energy = q_prod * erfc_alpha_r * inv_r;

    // Force magnitude: |F| = q_i * q_j * [erfc(α*r)/r³ + 2α/√π * exp(-α²r²)/r²]
    let exp_term = (-alpha_r_sq).exp();
    let force_mag = q_prod * inv_r * (erfc_alpha_r / r2 +
        2.0 * pme.alpha * INV_SQRT_PI * exp_term / r);

    let force = r_vec * (force_mag as f32);

    (force, energy)
}

/// Complementary error function: erfc(x) = 1 - erf(x)
///
/// Uses a rational approximation for efficiency.
/// Accurate to ~1e-7 for all x.
#[inline]
fn erfc(x: f64) -> f64 {
    // For x < 0, use erfc(-x) = 2 - erfc(x)
    let x_abs = x.abs();

    // Rational approximation (Abramowitz & Stegun 7.1.26)
    // erfc(x) ≈ t * exp(-x²) * P(t) where t = 1/(1 + p*x)
    let p = 0.3275911;
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;

    let t = 1.0 / (1.0 + p * x_abs);
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    let poly = a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5;
    let result = poly * (-x_abs * x_abs).exp();

    if x >= 0.0 {
        result
    } else {
        2.0 - result
    }
}

//
// PME Reciprocal Space Implementation
//

/// B-spline charge assignment
///
/// Assigns particle charges to a 3D grid using B-spline interpolation.
/// The spline order determines smoothness (4 = cubic, most common).
///
/// # Algorithm (Essmann et al., 1995)
/// 1. For each particle at position r:
/// 2. Find nearest grid point
/// 3. Compute B-spline weights for nearby points
/// 4. Distribute charge according to weights
fn assign_charges_to_grid(
    positions: &[Vec3],
    charges: &[f64],
    box_vectors: &Mat3,
    grid_x: usize,
    grid_y: usize,
    grid_z: usize,
    spline_order: usize,
) -> Vec<Vec<Vec<f64>>> {
    let mut grid = vec![vec![vec![0.0; grid_z]; grid_y]; grid_x];

    // Box dimensions
    let box_x = box_vectors.x_axis.x as f64;
    let box_y = box_vectors.y_axis.y as f64;
    let box_z = box_vectors.z_axis.z as f64;

    // Grid spacing
    let h_x = box_x / grid_x as f64;
    let h_y = box_y / grid_y as f64;
    let h_z = box_z / grid_z as f64;

    for (i, pos) in positions.iter().enumerate() {
        let charge = charges[i];

        // Fractional coordinates on grid
        let u_x = (pos.x as f64 / h_x) % grid_x as f64;
        let u_y = (pos.y as f64 / h_y) % grid_y as f64;
        let u_z = (pos.z as f64 / h_z) % grid_z as f64;

        // Nearest grid point (for even spline order, use grid center)
        let k_x = if spline_order % 2 == 0 {
            u_x.floor() as i32
        } else {
            u_x.round() as i32
        };
        let k_y = if spline_order % 2 == 0 {
            u_y.floor() as i32
        } else {
            u_y.round() as i32
        };
        let k_z = if spline_order % 2 == 0 {
            u_z.floor() as i32
        } else {
            u_z.round() as i32
        };

        // Compute B-spline weights and assign charge
        let lower = -(spline_order as i32 - 1) / 2;
        let upper = spline_order as i32 / 2;

        for dx in lower..=upper {
            let wx = b_spline(spline_order, u_x - (k_x + dx) as f64);
            for dy in lower..=upper {
                let wy = b_spline(spline_order, u_y - (k_y + dy) as f64);
                for dz in lower..=upper {
                    let wz = b_spline(spline_order, u_z - (k_z + dz) as f64);

                    let grid_i = ((k_x + dx).rem_euclid(grid_x as i32)) as usize;
                    let grid_j = ((k_y + dy).rem_euclid(grid_y as i32)) as usize;
                    let grid_k = ((k_z + dz).rem_euclid(grid_z as i32)) as usize;

                    grid[grid_i][grid_j][grid_k] += charge * wx * wy * wz;
                }
            }
        }
    }

    grid
}

/// Cardinal B-spline of order n evaluated at x
///
/// B-splines provide smooth interpolation for charge assignment.
/// Order 4 (cubic) is the most commonly used in PME.
fn b_spline(order: usize, x: f64) -> f64 {
    let u = x.abs();

    match order {
        2 => {
            // Linear B-spline
            if u < 1.0 {
                1.0 - u
            } else {
                0.0
            }
        }
        3 => {
            // Quadratic B-spline
            if u < 0.5 {
                0.75 - u * u
            } else if u < 1.5 {
                let t = 1.5 - u;
                0.5 * t * t
            } else {
                0.0
            }
        }
        4 => {
            // Cubic B-spline (most common for PME)
            if u < 1.0 {
                let u2 = u * u;
                let u3 = u2 * u;
                2.0 / 3.0 - u2 + 0.5 * u3
            } else if u < 2.0 {
                let t = 2.0 - u;
                t * t * t / 6.0
            } else {
                0.0
            }
        }
        5 => {
            // Quartic B-spline
            if u < 0.5 {
                let u2 = u * u;
                115.0 / 192.0 - 5.0 / 8.0 * u2 + 1.0 / 4.0 * u2 * u2
            } else if u < 1.5 {
                let t = u - 0.5;
                let t2 = t * t;
                19.0 / 96.0 + 1.0 / 4.0 * t + 1.0 / 4.0 * t2 - 1.0 / 6.0 * t2 * t - 1.0 / 6.0 * t2 * t2
            } else if u < 2.5 {
                let t = 2.5 - u;
                t * t * t * t / 24.0
            } else {
                0.0
            }
        }
        _ => {
            eprintln!("Warning: B-spline order {} not implemented, using cubic", order);
            b_spline(4, x)
        }
    }
}

/// PME influence function (optimal Gaussian in Fourier space)
///
/// This function modulates the Fourier-transformed charge density
/// to account for the finite grid spacing and spline interpolation.
///
/// # Formula
/// G(k) = exp(-π²k²/α²) / (V * k²) * M(k)
///
/// where M(k) is the Fourier transform of the interpolation function.
fn pme_influence_function(
    k_x: f64,
    k_y: f64,
    k_z: f64,
    alpha: f64,
    volume: f64,
    spline_order: usize,
) -> f64 {
    let k2 = k_x * k_x + k_y * k_y + k_z * k_z;

    if k2 < 1e-10 {
        return 0.0; // Skip k=0 term
    }

    let pi_k_over_alpha = std::f64::consts::PI * k2.sqrt() / alpha;

    // Gaussian factor: exp(-π²k²/α²)
    let gaussian = (-pi_k_over_alpha * pi_k_over_alpha).exp();

    // Spline modulation factor (simplified for cubic splines)
    // Full implementation would use the analytical Fourier transform of B-splines
    let modulation = 1.0; // Placeholder - would need proper M(k) calculation

    gaussian / (volume * k2) * modulation
}

/// Complete PME calculation (simplified structure)
///
/// This function shows the structure of a full PME calculation.
/// For production use, integrate with an FFT library like `rustfft`.
///
/// # Steps
/// 1. Assign charges to grid (B-splines)
/// 2. FFT forward: charge grid → Fourier space
/// 3. Apply influence function and compute reciprocal energy
/// 4. FFT inverse: Fourier space → force grid
/// 5. Interpolate forces back to atoms
///
/// # Note
/// This is a **structural template**. Full FFT implementation requires
/// the `rustfft` crate or similar. Add to Cargo.toml:
/// ```toml
/// rustfft = "6.1"
/// ```
pub fn pme_reciprocal_space(
    positions: &[Vec3],
    charges: &[f64],
    box_vectors: &Mat3,
    pme: &PMEParameters,
) -> (Vec<Vec3>, f64) {
    let n_atoms = positions.len();

    // 1. Assign charges to grid
    let charge_grid = assign_charges_to_grid(
        positions,
        charges,
        box_vectors,
        pme.grid_x,
        pme.grid_y,
        pme.grid_z,
        pme.spline_order,
    );

    // 2. FFT forward transform (requires FFT library)
    // In production: use rustfft to transform charge_grid
    // let mut planner = FftPlanner::new();
    // let fft = planner.plan_fft_forward(pme.grid_x);
    // ... apply 3D FFT ...

    // 3. Apply influence function in Fourier space
    let volume = box_vectors.determinant() as f64;
    let box_x = box_vectors.x_axis.x as f64;
    let box_y = box_vectors.y_axis.y as f64;
    let box_z = box_vectors.z_axis.z as f64;

    let mut reciprocal_energy = 0.0;

    // Loop over k-vectors (simplified - would be done in FFT space)
    for i_x in 0..pme.grid_x {
        let k_x = 2.0 * std::f64::consts::PI * i_x as f64 / box_x;
        for i_y in 0..pme.grid_y {
            let k_y = 2.0 * std::f64::consts::PI * i_y as f64 / box_y;
            for i_z in 0..pme.grid_z {
                let k_z = 2.0 * std::f64::consts::PI * i_z as f64 / box_z;

                // Get Fourier coefficient (from FFT, here using grid directly as placeholder)
                let rho_k = charge_grid[i_x][i_y][i_z];

                // Influence function
                let g_k = pme_influence_function(k_x, k_y, k_z, pme.alpha, volume, pme.spline_order);

                // Reciprocal space energy contribution
                reciprocal_energy += g_k * rho_k * rho_k;
            }
        }
    }

    // 4. FFT inverse transform for forces (requires FFT library)
    // ... inverse FFT to get force grid ...

    // 5. Interpolate forces back to atoms
    let forces = vec![Vec3::ZERO; n_atoms];
    // ... interpolate from force grid to atomic forces ...

    (forces, reciprocal_energy)
}

/// PME self-energy correction
///
/// The PME sum includes spurious self-interactions that must be removed.
///
/// E_self = -α/√π * Σ_i q_i²
pub fn pme_self_energy(charges: &[f64], alpha: f64) -> f64 {
    let mut self_energy = 0.0;
    for &q in charges {
        self_energy += q * q;
    }
    -alpha * INV_SQRT_PI * self_energy
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rf_parameters_conducting_boundary() {
        let rf = ReactionFieldParameters::gromos_default(1.4);
        assert_eq!(rf.cutoff, 1.4);
        assert_eq!(rf.epsilon, 1.0);
        assert_eq!(rf.epsilon_rf, 0.0);
        assert_eq!(rf.crf, -1.0);
    }

    #[test]
    fn test_rf_parameters_aqueous() {
        let rf = ReactionFieldParameters::aqueous(1.4);
        assert_eq!(rf.cutoff, 1.4);
        assert_eq!(rf.epsilon_rf, 78.4);
        // CRF should be negative but not -1
        assert!(rf.crf < 0.0 && rf.crf > -1.0);
    }

    #[test]
    fn test_rf_interaction() {
        let rf = ReactionFieldParameters::gromos_default(1.4);
        let r_vec = Vec3::new(1.0, 0.0, 0.0);
        let q_i = 1.0;  // +1 e
        let q_j = -1.0; // -1 e

        let (force, energy) = reaction_field_interaction(r_vec, q_i, q_j, &rf);

        // Energy should be negative (attractive)
        assert!(energy < 0.0);

        // Force should point from + to - (attractive)
        assert!(force.x < 0.0);
    }

    #[test]
    fn test_erfc_function() {
        // Test known values
        assert!((erfc(0.0) - 1.0).abs() < 1e-6);
        assert!((erfc(1.0) - 0.1572992).abs() < 1e-5);
        assert!((erfc(2.0) - 0.0046777).abs() < 1e-5);

        // Test symmetry
        let x = 1.5;
        assert!((erfc(-x) - (2.0 - erfc(x))).abs() < 1e-6);
    }

    #[test]
    fn test_pme_real_space() {
        let pme = PMEParameters::default();
        let r_vec = Vec3::new(0.5, 0.0, 0.0);
        let q_i = 1.0;
        let q_j = -1.0;

        let (force, energy) = pme_real_space_interaction(r_vec, q_i, q_j, &pme);

        // Energy should be negative (attractive) but smaller than full Coulomb
        assert!(energy < 0.0);

        // Force should be attractive
        assert!(force.x < 0.0);
    }
}
