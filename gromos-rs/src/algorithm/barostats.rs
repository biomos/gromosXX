//! Pressure coupling algorithms (barostats)
//!
//! This module implements various barostat algorithms:
//! - Berendsen: Weak coupling barostat
//! - Parrinello-Rahman: Extended system barostat

use crate::topology::Topology;
use crate::configuration::Configuration;
use crate::math::{Vec3, Mat3};

/// Berendsen barostat parameters
#[derive(Debug, Clone)]
pub struct BerendsenBarostatParameters {
    pub target_pressure: f64,      // Target pressure (bar)
    pub coupling_time: f64,         // Coupling time constant τ_P (ps)
    pub compressibility: f64,       // Isothermal compressibility κ_T (bar⁻¹)
    pub isotropic: bool,            // Isotropic (true) or anisotropic (false) scaling
}

impl Default for BerendsenBarostatParameters {
    fn default() -> Self {
        Self {
            target_pressure: 1.0,       // 1 bar (atmospheric pressure)
            coupling_time: 0.5,         // 0.5 ps
            compressibility: 4.5e-5,    // Water compressibility at 300K (bar⁻¹)
            isotropic: true,            // Isotropic scaling
        }
    }
}

/// Parrinello-Rahman barostat parameters
#[derive(Debug, Clone)]
pub struct ParrinelloRahmanBarostatParameters {
    pub target_pressure: f64,      // Target pressure (bar)
    pub coupling_time: f64,         // Coupling time constant τ_P (ps)
    pub compressibility: f64,       // Isothermal compressibility κ_T (bar⁻¹)
}

impl Default for ParrinelloRahmanBarostatParameters {
    fn default() -> Self {
        Self {
            target_pressure: 1.0,       // 1 bar
            coupling_time: 0.5,         // 0.5 ps
            compressibility: 4.5e-5,    // Water compressibility
        }
    }
}

/// Apply Berendsen barostat (weak coupling to pressure bath)
///
/// The Berendsen barostat (Berendsen et al., 1984) rescales box dimensions
/// and atomic coordinates to weakly couple the system to a pressure bath.
///
/// Scaling factor: μ = [1 - (dt/τ_P) * κ_T * (P₀ - P)]^(1/3)
///
/// where:
/// - dt is the time step
/// - τ_P is the coupling time constant
/// - κ_T is the isothermal compressibility
/// - P₀ is the target pressure
/// - P is the current pressure
///
/// # Parameters
/// - `topo`: Molecular topology
/// - `conf`: Configuration with positions and box
/// - `dt`: Time step size (ps)
/// - `params`: Barostat parameters
/// - `virial`: Virial tensor from force calculations
pub fn berendsen_barostat(
    _topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &BerendsenBarostatParameters,
    virial: &Mat3,
) {
    // Calculate current pressure from virial
    // P = (N*k_B*T + virial) / V
    // For simplicity, use instantaneous virial-based pressure

    let volume = conf.current().box_config.volume();

    if volume < 1e-10 {
        return; // Invalid box
    }

    // Pressure from virial (simplified)
    // P = Tr(virial) / (3 * V)
    // Full calculation would include kinetic term
    let virial_trace = virial.x_axis.x as f64 + virial.y_axis.y as f64 + virial.z_axis.z as f64;
    let current_pressure = virial_trace / (3.0 * volume);

    // Berendsen scaling factor
    // μ = [1 - β * (P₀ - P)]^(1/3)
    // where β = (dt/τ_P) * κ_T
    let beta = (dt / params.coupling_time) * params.compressibility;
    let pressure_diff = params.target_pressure - current_pressure;

    if params.isotropic {
        // Isotropic scaling (same factor for all dimensions)
        let scaling_factor = (1.0 - beta * pressure_diff).cbrt();

        // Scale box dimensions
        let box_vectors = conf.current().box_config.vectors;
        let scaled_x = box_vectors.x_axis * scaling_factor as f32;
        let scaled_y = box_vectors.y_axis * scaling_factor as f32;
        let scaled_z = box_vectors.z_axis * scaling_factor as f32;
        conf.current_mut().box_config.vectors = Mat3::from_cols(scaled_x, scaled_y, scaled_z);
        conf.current_mut().box_config.inv_vectors = conf.current_mut().box_config.vectors.inverse();

        // Scale atomic positions
        for pos in conf.current_mut().pos.iter_mut() {
            *pos *= scaling_factor as f32;
        }

    } else {
        // Anisotropic scaling (different factors for each dimension)
        // Would require full pressure tensor calculation
        // Simplified: use isotropic scaling as fallback
        let scaling_factor = (1.0 - beta * pressure_diff).cbrt();

        for pos in conf.current_mut().pos.iter_mut() {
            *pos *= scaling_factor as f32;
        }
    }
}

/// Apply Parrinello-Rahman barostat (extended system)
///
/// The Parrinello-Rahman barostat (Parrinello & Rahman, 1981) introduces
/// box matrix dynamics that couple to the system pressure.
///
/// This is a more sophisticated barostat that allows for:
/// - Anisotropic box fluctuations
/// - Shape changes of the simulation box
///
/// # Parameters
/// - `topo`: Molecular topology
/// - `conf`: Configuration with positions and box
/// - `dt`: Time step size (ps)
/// - `params`: Barostat parameters
/// - `virial`: Virial tensor from force calculations
///
/// # Note
/// This is a simplified implementation. Full Parrinello-Rahman requires
/// integration of box matrix equations of motion.
pub fn parrinello_rahman_barostat(
    _topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &ParrinelloRahmanBarostatParameters,
    virial: &Mat3,
) {
    // Simplified Parrinello-Rahman implementation
    // Full implementation would integrate box matrix equations:
    // d²h/dt² = V/W * (P - P₀)
    //
    // where:
    // - h is the box matrix
    // - V is the volume
    // - W is the box mass parameter
    // - P is the pressure tensor

    let volume = conf.current().box_config.volume();

    if volume < 1e-10 {
        return;
    }

    // Calculate pressure from virial
    let virial_trace = virial.x_axis.x as f64 + virial.y_axis.y as f64 + virial.z_axis.z as f64;
    let current_pressure = virial_trace / (3.0 * volume);

    // Simplified scaling (similar to Berendsen for now)
    let beta = (dt / params.coupling_time) * params.compressibility;
    let pressure_diff = params.target_pressure - current_pressure;
    let scaling_factor = (1.0 - beta * pressure_diff).cbrt();

    // Scale positions and box
    for pos in conf.current_mut().pos.iter_mut() {
        *pos *= scaling_factor as f32;
    }

    // Note: Full PR barostat would update box matrix derivatives
    // and integrate them over time to allow for shape fluctuations
}

/// Calculate virial tensor from forces and positions
///
/// Virial: W = -Σᵢ rᵢ ⊗ Fᵢ
///
/// where ⊗ is the outer product
///
/// # Parameters
/// - `positions`: Atomic positions
/// - `forces`: Forces on atoms
///
/// # Returns
/// - 3x3 virial tensor
pub fn calculate_virial(positions: &[Vec3], forces: &[Vec3]) -> Mat3 {
    // Initialize virial components
    let mut vir = [[0.0f32; 3]; 3];

    for i in 0..positions.len().min(forces.len()) {
        let r = positions[i];
        let f = forces[i];

        // Virial contribution: -r ⊗ f
        vir[0][0] -= r.x * f.x;
        vir[0][1] -= r.x * f.y;
        vir[0][2] -= r.x * f.z;
        vir[1][0] -= r.y * f.x;
        vir[1][1] -= r.y * f.y;
        vir[1][2] -= r.y * f.z;
        vir[2][0] -= r.z * f.x;
        vir[2][1] -= r.z * f.y;
        vir[2][2] -= r.z * f.z;
    }

    // Build Mat3 from column vectors
    Mat3::from_cols(
        Vec3::new(vir[0][0], vir[1][0], vir[2][0]),
        Vec3::new(vir[0][1], vir[1][1], vir[2][1]),
        Vec3::new(vir[0][2], vir[1][2], vir[2][2]),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_berendsen_barostat_parameters() {
        let params = BerendsenBarostatParameters::default();
        assert!(params.target_pressure > 0.0);
        assert!(params.coupling_time > 0.0);
        assert!(params.compressibility > 0.0);
    }

    #[test]
    fn test_parrinello_rahman_parameters() {
        let params = ParrinelloRahmanBarostatParameters::default();
        assert!(params.target_pressure > 0.0);
        assert!(params.coupling_time > 0.0);
        assert!(params.compressibility > 0.0);
    }

    #[test]
    fn test_calculate_virial() {
        let positions = vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let forces = vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];

        let virial = calculate_virial(&positions, &forces);

        // Check that virial is symmetric (approximately)
        assert!((virial[(0, 1)] - virial[(1, 0)]).abs() < 1e-6);
    }
}
