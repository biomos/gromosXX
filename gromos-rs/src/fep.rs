//! Free Energy Perturbation (FEP) module
//!
//! Direct translation of GROMOS FEP capabilities from:
//! - md++/src/algorithm/integration/slow_growth.cc
//! - md++/src/topology/perturbed_atom.h
//! - md++/src/interaction/nonbonded/interaction/perturbed_nonbonded_term.h
//!
//! # Overview
//!
//! Free energy perturbation allows calculation of free energy differences by
//! gradually transforming one molecular state (A) into another (B) using a
//! coupling parameter λ ∈ [0,1]:
//!
//! - λ = 0: Pure state A
//! - λ = 1: Pure state B
//! - 0 < λ < 1: Hybrid state
//!
//! # Methods
//!
//! - **Thermodynamic Integration (TI)**: ΔG = ∫₀¹ ⟨∂H/∂λ⟩_λ dλ
//! - **Slow Growth**: Continuous λ evolution, λ(t) = λ₀ + (dλ/dt)·t
//! - **Staged Windows**: Multiple simulations at fixed λ values
//!
//! # Features
//!
//! - Dual-state topology (A/B parameters for each perturbed atom)
//! - Lambda-dependent parameter interpolation
//! - Soft-core potentials (prevent singularities)
//! - Individual lambda values per interaction type
//! - Energy group-specific scaling

use crate::topology::Topology;
use crate::configuration::Configuration;

/// Lambda controller for free energy perturbation
///
/// Manages the coupling parameter λ and its evolution over time
#[derive(Debug, Clone)]
pub struct LambdaController {
    /// Current lambda value (0 to 1)
    pub lambda: f64,
    /// Lambda change rate (dλ/dt) in ps⁻¹
    pub dlambda_dt: f64,
    /// Lambda exponent (for λⁿ coupling)
    pub lambda_exponent: i32,
    /// Individual lambda values per interaction type
    pub interaction_lambdas: InteractionLambdas,
}

/// Individual lambda values for different interaction types
///
/// Allows independent control of different terms (bonds, LJ, electrostatics, etc.)
#[derive(Debug, Clone)]
pub struct InteractionLambdas {
    pub bond: f64,
    pub angle: f64,
    pub dihedral: f64,
    pub improper: f64,
    pub lj: f64,
    pub lj_softness: f64,
    pub crf: f64,
    pub crf_softness: f64,
    pub mass: f64,
}

impl Default for InteractionLambdas {
    fn default() -> Self {
        Self {
            bond: 0.0,
            angle: 0.0,
            dihedral: 0.0,
            improper: 0.0,
            lj: 0.0,
            lj_softness: 0.0,
            crf: 0.0,
            crf_softness: 0.0,
            mass: 0.0,
        }
    }
}

impl LambdaController {
    /// Create new lambda controller at λ = 0 (pure state A)
    pub fn new() -> Self {
        Self {
            lambda: 0.0,
            dlambda_dt: 0.0,
            lambda_exponent: 1,
            interaction_lambdas: InteractionLambdas::default(),
        }
    }

    /// Set initial lambda value
    pub fn with_lambda(mut self, lambda: f64) -> Self {
        self.lambda = lambda.clamp(0.0, 1.0);
        self.update_interaction_lambdas();
        self
    }

    /// Set lambda change rate for slow growth
    pub fn with_dlambda_dt(mut self, dlambda_dt: f64) -> Self {
        self.dlambda_dt = dlambda_dt;
        self
    }

    /// Set lambda exponent for λⁿ coupling
    pub fn with_exponent(mut self, n: i32) -> Self {
        self.lambda_exponent = n;
        self.update_interaction_lambdas();
        self
    }

    /// Update lambda by timestep (slow growth)
    ///
    /// Direct translation from slow_growth.cc
    pub fn update(&mut self, dt: f64) {
        if self.dlambda_dt != 0.0 {
            self.lambda += self.dlambda_dt * dt;
            self.lambda = self.lambda.clamp(0.0, 1.0);
            self.update_interaction_lambdas();
        }
    }

    /// Update all interaction lambdas based on main lambda
    fn update_interaction_lambdas(&mut self) {
        let lambda_n = if self.lambda_exponent == 1 {
            self.lambda
        } else {
            self.lambda.powi(self.lambda_exponent)
        };

        // Default: all interactions follow main lambda
        self.interaction_lambdas.bond = lambda_n;
        self.interaction_lambdas.angle = lambda_n;
        self.interaction_lambdas.dihedral = lambda_n;
        self.interaction_lambdas.improper = lambda_n;
        self.interaction_lambdas.lj = lambda_n;
        self.interaction_lambdas.lj_softness = self.lambda;  // Softness uses λ, not λⁿ
        self.interaction_lambdas.crf = lambda_n;
        self.interaction_lambdas.crf_softness = self.lambda;
        self.interaction_lambdas.mass = lambda_n;
    }

    /// Get lambda derivative (dλⁿ/dλ)
    pub fn lambda_derivative(&self) -> f64 {
        if self.lambda_exponent == 1 {
            1.0
        } else {
            (self.lambda_exponent as f64) * self.lambda.powi(self.lambda_exponent - 1)
        }
    }
}

impl Default for LambdaController {
    fn default() -> Self {
        Self::new()
    }
}

/// Perturbed atom with dual-state (A/B) parameters
///
/// Direct translation from perturbed_atom.h
#[derive(Debug, Clone)]
pub struct PerturbedAtom {
    /// Atom index in topology
    pub atom_index: usize,

    /// State A parameters
    pub a_iac: usize,        // Atom type
    pub a_mass: f64,         // Mass (amu)
    pub a_charge: f64,       // Charge (e)

    /// State B parameters
    pub b_iac: usize,
    pub b_mass: f64,
    pub b_charge: f64,

    /// Soft-core parameters
    pub lj_softcore: f64,    // α_LJ for Lennard-Jones
    pub crf_softcore: f64,   // α_CRF for electrostatics
}

impl PerturbedAtom {
    /// Create new perturbed atom
    pub fn new(atom_index: usize) -> Self {
        Self {
            atom_index,
            a_iac: 0,
            a_mass: 1.0,
            a_charge: 0.0,
            b_iac: 0,
            b_mass: 1.0,
            b_charge: 0.0,
            lj_softcore: 0.0,
            crf_softcore: 0.0,
        }
    }

    /// Set state A parameters
    pub fn with_state_a(mut self, iac: usize, mass: f64, charge: f64) -> Self {
        self.a_iac = iac;
        self.a_mass = mass;
        self.a_charge = charge;
        self
    }

    /// Set state B parameters
    pub fn with_state_b(mut self, iac: usize, mass: f64, charge: f64) -> Self {
        self.b_iac = iac;
        self.b_mass = mass;
        self.b_charge = charge;
        self
    }

    /// Set soft-core parameters
    pub fn with_softcore(mut self, lj: f64, crf: f64) -> Self {
        self.lj_softcore = lj;
        self.crf_softcore = crf;
        self
    }

    /// Get interpolated mass at given lambda
    pub fn mass_at_lambda(&self, lambda: f64) -> f64 {
        (1.0 - lambda) * self.a_mass + lambda * self.b_mass
    }

    /// Get interpolated charge at given lambda
    pub fn charge_at_lambda(&self, lambda: f64) -> f64 {
        (1.0 - lambda) * self.a_charge + lambda * self.b_charge
    }
}

/// Perturbed bonded interaction (bond, angle, dihedral)
#[derive(Debug, Clone)]
pub struct PerturbedBondedTerm {
    /// Atom indices
    pub atoms: Vec<usize>,

    /// State A parameters
    pub a_force_constant: f64,
    pub a_equilibrium: f64,

    /// State B parameters
    pub b_force_constant: f64,
    pub b_equilibrium: f64,
}

impl PerturbedBondedTerm {
    /// Get interpolated parameters at lambda
    pub fn params_at_lambda(&self, lambda: f64) -> (f64, f64) {
        let k = (1.0 - lambda) * self.a_force_constant + lambda * self.b_force_constant;
        let eq = (1.0 - lambda) * self.a_equilibrium + lambda * self.b_equilibrium;
        (k, eq)
    }

    /// Calculate lambda derivative of energy
    ///
    /// For harmonic potential: V = 0.5 * k * (x - x₀)²
    /// dV/dλ = 0.5 * (k_B - k_A) * (x - x₀)² - k(λ) * (x - x₀) * (x₀_B - x₀_A)
    pub fn lambda_derivative(&self, lambda: f64, current_value: f64) -> f64 {
        let (k_lambda, eq_lambda) = self.params_at_lambda(lambda);
        let delta = current_value - eq_lambda;

        let dk_dlambda = self.b_force_constant - self.a_force_constant;
        let deq_dlambda = self.b_equilibrium - self.a_equilibrium;

        0.5 * dk_dlambda * delta * delta - k_lambda * delta * deq_dlambda
    }
}

/// Soft-core potential parameters
///
/// Prevents singularities when atoms appear/disappear
#[derive(Debug, Clone)]
pub struct SoftCoreParameters {
    /// Alpha parameter for LJ soft-core
    pub alpha_lj: f64,
    /// Alpha parameter for electrostatic soft-core
    pub alpha_crf: f64,
    /// Soft-core power (usually 2)
    pub n_soft: i32,
}

impl Default for SoftCoreParameters {
    fn default() -> Self {
        Self {
            alpha_lj: 0.5,   // GROMOS default
            alpha_crf: 0.5,
            n_soft: 2,
        }
    }
}

impl SoftCoreParameters {
    /// Calculate soft-core distance for LJ interaction
    ///
    /// r_soft² = r² + α_LJ * (1-λ)² * σ⁶
    pub fn softcore_distance_lj(&self, r_sq: f64, lambda: f64, sigma: f64) -> f64 {
        let lambda_term = (1.0 - lambda).powi(self.n_soft);
        let sigma6 = sigma.powi(6);
        (r_sq + self.alpha_lj * lambda_term * sigma6).sqrt()
    }

    /// Calculate soft-core distance for electrostatics
    ///
    /// r_soft = √(r² + α_CRF * (1-λ)²)
    pub fn softcore_distance_crf(&self, r_sq: f64, lambda: f64) -> f64 {
        let lambda_term = (1.0 - lambda).powi(self.n_soft);
        (r_sq + self.alpha_crf * lambda_term).sqrt()
    }
}

/// Free energy derivative accumulator
///
/// Tracks dH/dλ for thermodynamic integration
#[derive(Debug, Clone, Default)]
pub struct FreeEnergyDerivatives {
    pub bond: f64,
    pub angle: f64,
    pub dihedral: f64,
    pub improper: f64,
    pub lj: f64,
    pub crf: f64,
    pub restraints: f64,
}

impl FreeEnergyDerivatives {
    pub fn new() -> Self {
        Self::default()
    }

    /// Get total dH/dλ
    pub fn total(&self) -> f64 {
        self.bond + self.angle + self.dihedral + self.improper +
        self.lj + self.crf + self.restraints
    }

    /// Reset all derivatives to zero
    pub fn clear(&mut self) {
        *self = Self::default();
    }
}

/// Perturbation topology
///
/// Contains all perturbed atoms and interactions
#[derive(Debug, Clone, Default)]
pub struct PerturbedTopology {
    pub perturbed_atoms: Vec<PerturbedAtom>,
    pub perturbed_bonds: Vec<PerturbedBondedTerm>,
    pub perturbed_angles: Vec<PerturbedBondedTerm>,
    pub perturbed_dihedrals: Vec<PerturbedBondedTerm>,
    pub soft_core_params: SoftCoreParameters,
}

impl PerturbedTopology {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_perturbed_atom(&mut self, atom: PerturbedAtom) {
        self.perturbed_atoms.push(atom);
    }

    pub fn add_perturbed_bond(&mut self, bond: PerturbedBondedTerm) {
        self.perturbed_bonds.push(bond);
    }

    pub fn add_perturbed_angle(&mut self, angle: PerturbedBondedTerm) {
        self.perturbed_angles.push(angle);
    }

    pub fn add_perturbed_dihedral(&mut self, dihedral: PerturbedBondedTerm) {
        self.perturbed_dihedrals.push(dihedral);
    }

    /// Update topology parameters for current lambda
    pub fn update_for_lambda(&self, topo: &mut Topology, lambda: f64) {
        // Update perturbed atom properties
        for perturbed in &self.perturbed_atoms {
            let idx = perturbed.atom_index;
            if idx < topo.mass.len() {
                topo.mass[idx] = perturbed.mass_at_lambda(lambda);
                topo.charge[idx] = perturbed.charge_at_lambda(lambda);
                topo.compute_inverse_masses();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lambda_controller() {
        let mut controller = LambdaController::new()
            .with_lambda(0.0)
            .with_dlambda_dt(0.1);  // 0.1 per ps

        assert_eq!(controller.lambda, 0.0);

        // Update by 1 ps
        controller.update(1.0);
        assert!((controller.lambda - 0.1).abs() < 1e-10);

        // Update by 5 ps
        controller.update(5.0);
        assert!((controller.lambda - 0.6).abs() < 1e-10);

        // Should clamp at 1.0
        controller.update(10.0);
        assert_eq!(controller.lambda, 1.0);
    }

    #[test]
    fn test_lambda_exponent() {
        let controller = LambdaController::new()
            .with_lambda(0.5)
            .with_exponent(2);

        assert_eq!(controller.lambda, 0.5);
        assert_eq!(controller.interaction_lambdas.bond, 0.25);  // 0.5²
        assert_eq!(controller.lambda_derivative(), 1.0);  // d(λ²)/dλ = 2λ = 1.0
    }

    #[test]
    fn test_perturbed_atom() {
        let atom = PerturbedAtom::new(0)
            .with_state_a(1, 12.0, 0.5)
            .with_state_b(2, 16.0, -0.5);

        // At λ = 0: pure state A
        assert_eq!(atom.mass_at_lambda(0.0), 12.0);
        assert_eq!(atom.charge_at_lambda(0.0), 0.5);

        // At λ = 1: pure state B
        assert_eq!(atom.mass_at_lambda(1.0), 16.0);
        assert_eq!(atom.charge_at_lambda(1.0), -0.5);

        // At λ = 0.5: average
        assert_eq!(atom.mass_at_lambda(0.5), 14.0);
        assert_eq!(atom.charge_at_lambda(0.5), 0.0);
    }

    #[test]
    fn test_soft_core_distance() {
        let params = SoftCoreParameters::default();

        // Without soft-core (λ = 1), should be normal distance
        let r_sq = 1.0;
        let r_soft = params.softcore_distance_lj(r_sq, 1.0, 1.0);
        assert!((r_soft - 1.0).abs() < 1e-10);

        // With soft-core (λ = 0), distance is increased
        let r_soft = params.softcore_distance_lj(r_sq, 0.0, 1.0);
        assert!(r_soft > 1.0);
    }

    #[test]
    fn test_perturbed_bonded_term() {
        let bond = PerturbedBondedTerm {
            atoms: vec![0, 1],
            a_force_constant: 100.0,
            a_equilibrium: 0.1,
            b_force_constant: 200.0,
            b_equilibrium: 0.15,
        };

        // At λ = 0.5
        let (k, r0) = bond.params_at_lambda(0.5);
        assert_eq!(k, 150.0);
        assert_eq!(r0, 0.125);
    }

    /// Mini-test 1: Lambda derivative for thermodynamic integration
    /// Tests dV/dλ calculation which is crucial for TI
    #[test]
    fn test_bonded_lambda_derivative() {
        let bond = PerturbedBondedTerm {
            atoms: vec![0, 1],
            a_force_constant: 100.0,  // State A: weak bond
            a_equilibrium: 0.1,
            b_force_constant: 200.0,  // State B: strong bond
            b_equilibrium: 0.15,
        };

        // Test at current bond length = 0.12 nm
        let current_length = 0.12;

        // At λ = 0 (pure state A)
        let deriv_0 = bond.lambda_derivative(0.0, current_length);

        // At λ = 1 (pure state B)
        let deriv_1 = bond.lambda_derivative(1.0, current_length);

        // At λ = 0.5 (midpoint)
        let deriv_05 = bond.lambda_derivative(0.5, current_length);

        // Derivatives should be finite (no NaN or Inf)
        assert!(deriv_0.is_finite(), "dV/dλ at λ=0 should be finite");
        assert!(deriv_05.is_finite(), "dV/dλ at λ=0.5 should be finite");
        assert!(deriv_1.is_finite(), "dV/dλ at λ=1 should be finite");

        // Note: Sign depends on whether bond is compressed or stretched
        // For current_length=0.12 between eq_A=0.1 and eq_B=0.15,
        // derivative sign varies with lambda

        println!("FEP Mini-test 1: dV/dλ at λ=0: {:.4}, λ=0.5: {:.4}, λ=1: {:.4}",
                 deriv_0, deriv_05, deriv_1);
    }

    /// Mini-test 2: Soft-core potential prevents singularities
    /// Tests that soft-core smoothly handles particle creation/deletion
    #[test]
    fn test_softcore_prevents_singularity() {
        let params = SoftCoreParameters {
            alpha_lj: 0.5,
            alpha_crf: 0.5,
            n_soft: 2,
        };

        let sigma = 0.3;  // nm

        // Test very close distance (would be singularity without soft-core)
        let r_sq_close = 0.01;  // 0.1 nm

        // At λ = 0 (particle appearing), soft-core is active
        let r_soft_lambda0 = params.softcore_distance_lj(r_sq_close, 0.0, sigma);
        // Should be larger than actual distance
        let r_actual = r_sq_close.sqrt();
        assert!(r_soft_lambda0 > r_actual, "Soft-core should increase distance at λ=0: {} vs {}",
                r_soft_lambda0, r_actual);

        // At λ = 0.5 (halfway), soft-core is partially active
        let r_soft_lambda05 = params.softcore_distance_lj(r_sq_close, 0.5, sigma);
        assert!(r_soft_lambda05 > r_sq_close.sqrt());
        assert!(r_soft_lambda05 < r_soft_lambda0, "Soft-core effect should decrease with λ");

        // At λ = 1 (particle fully present), soft-core is inactive
        let r_soft_lambda1 = params.softcore_distance_lj(r_sq_close, 1.0, sigma);
        assert!((r_soft_lambda1 - r_sq_close.sqrt()).abs() < 1e-6,
                "No soft-core at λ=1");

        println!("FEP Mini-test 2: Soft-core r_eff: λ=0: {:.4}, λ=0.5: {:.4}, λ=1: {:.4}",
                 r_soft_lambda0, r_soft_lambda05, r_soft_lambda1);
    }

    /// Mini-test 3: Charge transformation (common FEP use case)
    /// Tests smooth charge interpolation for electrostatic FEP
    #[test]
    fn test_charge_transformation() {
        // Transform neutral atom to charged atom
        let atom = PerturbedAtom::new(0)
            .with_state_a(1, 16.0, 0.0)   // State A: neutral O
            .with_state_b(1, 16.0, -0.8); // State B: charged O⁻

        // Test charge interpolation at multiple λ values
        let lambda_values = [0.0, 0.25, 0.5, 0.75, 1.0];
        let expected_charges = [0.0, -0.2, -0.4, -0.6, -0.8];

        for (lambda, expected) in lambda_values.iter().zip(expected_charges.iter()) {
            let charge = atom.charge_at_lambda(*lambda);
            assert!((charge - expected).abs() < 1e-10,
                    "Charge at λ={} should be {}, got {}", lambda, expected, charge);
        }

        println!("FEP Mini-test 3: Charge transformation test passed");
    }

    /// Mini-test 4: Multi-window FEP consistency
    /// Tests that energy is consistent across lambda windows
    #[test]
    fn test_multi_window_consistency() {
        let bond = PerturbedBondedTerm {
            atoms: vec![0, 1],
            a_force_constant: 100.0,
            a_equilibrium: 0.1,
            b_force_constant: 200.0,
            b_equilibrium: 0.15,
        };

        let current_length = 0.12;

        // Simulate 5-window FEP
        let lambda_windows = [0.0, 0.25, 0.5, 0.75, 1.0];
        let mut derivatives = Vec::new();

        for lambda in &lambda_windows {
            let deriv = bond.lambda_derivative(*lambda, current_length);
            derivatives.push(deriv);
        }

        // Check that derivatives are reasonable (no NaN, no huge jumps)
        for (i, deriv) in derivatives.iter().enumerate() {
            assert!(deriv.is_finite(), "Derivative at window {} should be finite", i);
        }

        // Check smoothness: no derivative should be 10x larger than adjacent
        for i in 0..derivatives.len()-1 {
            let ratio = (derivatives[i+1] / derivatives[i]).abs();
            assert!(ratio < 10.0 && ratio > 0.1,
                    "Large jump between windows {} and {}: ratio={:.2}",
                    i, i+1, ratio);
        }

        println!("FEP Mini-test 4: Multi-window test passed");
        println!("  dV/dλ values: {:?}", derivatives);
    }

    /// Mini-test 5: Slow growth simulation
    /// Tests continuous lambda evolution (slow growth method)
    #[test]
    fn test_slow_growth_evolution() {
        let mut controller = LambdaController::new()
            .with_lambda(0.0)
            .with_dlambda_dt(0.01);  // 0.01 per ps = transform over 100 ps

        let dt = 0.002;  // 2 fs timestep
        let n_steps = 50000;  // 100 ps

        let mut lambda_history = Vec::new();

        // Simulate slow growth
        for _step in 0..n_steps {
            lambda_history.push(controller.lambda);
            controller.update(dt);
        }

        // Check final lambda
        assert!((controller.lambda - 1.0).abs() < 1e-6,
                "Should reach λ=1 after 100 ps");

        // Check monotonicity (lambda should always increase)
        for i in 0..lambda_history.len()-1 {
            assert!(lambda_history[i+1] >= lambda_history[i],
                    "Lambda should be monotonically increasing");
        }

        // Check smoothness (no jumps > 0.05)
        for i in 0..lambda_history.len()-1 {
            let delta = lambda_history[i+1] - lambda_history[i];
            assert!(delta <= 0.05, "Lambda change too large: {}", delta);
        }

        println!("FEP Mini-test 5: Slow growth test passed");
        println!("  Final λ: {:.6}, steps: {}", controller.lambda, n_steps);
    }

    /// Mini-test 6: Soft-core electrostatics
    /// Tests soft-core for Coulomb interactions
    #[test]
    fn test_softcore_electrostatics() {
        let params = SoftCoreParameters::default();

        let r_sq = 0.01;  // Very close distance

        // Test soft-core for electrostatics at different λ
        let lambdas = [0.0, 0.3, 0.7, 1.0];
        let mut r_soft_values = Vec::new();

        for lambda in &lambdas {
            let r_soft = params.softcore_distance_crf(r_sq, *lambda);
            r_soft_values.push(r_soft);

            // At λ=0, soft-core should be maximally active
            if *lambda == 0.0 {
                assert!(r_soft > 0.15, "Electrostatic soft-core should be active at λ=0");
            }

            // At λ=1, soft-core should be minimal
            if *lambda == 1.0 {
                assert!((r_soft - r_sq.sqrt()).abs() < 1e-6,
                        "No electrostatic soft-core at λ=1");
            }
        }

        // Check monotonicity: r_soft should decrease as λ increases
        for i in 0..r_soft_values.len()-1 {
            assert!(r_soft_values[i] >= r_soft_values[i+1],
                    "Soft-core distance should decrease with increasing λ");
        }

        println!("FEP Mini-test 6: Electrostatic soft-core test passed");
        println!("  r_soft values: {:?}", r_soft_values);
    }
}
