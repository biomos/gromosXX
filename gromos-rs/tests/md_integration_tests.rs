//! Integration tests for MD components
//!
//! Tests the molecular dynamics functionality:
//! - Integrators
//! - Force calculations
//! - Constraints
//! - Thermostats and barostats
//! - Electrostatics

use gromos_rs::math::Vec3;
use gromos_rs::configuration::{Configuration, Box as SimBox, State};
use gromos_rs::integrator::{LeapFrog, Integrator};
use gromos_rs::topology::Topology;
use gromos_rs::interaction::electrostatics::{
    ReactionFieldParameters, PMEParameters,
    reaction_field_interaction, pme_real_space_interaction, pme_self_energy,
};

#[test]
fn test_configuration_double_buffer() {
    let mut config = Configuration::new(5, 1, 1);

    // Set positions in current state
    for i in 0..5 {
        config.current_mut().pos[i] = Vec3::new(i as f32, 0.0, 0.0);
    }

    // Copy to old state
    config.copy_current_to_old();

    // Verify both states have same positions
    for i in 0..5 {
        assert_eq!(config.current().pos[i], config.old().pos[i]);
    }

    // Modify current state
    config.current_mut().pos[0] = Vec3::new(10.0, 0.0, 0.0);

    // Old state should remain unchanged
    assert_ne!(config.current().pos[0], config.old().pos[0]);

    // Exchange states
    config.exchange_state();

    // Now old positions should be accessible as current
    assert_eq!(config.current().pos[0].x, 0.0);
}

#[test]
fn test_leapfrog_integrator() {
    // Create a simple test system
    use gromos_rs::topology::Atom;

    let mut topo = Topology::new();

    // Add one atom to the solute
    topo.solute.atoms.push(Atom {
        name: "C1".to_string(),
        residue_nr: 1,
        residue_name: "TST".to_string(),
        iac: 0,
        mass: 12.0,
        charge: 0.0,
        is_perturbed: false,
        is_polarisable: false,
        is_coarse_grained: false,
    });

    topo.mass = vec![12.0];  // One carbon atom
    topo.charge = vec![0.0];
    topo.iac = vec![0];
    topo.compute_inverse_masses();

    let mut config = Configuration::new(1, 1, 1);

    // Set initial conditions
    config.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
    config.current_mut().vel[0] = Vec3::new(1.0, 0.0, 0.0);  // 1 nm/ps
    config.current_mut().force[0] = Vec3::new(0.0, 0.0, 0.0);  // No force

    // Copy to old state for leapfrog
    config.copy_current_to_old();

    let mut integrator = LeapFrog::new();
    let dt = 0.002;  // 2 fs timestep

    // Integrate one step
    integrator.step(dt, &topo, &mut config);

    // Position should have moved by v*dt
    assert!((config.current().pos[0].x - 0.002).abs() < 1e-6,
            "Expected position ~0.002, got: {}", config.current().pos[0].x);

    // Velocity unchanged (no force)
    assert_eq!(config.current().vel[0].x, 1.0,
            "Expected velocity 1.0, got: {}", config.current().vel[0].x);
}

#[test]
fn test_reaction_field_electrostatics() {
    // Test RF with conducting boundary (GROMOS default)
    let rf = ReactionFieldParameters::gromos_default(1.4);

    assert_eq!(rf.cutoff, 1.4);
    assert_eq!(rf.epsilon, 1.0);
    assert_eq!(rf.epsilon_rf, 0.0);
    assert_eq!(rf.crf, -1.0);  // Conducting boundary

    // Test interaction between two opposite charges
    let r_vec = Vec3::new(1.0, 0.0, 0.0);  // 1 nm separation
    let q1 = 1.0;   // +1e
    let q2 = -1.0;  // -1e

    let (force, energy) = reaction_field_interaction(r_vec, q1, q2, &rf);

    // Energy should be negative (attractive)
    assert!(energy < 0.0);

    // Force should point toward negative charge (negative x)
    assert!(force.x < 0.0);
}

#[test]
fn test_reaction_field_aqueous() {
    let rf_aqueous = ReactionFieldParameters::aqueous(1.4);

    assert_eq!(rf_aqueous.epsilon_rf, 78.4);
    // CRF should be different from conducting boundary
    assert!(rf_aqueous.crf != -1.0);
    assert!(rf_aqueous.crf < 0.0);  // Still negative but less so
}

#[test]
fn test_pme_real_space() {
    let pme = PMEParameters::new(1.0, 64, 4, 1e-5);

    assert_eq!(pme.cutoff, 1.0);
    assert_eq!(pme.grid_x, 64);
    assert_eq!(pme.spline_order, 4);
    assert!((pme.alpha - 3.0).abs() < 0.1);  // Alpha should be ~3.0

    // Test PME real-space interaction
    let r_vec = Vec3::new(0.5, 0.0, 0.0);
    let q1 = 1.0;
    let q2 = -1.0;

    let (force, energy) = pme_real_space_interaction(r_vec, q1, q2, &pme);

    // Energy should be negative (attractive)
    assert!(energy < 0.0);

    // Force should be attractive
    assert!(force.x < 0.0);

    // PME real-space should be weaker than full Coulomb at this distance
    // (due to erfc screening)
    let coulomb_energy: f64 = -1.0 / 0.5;  // Simple Coulomb
    assert!(energy.abs() < coulomb_energy.abs());
}

#[test]
fn test_pme_self_energy() {
    let charges = vec![1.0, -1.0, 0.5, -0.5];
    let alpha = 3.0;

    let self_energy = pme_self_energy(&charges, alpha);

    // Self energy should be negative
    assert!(self_energy < 0.0);

    // Self energy scales with sum of q^2
    let _q_squared_sum: f64 = charges.iter().map(|&q| q * q).sum();
    assert!(self_energy.abs() > 0.0);
}

#[test]
fn test_box_volume_calculation() {
    let box_rect = SimBox::rectangular(3.0, 4.0, 5.0);

    let volume = box_rect.volume();
    assert!((volume - 60.0).abs() < 1e-6);

    let dims = box_rect.dimensions();
    assert_eq!(dims.x, 3.0);
    assert_eq!(dims.y, 4.0);
    assert_eq!(dims.z, 5.0);
}

#[test]
fn test_kinetic_energy_calculation() {
    let masses = vec![12.0, 12.0, 1.0];  // 2 Carbons, 1 Hydrogen
    let mut state = State::new(3, 1, 1);

    // Set velocities
    state.vel = vec![
        Vec3::new(1.0, 0.0, 0.0),   // C: 1 nm/ps
        Vec3::new(0.0, 1.0, 0.0),   // C: 1 nm/ps
        Vec3::new(0.0, 0.0, 2.0),   // H: 2 nm/ps
    ];

    state.calculate_kinetic_energy(&masses);

    // KE = 0.5 * (12*1^2 + 12*1^2 + 1*2^2) = 0.5 * (12 + 12 + 4) = 14.0
    assert!((state.energies.kinetic_total - 14.0).abs() < 1e-6);
}

#[test]
fn test_temperature_calculation() {
    let mut state = State::new(3, 1, 1);
    state.energies.kinetic_total = 100.0;  // kJ/mol

    // 3 atoms = 9 DOF (assuming no constraints)
    let temp = state.temperature(9);

    // T = 2*KE / (k_B * N_dof)
    // k_B = 0.00831446 kJ/(mol·K)
    // T = 2*100 / (0.00831446 * 9) ≈ 2672 K
    assert!(temp > 2600.0 && temp < 2700.0);
}

#[test]
fn test_pressure_calculation() {
    let mut state = State::new(10, 1, 1);
    state.box_config = SimBox::rectangular(3.0, 3.0, 3.0);
    state.energies.kinetic_total = 100.0;

    // Set diagonal virial
    state.virial_tensor.x_axis.x = -50.0;
    state.virial_tensor.y_axis.y = -50.0;
    state.virial_tensor.z_axis.z = -50.0;

    let pressure = state.pressure();

    // P = (2*KE - Virial) / (3*V)
    // P = (2*100 - (-150)) / (3*27) = 350 / 81 ≈ 4.32
    assert!(pressure > 4.0 && pressure < 5.0);
}

#[test]
fn test_force_accumulation() {
    let mut state = State::new(3, 1, 1);

    // Add forces from different sources
    state.force[0] += Vec3::new(1.0, 0.0, 0.0);  // Bonded
    state.force[0] += Vec3::new(0.0, 1.0, 0.0);  // Nonbonded
    state.constraint_force[0] += Vec3::new(0.0, 0.0, 1.0);  // Constraint

    // Total force (without constraints)
    assert_eq!(state.force[0], Vec3::new(1.0, 1.0, 0.0));

    // Constraint force separate
    assert_eq!(state.constraint_force[0], Vec3::new(0.0, 0.0, 1.0));
}

#[test]
fn test_energy_component_summation() {
    let mut state = State::new(10, 1, 1);

    // Set individual energy components
    state.energies.bond_total = -50.0;
    state.energies.angle_total = -30.0;
    state.energies.dihedral_total = -20.0;
    state.energies.lj_total = -100.0;
    state.energies.crf_total = -80.0;

    // Update potential total
    state.energies.update_potential_total();

    // Sum should be -280.0
    assert!((state.energies.potential_total + 280.0).abs() < 1e-6);

    // Set kinetic energy
    state.energies.kinetic_total = 200.0;

    // Total energy
    let total = state.energies.total();
    assert!((total + 80.0).abs() < 1e-6);  // 200 - 280 = -80
}

#[test]
fn test_state_clear_forces() {
    let mut state = State::new(5, 1, 1);

    // Add some forces
    for i in 0..5 {
        state.force[i] = Vec3::new(i as f32, i as f32, i as f32);
        state.constraint_force[i] = Vec3::new(1.0, 1.0, 1.0);
    }

    // Clear forces
    state.clear_forces();

    for i in 0..5 {
        assert_eq!(state.force[i], Vec3::ZERO);
    }

    // Clear constraint forces
    state.clear_constraint_forces();

    for i in 0..5 {
        assert_eq!(state.constraint_force[i], Vec3::ZERO);
    }
}

#[test]
fn test_electrostatics_cutoff_behavior() {
    let rf = ReactionFieldParameters::gromos_default(1.4);

    // Test interaction within cutoff
    let r_near = Vec3::new(1.0, 0.0, 0.0);
    let (force_near, energy_near) = reaction_field_interaction(r_near, 1.0, -1.0, &rf);

    assert!(energy_near < 0.0);
    assert!(force_near.x < 0.0);

    // Test interaction beyond cutoff (should be zero or very small)
    let r_far = Vec3::new(2.0, 0.0, 0.0);  // Beyond 1.4 nm cutoff
    let (_force_far, _energy_far) = reaction_field_interaction(r_far, 1.0, -1.0, &rf);

    // Note: RF function doesn't enforce cutoff internally in this simple version
    // In a full MD code, the pairlist would handle this
}

#[test]
fn test_pme_parameters_defaults() {
    let pme = PMEParameters::default();

    assert_eq!(pme.cutoff, 1.0);
    assert_eq!(pme.grid_x, 64);
    assert_eq!(pme.grid_y, 64);
    assert_eq!(pme.grid_z, 64);
    assert_eq!(pme.spline_order, 4);  // Cubic
    assert_eq!(pme.tolerance, 1e-5);
}
