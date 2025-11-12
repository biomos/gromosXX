//! Interaction calculations (nonbonded, bonded, electrostatics, restraints)

pub mod nonbonded;
pub mod bonded;
pub mod electrostatics;
pub mod restraints;

// Re-export commonly used types and functions
pub use bonded::{
    ForceEnergy,
    ForceEnergyLambda,
    calculate_bond_forces_quartic,
    calculate_bond_forces_harmonic,
    calculate_angle_forces,
    calculate_dihedral_forces,
    calculate_improper_dihedral_forces,
    calculate_bonded_forces,
    calculate_perturbed_bond_forces,
    calculate_perturbed_angle_forces,
    calculate_perturbed_dihedral_forces,
};

pub use electrostatics::{
    ElectrostaticsMethod,
    ReactionFieldParameters,
    PMEParameters,
    reaction_field_interaction,
    pme_real_space_interaction,
    pme_reciprocal_space,
    pme_self_energy,
};

pub use restraints::{
    PositionRestraint,
    PositionRestraints,
    DistanceRestraint,
    DistanceRestraints,
    AngleRestraint,
    AngleRestraints,
    DihedralRestraint,
    DihedralRestraints,
};
