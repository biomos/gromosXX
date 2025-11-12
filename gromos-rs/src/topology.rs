//! Topology module - molecular structure and force field parameters
//!
//! This is a direct Rust translation of GROMOS topology structures from:
//! - md++/src/topology/topology.h
//! - md++/src/topology/solute.h
//! - md++/src/topology/solvent.h

use std::collections::HashSet;
use crate::math::Vec3;

/// Atom properties
#[derive(Debug, Clone)]
pub struct Atom {
    pub name: String,
    pub residue_nr: usize,
    pub residue_name: String,
    pub iac: usize,           // Integer atom code (atom type)
    pub mass: f64,
    pub charge: f64,
    pub is_perturbed: bool,
    pub is_polarisable: bool,
    pub is_coarse_grained: bool,
}

/// Two-body bonded term (bonds, harmonic constraints)
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    pub i: usize,
    pub j: usize,
    pub bond_type: usize,
}

/// Three-body bonded term (angles)
#[derive(Debug, Clone, Copy)]
pub struct Angle {
    pub i: usize,
    pub j: usize,  // Central atom
    pub k: usize,
    pub angle_type: usize,
}

/// Four-body bonded term (proper and improper dihedrals)
#[derive(Debug, Clone, Copy)]
pub struct Dihedral {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub dihedral_type: usize,
}

/// Cross-dihedral term (8 atoms)
#[derive(Debug, Clone, Copy)]
pub struct CrossDihedral {
    pub a: usize,
    pub b: usize,
    pub c: usize,
    pub d: usize,
    pub e: usize,
    pub f: usize,
    pub g: usize,
    pub h: usize,
    pub cross_dihedral_type: usize,
}

/// Perturbed bond (for FEP calculations)
#[derive(Debug, Clone, Copy)]
pub struct PerturbedBond {
    pub i: usize,
    pub j: usize,
    pub a_type: usize,  // State A bond type
    pub b_type: usize,  // State B bond type
}

/// Perturbed angle (for FEP calculations)
#[derive(Debug, Clone, Copy)]
pub struct PerturbedAngle {
    pub i: usize,
    pub j: usize,  // Central atom
    pub k: usize,
    pub a_type: usize,  // State A angle type
    pub b_type: usize,  // State B angle type
}

/// Perturbed dihedral (for FEP calculations)
#[derive(Debug, Clone, Copy)]
pub struct PerturbedDihedral {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub a_type: usize,  // State A dihedral type
    pub b_type: usize,  // State B dihedral type
}

/// Bond force field parameters (GROMOS format)
#[derive(Debug, Clone, Copy)]
pub struct BondParameters {
    pub k_quartic: f64,   // Quartic force constant
    pub k_harmonic: f64,  // Harmonic force constant
    pub r0: f64,          // Equilibrium bond length
}

/// Angle force field parameters (GROMOS format)
#[derive(Debug, Clone, Copy)]
pub struct AngleParameters {
    pub k_cosine: f64,   // Force constant (harmonic in cosine)
    pub k_harmonic: f64, // Force constant (harmonic in angle)
    pub theta0: f64,     // Equilibrium angle in radians
}

/// Dihedral force field parameters
#[derive(Debug, Clone, Copy)]
pub struct DihedralParameters {
    pub k: f64,      // Force constant
    pub pd: f64,     // Phase shift
    pub cospd: f64,  // Cos(phase)
    pub m: i32,      // Multiplicity
}

/// Improper dihedral force field parameters
#[derive(Debug, Clone, Copy)]
pub struct ImproperDihedralParameters {
    pub k: f64,   // Force constant
    pub q0: f64,  // Equilibrium improper dihedral angle
}

/// Lennard-Jones parameters
#[derive(Debug, Clone, Copy)]
pub struct LJParameters {
    pub c6: f64,    // C6 coefficient (attractive: -C6/r^6)
    pub c12: f64,   // C12 coefficient (repulsive: C12/r^12)
    pub cs6: f64,   // Softcore C6 (for free energy calculations)
    pub cs12: f64,  // Softcore C12
}

impl LJParameters {
    /// Create standard LJ parameters (without softcore)
    pub fn new(c6: f64, c12: f64) -> Self {
        Self {
            c6,
            c12,
            cs6: c6,
            cs12: c12,
        }
    }

    /// Calculate sigma (size parameter) from C6 and C12
    /// sigma = (C12/C6)^(1/6)
    pub fn sigma(&self) -> f64 {
        if self.c6 > 0.0 {
            (self.c12 / self.c6).powf(1.0 / 6.0)
        } else {
            0.0
        }
    }

    /// Calculate epsilon (well depth) from C6 and C12
    /// epsilon = C6^2 / (4 * C12)
    pub fn epsilon(&self) -> f64 {
        if self.c12 > 0.0 {
            self.c6 * self.c6 / (4.0 * self.c12)
        } else {
            0.0
        }
    }
}

impl Default for LJParameters {
    fn default() -> Self {
        Self {
            c6: 0.0,
            c12: 0.0,
            cs6: 0.0,
            cs12: 0.0,
        }
    }
}

/// Exclusion list for an atom (atoms that don't interact via nonbonded forces)
///
/// Uses HashSet for fast O(1) lookups. In GROMOS C++, this is a sorted vector
/// with binary search (O(log n)), but HashSet is more idiomatic in Rust.
pub type Exclusions = HashSet<usize>;

/// Chargegroup - atoms whose charge sum is zero (or small)
///
/// Chargegroups are used for cutoff optimizations:
/// - Distance calculated between chargegroup centers
/// - All atoms in both groups included if distance < cutoff
#[derive(Debug, Clone)]
pub struct ChargeGroup {
    pub atoms: Vec<usize>,
}

impl ChargeGroup {
    /// Calculate center of geometry for this chargegroup
    pub fn center_of_geometry(&self, positions: &[Vec3]) -> Vec3 {
        if self.atoms.is_empty() {
            return Vec3::ZERO;
        }

        let sum: Vec3 = self.atoms.iter()
            .map(|&i| positions[i])
            .sum();

        sum / self.atoms.len() as f32
    }
}

/// Solute molecule structure
#[derive(Debug, Clone)]
pub struct Solute {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub proper_dihedrals: Vec<Dihedral>,
    pub improper_dihedrals: Vec<Dihedral>,
    pub cross_dihedrals: Vec<CrossDihedral>,
}

impl Solute {
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            proper_dihedrals: Vec::new(),
            improper_dihedrals: Vec::new(),
            cross_dihedrals: Vec::new(),
        }
    }

    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }
}

/// Perturbed solute for FEP calculations
///
/// Contains dual-topology (A/B state) bonded terms for free energy perturbation
#[derive(Debug, Clone)]
pub struct PerturbedSolute {
    pub bonds: Vec<PerturbedBond>,
    pub angles: Vec<PerturbedAngle>,
    pub proper_dihedrals: Vec<PerturbedDihedral>,
    pub improper_dihedrals: Vec<PerturbedDihedral>,
}

impl PerturbedSolute {
    pub fn new() -> Self {
        Self {
            bonds: Vec::new(),
            angles: Vec::new(),
            proper_dihedrals: Vec::new(),
            improper_dihedrals: Vec::new(),
        }
    }
}

impl Default for PerturbedSolute {
    fn default() -> Self {
        Self::new()
    }
}

/// Solvent molecule structure (typically water)
#[derive(Debug, Clone)]
pub struct Solvent {
    pub name: String,
    pub atoms: Vec<Atom>,
    pub num_molecules: usize,
}

impl Solvent {
    pub fn new(name: String) -> Self {
        Self {
            name,
            atoms: Vec::new(),
            num_molecules: 0,
        }
    }

    pub fn atoms_per_molecule(&self) -> usize {
        self.atoms.len()
    }

    pub fn total_atoms(&self) -> usize {
        self.atoms.len() * self.num_molecules
    }
}

/// Main topology structure containing all molecular information
#[derive(Debug, Clone)]
pub struct Topology {
    // Solute and solvent
    pub solute: Solute,
    pub perturbed_solute: PerturbedSolute,  // FEP dual-topology terms
    pub solvents: Vec<Solvent>,

    // Per-atom properties (flat arrays for all atoms)
    pub iac: Vec<usize>,           // Integer atom codes (atom types)
    pub mass: Vec<f64>,            // Atomic masses
    pub inverse_mass: Vec<f64>,    // Precomputed 1/mass for efficiency
    pub charge: Vec<f64>,          // Atomic charges

    // Exclusions (atoms that don't interact via nonbonded)
    pub exclusions: Vec<Exclusions>,  // exclusions[i] = set of atoms excluded from i
    pub one_four_pairs: Vec<Vec<usize>>, // 1-4 pairs (special scaling)

    // Chargegroups
    pub chargegroups: Vec<ChargeGroup>,
    pub atom_to_chargegroup: Vec<usize>,  // atom -> chargegroup index

    // Temperature and pressure coupling groups
    pub temperature_groups: Vec<Vec<usize>>,  // Atoms in each T-coupling group
    pub pressure_groups: Vec<Vec<usize>>,     // Atoms in each P-coupling group
    pub atom_to_temperature_group: Vec<usize>,
    pub atom_to_pressure_group: Vec<usize>,

    // Energy groups (for energy monitoring)
    pub energy_groups: Vec<Vec<usize>>,
    pub atom_to_energy_group: Vec<usize>,

    // Molecule boundaries (for molecule-based operations)
    pub molecules: Vec<std::ops::Range<usize>>,  // molecule_start..molecule_end

    // Force field parameters
    pub bond_parameters: Vec<BondParameters>,
    pub angle_parameters: Vec<AngleParameters>,
    pub dihedral_parameters: Vec<DihedralParameters>,
    pub improper_dihedral_parameters: Vec<ImproperDihedralParameters>,
    pub lj_parameters: Vec<Vec<LJParameters>>,  // [type_i][type_j] matrix
}

impl Topology {
    pub fn new() -> Self {
        Self {
            solute: Solute::new(),
            perturbed_solute: PerturbedSolute::new(),
            solvents: Vec::new(),
            iac: Vec::new(),
            mass: Vec::new(),
            inverse_mass: Vec::new(),
            charge: Vec::new(),
            exclusions: Vec::new(),
            one_four_pairs: Vec::new(),
            chargegroups: Vec::new(),
            atom_to_chargegroup: Vec::new(),
            temperature_groups: Vec::new(),
            pressure_groups: Vec::new(),
            atom_to_temperature_group: Vec::new(),
            atom_to_pressure_group: Vec::new(),
            energy_groups: Vec::new(),
            atom_to_energy_group: Vec::new(),
            molecules: Vec::new(),
            bond_parameters: Vec::new(),
            angle_parameters: Vec::new(),
            dihedral_parameters: Vec::new(),
            improper_dihedral_parameters: Vec::new(),
            lj_parameters: Vec::new(),
        }
    }

    /// Total number of atoms in the system
    pub fn num_atoms(&self) -> usize {
        self.solute.num_atoms() + self.solvents.iter().map(|s| s.total_atoms()).sum::<usize>()
    }

    /// Number of atom types
    pub fn num_atom_types(&self) -> usize {
        self.lj_parameters.len()
    }

    /// Get LJ parameters for atom pair (i, j)
    #[inline]
    pub fn lj_parameter(&self, type_i: usize, type_j: usize) -> &LJParameters {
        &self.lj_parameters[type_i][type_j]
    }

    /// Check if atoms i and j are excluded from nonbonded interactions
    #[inline]
    pub fn is_excluded(&self, i: usize, j: usize) -> bool {
        self.exclusions[i].contains(&j) || self.exclusions[j].contains(&i)
    }

    /// Get chargegroup index for atom i
    #[inline]
    pub fn chargegroup(&self, i: usize) -> usize {
        self.atom_to_chargegroup[i]
    }

    /// Get temperature group index for atom i
    #[inline]
    pub fn temperature_group(&self, i: usize) -> usize {
        self.atom_to_temperature_group[i]
    }

    /// Get energy group index for atom i
    #[inline]
    pub fn energy_group(&self, i: usize) -> usize {
        self.atom_to_energy_group[i]
    }

    /// Build exclusions from bonded topology
    ///
    /// Excludes:
    /// - 1-2 bonded neighbors
    /// - 1-3 angle neighbors
    /// - Optionally 1-4 dihedral neighbors (or use special 1-4 scaling)
    pub fn build_exclusions(&mut self, exclude_14: bool) {
        let n_atoms = self.num_atoms();
        self.exclusions = vec![HashSet::new(); n_atoms];

        // Exclude bonded neighbors (1-2)
        for bond in &self.solute.bonds {
            self.exclusions[bond.i].insert(bond.j);
            self.exclusions[bond.j].insert(bond.i);
        }

        // Exclude angle neighbors (1-3)
        for angle in &self.solute.angles {
            self.exclusions[angle.i].insert(angle.k);
            self.exclusions[angle.k].insert(angle.i);
        }

        // Optionally exclude or mark 1-4 pairs
        if exclude_14 {
            for dihedral in &self.solute.proper_dihedrals {
                self.exclusions[dihedral.i].insert(dihedral.l);
                self.exclusions[dihedral.l].insert(dihedral.i);
            }
        } else {
            // Store as special 1-4 pairs
            self.one_four_pairs = vec![Vec::new(); n_atoms];
            for dihedral in &self.solute.proper_dihedrals {
                self.one_four_pairs[dihedral.i].push(dihedral.l);
                self.one_four_pairs[dihedral.l].push(dihedral.i);
            }
        }
    }

    /// Initialize LJ parameter matrix from combining rules
    ///
    /// GROMOS uses geometric mean combining rules:
    /// - C6_ij = sqrt(C6_ii * C6_jj)
    /// - C12_ij = sqrt(C12_ii * C12_jj)
    pub fn build_lj_matrix(&mut self, lj_types: &[LJParameters]) {
        let n_types = lj_types.len();
        self.lj_parameters = vec![vec![LJParameters::new(0.0, 0.0); n_types]; n_types];

        for i in 0..n_types {
            for j in 0..n_types {
                let c6 = (lj_types[i].c6 * lj_types[j].c6).sqrt();
                let c12 = (lj_types[i].c12 * lj_types[j].c12).sqrt();
                self.lj_parameters[i][j] = LJParameters::new(c6, c12);
            }
        }
    }

    /// Resize all per-atom arrays to match number of atoms
    pub fn resize_atom_arrays(&mut self) {
        let n_atoms = self.num_atoms();
        self.iac.resize(n_atoms, 0);
        self.mass.resize(n_atoms, 0.0);
        self.inverse_mass.resize(n_atoms, 0.0);
        self.charge.resize(n_atoms, 0.0);
        self.exclusions.resize(n_atoms, HashSet::new());
        self.one_four_pairs.resize(n_atoms, Vec::new());
        self.atom_to_chargegroup.resize(n_atoms, 0);
        self.atom_to_temperature_group.resize(n_atoms, 0);
        self.atom_to_pressure_group.resize(n_atoms, 0);
        self.atom_to_energy_group.resize(n_atoms, 0);
    }

    /// Compute inverse masses for all atoms (for efficient integration)
    pub fn compute_inverse_masses(&mut self) {
        self.inverse_mass = self.mass.iter()
            .map(|&m| if m > 0.0 { 1.0 / m } else { 0.0 })
            .collect();
    }
}

impl Default for Topology {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lj_parameters() {
        // Typical values for sp3 carbon
        let c6 = 0.0022175;  // kJ/mol nm^6
        let c12 = 0.0014297; // kJ/mol nm^12
        let lj = LJParameters::new(c6, c12);

        // Check sigma (size parameter)
        let sigma = lj.sigma();
        assert!((sigma - 0.375).abs() < 0.01, "Sigma should be ~0.375 nm");

        // Check epsilon (well depth)
        let epsilon = lj.epsilon();
        assert!((epsilon - 0.866).abs() < 0.1, "Epsilon should be ~0.866 kJ/mol");
    }

    #[test]
    fn test_topology_creation() {
        let mut topo = Topology::new();

        // Add a simple water molecule
        let mut water = Solvent::new("SOL".to_string());
        water.atoms.push(Atom {
            name: "OW".to_string(),
            residue_nr: 1,
            residue_name: "SOL".to_string(),
            iac: 0,
            mass: 15.9994,
            charge: -0.82,
            is_perturbed: false,
            is_polarisable: false,
            is_coarse_grained: false,
        });
        water.num_molecules = 100;

        topo.solvents.push(water);

        assert_eq!(topo.num_atoms(), 100);
    }

    #[test]
    fn test_exclusions() {
        let mut topo = Topology::new();
        topo.solute.bonds.push(Bond { i: 0, j: 1, bond_type: 0 });
        topo.solute.bonds.push(Bond { i: 1, j: 2, bond_type: 0 });
        topo.solute.angles.push(Angle { i: 0, j: 1, k: 2, angle_type: 0 });

        topo.resize_atom_arrays();
        topo.build_exclusions(false);

        // 0-1 should be excluded (bonded)
        assert!(topo.is_excluded(0, 1));

        // 0-2 should be excluded (1-3 neighbors)
        assert!(topo.is_excluded(0, 2));
    }

    #[test]
    fn test_lj_matrix() {
        let mut topo = Topology::new();

        // Two atom types
        let lj_types = vec![
            LJParameters::new(0.001, 0.0001),  // Type 0
            LJParameters::new(0.002, 0.0002),  // Type 1
        ];

        topo.build_lj_matrix(&lj_types);

        // Check diagonal elements
        assert!((topo.lj_parameter(0, 0).c6 - 0.001).abs() < 1e-10);
        assert!((topo.lj_parameter(1, 1).c6 - 0.002).abs() < 1e-10);

        // Check off-diagonal (geometric mean)
        let c6_01 = (0.001_f64 * 0.002_f64).sqrt();
        assert!((topo.lj_parameter(0, 1).c6 - c6_01).abs() < 1e-10);
    }
}
