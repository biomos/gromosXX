//! Atom selection system for GROMOS trajectory analysis
//!
//! Provides unified atom selection syntax matching GROMOS++ AtomSpecifier.
//!
//! # Selection Syntax
//!
//! - `all` - All atoms
//! - `1-100` - Atoms 1 through 100 (1-based)
//! - `1,5,10-20` - Atoms 1, 5, and 10 through 20
//! - `1:1-10` - Molecule 1, atoms 1-10
//! - `1:a` - All atoms in molecule 1
//! - `s:1-10` - Solvent molecules 1-10
//! - `r:1-5` - Residues 1-5
//! - `a:CA` - All atoms named CA
//! - `a:CA,CB,N` - All atoms named CA, CB, or N
//!
//! # Examples
//!
//! ```rust,ignore
//! use gromos_rs::selection::AtomSelection;
//! use gromos_rs::topology::Topology;
//!
//! // Create selection from string
//! let selection = AtomSelection::from_string("1:1-10", &topology)?;
//!
//! // Get selected atom indices (0-based)
//! let indices = selection.indices();
//!
//! // Iterate over selected atoms
//! for idx in selection.iter() {
//!     println!("Atom {}: {}", idx, topology.atoms[idx].name);
//! }
//! ```

use crate::topology::Topology;
use std::collections::HashSet;

/// Atom selection error
#[derive(Debug, Clone)]
pub enum SelectionError {
    ParseError(String),
    InvalidRange(String),
    InvalidMolecule(usize),
    InvalidResidue(usize),
    InvalidAtomName(String),
    EmptySelection,
}

impl std::fmt::Display for SelectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SelectionError::ParseError(msg) => write!(f, "Parse error: {}", msg),
            SelectionError::InvalidRange(msg) => write!(f, "Invalid range: {}", msg),
            SelectionError::InvalidMolecule(m) => write!(f, "Invalid molecule: {}", m),
            SelectionError::InvalidResidue(r) => write!(f, "Invalid residue: {}", r),
            SelectionError::InvalidAtomName(name) => write!(f, "Invalid atom name: {}", name),
            SelectionError::EmptySelection => write!(f, "Empty selection"),
        }
    }
}

impl std::error::Error for SelectionError {}

/// Atom selection
#[derive(Debug, Clone)]
pub struct AtomSelection {
    /// Selected atom indices (0-based)
    indices: Vec<usize>,
    /// Total number of atoms in system
    n_atoms: usize,
}

impl AtomSelection {
    /// Create selection from all atoms
    pub fn all(n_atoms: usize) -> Self {
        Self {
            indices: (0..n_atoms).collect(),
            n_atoms,
        }
    }

    /// Create selection from atom indices (0-based)
    pub fn from_indices(indices: Vec<usize>, n_atoms: usize) -> Result<Self, SelectionError> {
        // Validate indices
        for &idx in &indices {
            if idx >= n_atoms {
                return Err(SelectionError::InvalidRange(format!(
                    "Atom index {} >= {} total atoms",
                    idx, n_atoms
                )));
            }
        }

        if indices.is_empty() {
            return Err(SelectionError::EmptySelection);
        }

        Ok(Self { indices, n_atoms })
    }

    /// Parse selection string
    ///
    /// Supported formats:
    /// - `all` - All atoms
    /// - `1-100` - Atoms 1-100 (1-based)
    /// - `1,5,10-20` - Individual atoms and ranges
    /// - `1:1-10` - Molecule 1, atoms 1-10
    /// - `1:a` - All atoms in molecule 1
    /// - `s:1-10` - Solvent molecules 1-10
    /// - `r:1-5` - Residues 1-5
    /// - `a:CA` - All atoms named CA
    pub fn from_string(spec: &str, topology: &Topology) -> Result<Self, SelectionError> {
        let spec = spec.trim();

        if spec == "all" {
            return Ok(Self::all(topology.num_atoms()));
        }

        // Check for prefix (molecule:, solvent:, residue:, atom:)
        if spec.contains(':') {
            Self::parse_prefixed(spec, topology)
        } else {
            // Simple atom list
            Self::parse_atom_list(spec, topology.num_atoms())
        }
    }

    /// Parse simple atom list: "1-10", "1,5,10-20"
    fn parse_atom_list(spec: &str, n_atoms: usize) -> Result<Self, SelectionError> {
        let mut indices = HashSet::new();

        for part in spec.split(',') {
            let part = part.trim();

            if part.contains('-') {
                // Range: "1-10"
                let range: Vec<&str> = part.split('-').collect();
                if range.len() != 2 {
                    return Err(SelectionError::ParseError(format!("Invalid range: {}", part)));
                }

                let start: usize = range[0]
                    .trim()
                    .parse()
                    .map_err(|_| SelectionError::ParseError(format!("Invalid start: {}", range[0])))?;
                let end: usize = range[1]
                    .trim()
                    .parse()
                    .map_err(|_| SelectionError::ParseError(format!("Invalid end: {}", range[1])))?;

                if start < 1 || end < 1 || start > end {
                    return Err(SelectionError::InvalidRange(format!("{}-{}", start, end)));
                }

                // Convert 1-based to 0-based
                for i in (start - 1)..end {
                    if i < n_atoms {
                        indices.insert(i);
                    }
                }
            } else {
                // Single atom
                let atom: usize = part
                    .parse()
                    .map_err(|_| SelectionError::ParseError(format!("Invalid atom: {}", part)))?;

                if atom < 1 {
                    return Err(SelectionError::InvalidRange("Atom numbers must be >= 1".to_string()));
                }

                // Convert 1-based to 0-based
                if atom - 1 < n_atoms {
                    indices.insert(atom - 1);
                }
            }
        }

        let mut indices_vec: Vec<usize> = indices.into_iter().collect();
        indices_vec.sort_unstable();

        Self::from_indices(indices_vec, n_atoms)
    }

    /// Parse prefixed selection: "1:1-10", "s:1-10", "r:1-5", "a:CA"
    fn parse_prefixed(spec: &str, topology: &Topology) -> Result<Self, SelectionError> {
        let parts: Vec<&str> = spec.splitn(2, ':').collect();
        if parts.len() != 2 {
            return Err(SelectionError::ParseError(format!("Invalid format: {}", spec)));
        }

        let prefix = parts[0].trim();
        let selector = parts[1].trim();

        match prefix {
            "a" => Self::parse_atom_name(selector, topology),
            "r" => Self::parse_residue(selector, topology),
            "s" => Self::parse_solvent(selector, topology),
            mol_str => {
                // Molecule number
                let mol: usize = mol_str
                    .parse()
                    .map_err(|_| SelectionError::ParseError(format!("Invalid molecule: {}", mol_str)))?;
                Self::parse_molecule(mol, selector, topology)
            }
        }
    }

    /// Parse atom name selection: "a:CA", "a:CA,CB,N"
    fn parse_atom_name(names: &str, topology: &Topology) -> Result<Self, SelectionError> {
        let name_list: Vec<&str> = names.split(',').map(|s| s.trim()).collect();
        let mut indices = Vec::new();

        for (idx, atom) in topology.solute.atoms.iter().enumerate() {
            if name_list.contains(&atom.name.as_str()) {
                indices.push(idx);
            }
        }

        if indices.is_empty() {
            return Err(SelectionError::InvalidAtomName(names.to_string()));
        }

        Ok(Self {
            indices,
            n_atoms: topology.num_atoms(),
        })
    }

    /// Parse residue selection: "r:1-5"
    fn parse_residue(spec: &str, topology: &Topology) -> Result<Self, SelectionError> {
        let residues = Self::parse_range_or_list(spec)?;
        let mut indices = Vec::new();

        for (idx, atom) in topology.solute.atoms.iter().enumerate() {
            if residues.contains(&atom.residue_nr) {
                indices.push(idx);
            }
        }

        if indices.is_empty() {
            return Err(SelectionError::InvalidResidue(residues[0]));
        }

        Ok(Self {
            indices,
            n_atoms: topology.num_atoms(),
        })
    }

    /// Parse solvent selection: "s:1-10"
    fn parse_solvent(_spec: &str, topology: &Topology) -> Result<Self, SelectionError> {
        // For now, select all solvent atoms
        // TODO: Implement proper solvent molecule selection
        let n_solute = topology.solute.atoms.len();
        let indices: Vec<usize> = (n_solute..topology.num_atoms()).collect();

        if indices.is_empty() {
            return Err(SelectionError::EmptySelection);
        }

        Ok(Self {
            indices,
            n_atoms: topology.num_atoms(),
        })
    }

    /// Parse molecule selection: "1:1-10", "1:a"
    fn parse_molecule(mol: usize, spec: &str, topology: &Topology) -> Result<Self, SelectionError> {
        // Validate molecule number
        if mol < 1 {
            return Err(SelectionError::InvalidMolecule(mol));
        }

        // For now, assume molecule 1 is the solute
        // TODO: Properly handle multiple molecules
        if mol != 1 {
            return Err(SelectionError::InvalidMolecule(mol));
        }

        if spec == "a" {
            // All atoms in molecule
            let indices: Vec<usize> = (0..topology.solute.atoms.len()).collect();
            Ok(Self {
                indices,
                n_atoms: topology.num_atoms(),
            })
        } else {
            // Atom range within molecule
            Self::parse_atom_list(spec, topology.solute.atoms.len())
        }
    }

    /// Parse range or list: "1-10" or "1,3,5-10"
    fn parse_range_or_list(spec: &str) -> Result<Vec<usize>, SelectionError> {
        let mut values = HashSet::new();

        for part in spec.split(',') {
            let part = part.trim();

            if part.contains('-') {
                let range: Vec<&str> = part.split('-').collect();
                if range.len() != 2 {
                    return Err(SelectionError::ParseError(format!("Invalid range: {}", part)));
                }

                let start: usize = range[0]
                    .trim()
                    .parse()
                    .map_err(|_| SelectionError::ParseError(format!("Invalid start: {}", range[0])))?;
                let end: usize = range[1]
                    .trim()
                    .parse()
                    .map_err(|_| SelectionError::ParseError(format!("Invalid end: {}", range[1])))?;

                if start > end {
                    return Err(SelectionError::InvalidRange(format!("{}-{}", start, end)));
                }

                for i in start..=end {
                    values.insert(i);
                }
            } else {
                let val: usize = part
                    .parse()
                    .map_err(|_| SelectionError::ParseError(format!("Invalid value: {}", part)))?;
                values.insert(val);
            }
        }

        let mut values_vec: Vec<usize> = values.into_iter().collect();
        values_vec.sort_unstable();
        Ok(values_vec)
    }

    /// Get selected atom indices (0-based)
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }

    /// Get number of selected atoms
    pub fn len(&self) -> usize {
        self.indices.len()
    }

    /// Check if selection is empty
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    /// Iterate over selected atom indices
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.indices.iter().copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::{Atom, Solute};

    fn create_test_topology() -> Topology {
        let mut topo = Topology::new();

        // Add 10 test atoms
        for i in 0..10 {
            topo.solute.atoms.push(Atom {
                name: if i % 3 == 0 { "CA".to_string() } else { "CB".to_string() },
                residue_nr: (i / 3) + 1,
                residue_name: "ALA".to_string(),
                iac: 1,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
            topo.iac.push(1);
            topo.mass.push(12.0);
            topo.charge.push(0.0);
        }

        topo
    }

    #[test]
    fn test_all_selection() {
        let topo = create_test_topology();
        let sel = AtomSelection::all(10);
        assert_eq!(sel.len(), 10);
        assert_eq!(sel.indices(), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    }

    #[test]
    fn test_atom_list() {
        let topo = create_test_topology();
        let sel = AtomSelection::from_string("1-5", &topo).unwrap();
        assert_eq!(sel.len(), 5);
        assert_eq!(sel.indices(), &[0, 1, 2, 3, 4]);
    }

    #[test]
    fn test_atom_list_with_commas() {
        let topo = create_test_topology();
        let sel = AtomSelection::from_string("1,3,5-7", &topo).unwrap();
        assert_eq!(sel.len(), 5);
        assert_eq!(sel.indices(), &[0, 2, 4, 5, 6]);
    }

    #[test]
    fn test_atom_name_selection() {
        let topo = create_test_topology();
        let sel = AtomSelection::from_string("a:CA", &topo).unwrap();
        assert_eq!(sel.len(), 4); // Atoms 0, 3, 6, 9
        assert_eq!(sel.indices(), &[0, 3, 6, 9]);
    }

    #[test]
    fn test_residue_selection() {
        let topo = create_test_topology();
        let sel = AtomSelection::from_string("r:1-2", &topo).unwrap();
        assert_eq!(sel.len(), 6); // Residues 1 and 2
    }

    #[test]
    fn test_molecule_all() {
        let topo = create_test_topology();
        let sel = AtomSelection::from_string("1:a", &topo).unwrap();
        assert_eq!(sel.len(), 10);
    }

    #[test]
    fn test_molecule_range() {
        let topo = create_test_topology();
        let sel = AtomSelection::from_string("1:1-5", &topo).unwrap();
        assert_eq!(sel.len(), 5);
        assert_eq!(sel.indices(), &[0, 1, 2, 3, 4]);
    }
}
