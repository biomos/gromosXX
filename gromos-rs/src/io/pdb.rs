//! PDB file format reader/writer
//!
//! Implements reading and writing of Protein Data Bank (PDB) format files.
//! Based on PDB format specification v3.3

use crate::math::Vec3;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// A single atom record from a PDB file
#[derive(Debug, Clone)]
pub struct PDBAtom {
    /// Serial number
    pub serial: usize,
    /// Atom name
    pub name: String,
    /// Residue name
    pub resname: String,
    /// Chain identifier
    pub chain: String,
    /// Residue sequence number
    pub resnum: i32,
    /// Coordinates (Angstroms in PDB, converted to nm)
    pub coord: Vec3,
    /// Occupancy
    pub occupancy: f64,
    /// Temperature factor (B-factor)
    pub bfactor: f64,
    /// Element symbol
    pub element: String,
}

impl PDBAtom {
    /// Parse a PDB ATOM or HETATM record line
    ///
    /// PDB format uses fixed-width columns:
    /// - Columns 1-6: Record type (ATOM/HETATM)
    /// - Columns 7-11: Atom serial number
    /// - Columns 13-16: Atom name
    /// - Columns 18-20: Residue name
    /// - Column 22: Chain ID
    /// - Columns 23-26: Residue sequence number
    /// - Columns 31-38: X coordinate (Angstrom)
    /// - Columns 39-46: Y coordinate (Angstrom)
    /// - Columns 47-54: Z coordinate (Angstrom)
    /// - Columns 55-60: Occupancy
    /// - Columns 61-66: B-factor
    /// - Columns 77-78: Element symbol
    pub fn from_pdb_line(line: &str) -> Result<Self, String> {
        if line.len() < 54 {
            return Err(format!("PDB line too short: {}", line.len()));
        }

        let record_type = line[0..6].trim();
        if record_type != "ATOM" && record_type != "HETATM" {
            return Err(format!("Not an ATOM/HETATM record: {}", record_type));
        }

        // Parse serial number (columns 7-11)
        let serial = line[6..11].trim().parse::<usize>()
            .map_err(|e| format!("Invalid serial number: {}", e))?;

        // Parse atom name (columns 13-16)
        let name = line[12..16].trim().to_string();

        // Parse residue name (columns 18-20)
        let resname = line[17..20].trim().to_string();

        // Parse chain ID (column 22)
        let chain = if line.len() > 21 {
            line[21..22].trim().to_string()
        } else {
            String::new()
        };

        // Parse residue number (columns 23-26)
        let resnum = line[22..26].trim().parse::<i32>()
            .map_err(|e| format!("Invalid residue number: {}", e))?;

        // Parse coordinates (columns 31-38, 39-46, 47-54)
        let x = line[30..38].trim().parse::<f64>()
            .map_err(|e| format!("Invalid X coordinate: {}", e))?;
        let y = line[38..46].trim().parse::<f64>()
            .map_err(|e| format!("Invalid Y coordinate: {}", e))?;
        let z = line[46..54].trim().parse::<f64>()
            .map_err(|e| format!("Invalid Z coordinate: {}", e))?;

        // Convert from Angstrom to nanometers
        let coord = Vec3::new((x * 0.1) as f32, (y * 0.1) as f32, (z * 0.1) as f32);

        // Parse occupancy (columns 55-60)
        let occupancy = if line.len() >= 60 {
            line[54..60].trim().parse::<f64>().unwrap_or(1.0)
        } else {
            1.0
        };

        // Parse B-factor (columns 61-66)
        let bfactor = if line.len() >= 66 {
            line[60..66].trim().parse::<f64>().unwrap_or(0.0) * 0.01 // Convert to nm²
        } else {
            0.0
        };

        // Parse element (columns 77-78)
        let element = if line.len() >= 78 {
            line[76..78].trim().to_string()
        } else {
            // Guess element from atom name
            name.chars().next()
                .filter(|c| c.is_alphabetic())
                .map(|c| c.to_string())
                .unwrap_or_default()
        };

        Ok(PDBAtom {
            serial,
            name,
            resname,
            chain,
            resnum,
            coord,
            occupancy,
            bfactor,
            element,
        })
    }

    /// Write atom as a PDB ATOM record
    pub fn to_pdb_line(&self, use_hetatm: bool) -> String {
        let record = if use_hetatm { "HETATM" } else { "ATOM  " };

        // Convert coordinates back to Angstroms
        let x = (self.coord.x as f64) * 10.0;
        let y = (self.coord.y as f64) * 10.0;
        let z = (self.coord.z as f64) * 10.0;

        format!(
            "{:6}{:5} {:^4} {:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}",
            record,
            self.serial,
            self.name,
            self.resname,
            self.chain,
            self.resnum,
            x, y, z,
            self.occupancy,
            self.bfactor * 100.0, // Convert back from nm²
            self.element
        )
    }
}

/// A residue containing multiple atoms
#[derive(Debug, Clone)]
pub struct PDBResidue {
    pub name: String,
    pub number: i32,
    pub chain: String,
    pub atoms: Vec<PDBAtom>,
}

/// A PDB structure
#[derive(Debug, Clone, Default)]
pub struct PDBStructure {
    pub title: Vec<String>,
    pub residues: Vec<PDBResidue>,
}

impl PDBStructure {
    /// Read a PDB file
    pub fn read_pdb<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let file = File::open(path)
            .map_err(|e| format!("Cannot open PDB file: {}", e))?;
        let reader = BufReader::new(file);

        let mut structure = PDBStructure::default();
        let mut current_residue: Option<PDBResidue> = None;

        for line in reader.lines() {
            let line = line.map_err(|e| format!("Error reading line: {}", e))?;

            if line.starts_with("TITLE") || line.starts_with("HEADER") {
                structure.title.push(line[6..].trim().to_string());
            } else if line.starts_with("ATOM") || line.starts_with("HETATM") {
                match PDBAtom::from_pdb_line(&line) {
                    Ok(atom) => {
                        // Check if we need to start a new residue
                        let start_new_residue = match &current_residue {
                            None => true,
                            Some(res) => res.number != atom.resnum || res.chain != atom.chain,
                        };

                        if start_new_residue {
                            // Save previous residue
                            if let Some(res) = current_residue.take() {
                                structure.residues.push(res);
                            }

                            // Start new residue
                            current_residue = Some(PDBResidue {
                                name: atom.resname.clone(),
                                number: atom.resnum,
                                chain: atom.chain.clone(),
                                atoms: vec![atom],
                            });
                        } else if let Some(ref mut res) = current_residue {
                            res.atoms.push(atom);
                        }
                    }
                    Err(e) => {
                        eprintln!("Warning: Skipping invalid PDB line: {}", e);
                    }
                }
            }
        }

        // Save last residue
        if let Some(res) = current_residue {
            structure.residues.push(res);
        }

        Ok(structure)
    }

    /// Write structure as PDB file
    pub fn write_pdb<P: AsRef<Path>>(&self, path: P) -> Result<(), String> {
        let mut file = File::create(path)
            .map_err(|e| format!("Cannot create PDB file: {}", e))?;

        // Write title
        for (i, title) in self.title.iter().enumerate() {
            writeln!(file, "TITLE   {:2} {}", i + 1, title)
                .map_err(|e| format!("Write error: {}", e))?;
        }

        // Write atoms
        for residue in &self.residues {
            for atom in &residue.atoms {
                writeln!(file, "{}", atom.to_pdb_line(false))
                    .map_err(|e| format!("Write error: {}", e))?;
            }
        }

        writeln!(file, "END")
            .map_err(|e| format!("Write error: {}", e))?;

        Ok(())
    }

    /// Get total number of atoms
    pub fn num_atoms(&self) -> usize {
        self.residues.iter().map(|r| r.atoms.len()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_pdb_atom() {
        let line = "ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  0.00           N";
        let atom = PDBAtom::from_pdb_line(line).unwrap();

        assert_eq!(atom.serial, 1);
        assert_eq!(atom.name, "N");
        assert_eq!(atom.resname, "MET");
        assert_eq!(atom.chain, "A");
        assert_eq!(atom.resnum, 1);
        assert!((atom.coord.x - 2.734).abs() < 1e-3); // 27.340 Å = 2.734 nm
        assert_eq!(atom.element, "N");
    }

    #[test]
    fn test_angstrom_to_nm_conversion() {
        let line = "ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00  0.00           C";
        let atom = PDBAtom::from_pdb_line(line).unwrap();

        assert!((atom.coord.x - 1.0).abs() < 1e-6); // 10 Å = 1 nm
        assert!((atom.coord.y - 2.0).abs() < 1e-6); // 20 Å = 2 nm
        assert!((atom.coord.z - 3.0).abs() < 1e-6); // 30 Å = 3 nm
    }
}
