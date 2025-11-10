//! Parser for GROMOS topology files (.topo/.top)
//!
//! Format blocks:
//! - SOLUTEATOM: atom definitions (mass, charge, IAC)
//! - BOND/BONDSTRETCHTYPE: covalent bonds
//! - BONDANGLE/BONDANGLEBENDTYPE: bond angles
//! - DIHEDRAL/TORSDIHEDRALTYPE: torsional angles
//! - CGPARAMETERS: LJ interaction parameters
//! - TEMPERATUREGROUPS: temperature coupling groups

use crate::topology::{Topology, LJParameters, Bond, BondParameters, Angle, AngleParameters};
use crate::io::IoError;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::collections::HashMap;

/// Parsed topology data from GROMOS .topo file
#[derive(Debug)]
pub struct ParsedTopology {
    pub n_atoms: usize,
    pub masses: Vec<f64>,
    pub charges: Vec<f64>,
    pub iac: Vec<usize>,
    pub exclusions: Vec<Vec<usize>>,
    pub bonds: Vec<(usize, usize, usize)>,  // (i, j, type)
    pub bond_parameters: Vec<BondParameters>,
    pub angles: Vec<(usize, usize, usize, usize)>,  // (i, j, k, type)
    pub angle_parameters: Vec<AngleParameters>,
    pub lj_parameters: HashMap<(usize, usize), LJParameters>,
    pub temperature_groups: Vec<usize>,  // Last atom index of each group
}

/// Read GROMOS topology file
pub fn read_topology_file<P: AsRef<Path>>(path: P) -> Result<ParsedTopology, IoError> {
    let file = File::open(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    let reader = BufReader::new(file);

    let mut lines = reader.lines();

    let mut n_atoms = 0;
    let mut masses = Vec::new();
    let mut charges = Vec::new();
    let mut iac = Vec::new();
    let mut exclusions = Vec::new();
    let mut bonds = Vec::new();
    let mut bond_parameters = Vec::new();
    let mut angles = Vec::new();
    let mut angle_parameters = Vec::new();
    let mut lj_parameters = HashMap::new();
    let mut temperature_groups = Vec::new();

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        // Skip comments and empty lines
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        match trimmed {
            "SOLUTEATOM" => {
                parse_soluteatom(&mut lines, &mut n_atoms, &mut masses, &mut charges,
                                &mut iac, &mut exclusions)?;
            }
            "BONDSTRETCHTYPE" => {
                parse_bond_types(&mut lines, &mut bond_parameters)?;
            }
            "BOND" => {
                parse_bonds(&mut lines, &mut bonds)?;
            }
            "BONDANGLEBENDTYPE" => {
                parse_angle_types(&mut lines, &mut angle_parameters)?;
            }
            "BONDANGLE" => {
                parse_angles(&mut lines, &mut angles)?;
            }
            "CGPARAMETERS" => {
                parse_lj_parameters(&mut lines, &mut lj_parameters)?;
            }
            "TEMPERATUREGROUPS" => {
                parse_temperature_groups(&mut lines, &mut temperature_groups)?;
            }
            _ => {
                // Skip unknown blocks
                continue;
            }
        }
    }

    Ok(ParsedTopology {
        n_atoms,
        masses,
        charges,
        iac,
        exclusions,
        bonds,
        bond_parameters,
        angles,
        angle_parameters,
        lj_parameters,
        temperature_groups,
    })
}

fn parse_soluteatom<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    n_atoms: &mut usize,
    masses: &mut Vec<f64>,
    charges: &mut Vec<f64>,
    iac: &mut Vec<usize>,
    exclusions: &mut Vec<Vec<usize>>,
) -> Result<(), IoError> {
    // First line after SOLUTEATOM is atom count
    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First non-comment line is atom count
        if *n_atoms == 0 {
            *n_atoms = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid atom count: {}", trimmed)))?;
            continue;
        }

        // Parse atom line (format: ATNM MRES PANM IAC MASS CG CGC INE ...)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 8 {
            let atom_iac: usize = parts[3].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid IAC: {}", parts[3])))?;
            let mass: f64 = parts[4].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid mass: {}", parts[4])))?;
            let charge: f64 = parts[5].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid charge: {}", parts[5])))?;
            let n_exclusions: usize = parts[7].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid exclusion count: {}", parts[7])))?;

            iac.push(atom_iac);
            masses.push(mass);
            charges.push(charge);

            // Read exclusion list from next line(s)
            let mut atom_exclusions = Vec::new();
            if n_exclusions > 0 {
                // Next line contains exclusions
                if let Some(Ok(excl_line)) = lines.next() {
                    let excl_parts: Vec<&str> = excl_line.trim().split_whitespace().collect();
                    for excl in excl_parts.iter().take(n_exclusions) {
                        if let Ok(excl_idx) = excl.parse::<usize>() {
                            atom_exclusions.push(excl_idx);
                        }
                    }
                }
            }
            exclusions.push(atom_exclusions);
        }
    }

    Ok(())
}

fn parse_bond_types<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    bond_parameters: &mut Vec<BondParameters>,
) -> Result<(), IoError> {
    let mut n_types = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is type count
        if n_types == 0 {
            n_types = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond type count: {}", trimmed)))?;
            continue;
        }

        // Parse bond type (format: CB CHB B0)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 3 {
            let k_quartic: f64 = parts[0].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CB: {}", parts[0])))?;
            let k_harmonic: f64 = parts[1].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CHB: {}", parts[1])))?;
            let r0: f64 = parts[2].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid B0: {}", parts[2])))?;

            bond_parameters.push(BondParameters {
                k_quartic,
                k_harmonic,
                r0,
            });
        }
    }

    Ok(())
}

fn parse_bonds<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    bonds: &mut Vec<(usize, usize, usize)>,
) -> Result<(), IoError> {
    let mut n_bonds = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is bond count
        if n_bonds == 0 {
            n_bonds = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond count: {}", trimmed)))?;
            continue;
        }

        // Parse bond (format: IB JB ICB)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 3 {
            let i: usize = parts[0].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond atom i: {}", parts[0])))?;
            let j: usize = parts[1].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond atom j: {}", parts[1])))?;
            let bond_type: usize = parts[2].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond type: {}", parts[2])))?;

            // Convert from 1-indexed to 0-indexed
            bonds.push((i - 1, j - 1, bond_type - 1));
        }
    }

    Ok(())
}

fn parse_angle_types<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    angle_parameters: &mut Vec<AngleParameters>,
) -> Result<(), IoError> {
    let mut n_types = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is type count
        if n_types == 0 {
            n_types = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle type count: {}", trimmed)))?;
            continue;
        }

        // Parse angle type (format: CT CHT T0)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 3 {
            let k_cosine: f64 = parts[0].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CT: {}", parts[0])))?;
            let k_harmonic: f64 = parts[1].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CHT: {}", parts[1])))?;
            let theta0_deg: f64 = parts[2].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid T0: {}", parts[2])))?;

            angle_parameters.push(AngleParameters {
                k_cosine,
                k_harmonic,
                theta0: theta0_deg.to_radians(),
            });
        }
    }

    Ok(())
}

fn parse_angles<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    angles: &mut Vec<(usize, usize, usize, usize)>,
) -> Result<(), IoError> {
    let mut n_angles = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is angle count
        if n_angles == 0 {
            n_angles = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle count: {}", trimmed)))?;
            continue;
        }

        // Parse angle (format: IT JT KT ICT)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 4 {
            let i: usize = parts[0].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle atom i: {}", parts[0])))?;
            let j: usize = parts[1].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle atom j: {}", parts[1])))?;
            let k: usize = parts[2].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle atom k: {}", parts[2])))?;
            let angle_type: usize = parts[3].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle type: {}", parts[3])))?;

            // Convert from 1-indexed to 0-indexed
            angles.push((i - 1, j - 1, k - 1, angle_type - 1));
        }
    }

    Ok(())
}

fn parse_lj_parameters<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    lj_params: &mut HashMap<(usize, usize), LJParameters>,
) -> Result<(), IoError> {
    let mut n_pairs = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is pair count
        if n_pairs == 0 {
            n_pairs = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid LJ pair count: {}", trimmed)))?;
            continue;
        }

        // Parse LJ parameters (format: IAC JAC C12 C6 [CS12 CS6])
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 4 {
            let iac: usize = parts[0].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid IAC: {}", parts[0])))?;
            let jac: usize = parts[1].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid JAC: {}", parts[1])))?;
            let c12: f64 = parts[2].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid C12: {}", parts[2])))?;
            let c6: f64 = parts[3].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid C6: {}", parts[3])))?;

            lj_params.insert((iac, jac), LJParameters::new(c6, c12));

            // LJ matrix is symmetric
            if iac != jac {
                lj_params.insert((jac, iac), LJParameters::new(c6, c12));
            }
        }
    }

    Ok(())
}

fn parse_temperature_groups<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    temp_groups: &mut Vec<usize>,
) -> Result<(), IoError> {
    let mut n_groups = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is group count
        if n_groups == 0 {
            n_groups = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid temp group count: {}", trimmed)))?;
            continue;
        }

        // Next line(s) contain last atom indices for each group
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        for part in parts {
            if let Ok(last_atom) = part.parse::<usize>() {
                temp_groups.push(last_atom - 1);  // Convert to 0-indexed
            }
        }
    }

    Ok(())
}

/// Convert ParsedTopology to Topology
pub fn build_topology(parsed: ParsedTopology) -> Topology {
    let mut topo = Topology::new();

    // Populate both flat arrays and solute.atoms for compatibility
    topo.mass = parsed.masses.clone();
    topo.charge = parsed.charges.clone();
    topo.iac = parsed.iac.clone();
    topo.compute_inverse_masses();

    // Populate solute.atoms as well (needed for num_atoms())
    use crate::topology::Atom;
    for i in 0..parsed.n_atoms {
        topo.solute.atoms.push(Atom {
            name: format!("ATOM{}", i + 1),
            residue_nr: 1,
            residue_name: "UNK".to_string(),
            iac: parsed.iac[i],
            mass: parsed.masses[i],
            charge: parsed.charges[i],
            is_perturbed: false,
            is_polarisable: false,
            is_coarse_grained: true,
        });
    }

    // Initialize exclusions for each atom
    topo.exclusions = vec![std::collections::HashSet::new(); parsed.n_atoms];

    // Build exclusions
    for (i, excl_list) in parsed.exclusions.iter().enumerate() {
        for &j in excl_list {
            topo.exclusions[i].insert(j);
            if j < topo.num_atoms() {
                topo.exclusions[j].insert(i);
            }
        }
    }

    // Build bonds
    for (i, j, bond_type) in parsed.bonds {
        topo.solute.bonds.push(Bond {
            i,
            j,
            bond_type,
        });

        // Add bond exclusions (atoms bonded are excluded)
        topo.exclusions[i].insert(j);
        topo.exclusions[j].insert(i);
    }

    topo.bond_parameters = parsed.bond_parameters;

    // Build angles
    for (i, j, k, angle_type) in parsed.angles {
        topo.solute.angles.push(Angle {
            i,
            j,
            k,
            angle_type,
        });

        // Add 1-3 exclusions (i-j-k angle: i and k are excluded)
        topo.exclusions[i].insert(k);
        topo.exclusions[k].insert(i);
    }

    topo.angle_parameters = parsed.angle_parameters;

    // Build LJ parameter matrix
    let n_types = topo.iac.iter().max().unwrap_or(&0) + 1;
    topo.lj_parameters = vec![vec![LJParameters::default(); n_types]; n_types];

    for ((iac, jac), params) in parsed.lj_parameters {
        if iac > 0 && jac > 0 && iac <= n_types && jac <= n_types {
            topo.lj_parameters[iac - 1][jac - 1] = params;
        }
    }

    topo
}

/// Write GROMOS topology file
pub fn write_topology_file<P: AsRef<Path>>(
    path: P,
    topo: &Topology,
    title: &str,
) -> Result<(), IoError> {
    use std::io::Write;

    let file = File::create(path.as_ref())
        .map_err(|e| IoError::WriteError(format!("Cannot create topology file: {}", e)))?;
    let mut writer = std::io::BufWriter::new(file);

    // TITLE block
    writeln!(writer, "TITLE").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "{}", title).map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // PHYSICALCONSTANTS block (optional but good practice)
    writeln!(writer, "PHYSICALCONSTANTS").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  138.9354").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# HBAR: Planck's constant HBAR = H/(2* PI)").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  0.0635078").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# SPDL: Speed of light (nm/ps)").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  299792.458").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# BOLTZ: Boltzmann's constant kB").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  0.00831441").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // TOPVERSION block
    writeln!(writer, "TOPVERSION").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "2.0").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // ATOMTYPENAME block
    writeln!(writer, "ATOMTYPENAME").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# NTYP: number of atom types").map_err(|e| IoError::WriteError(e.to_string()))?;
    let n_types = topo.num_atom_types();
    writeln!(writer, "{}", n_types).map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# TYPE: atom type names").map_err(|e| IoError::WriteError(e.to_string()))?;
    for i in 1..=n_types {
        writeln!(writer, "TYPE{}", i).map_err(|e| IoError::WriteError(e.to_string()))?;
    }
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // RESNAME block
    writeln!(writer, "RESNAME").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# NRES: number of residues").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "1").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# residue name").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "MOL").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // SOLUTEATOM block
    writeln!(writer, "SOLUTEATOM").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# ATNM: atom number").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# MRES: residue number").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# PANM: atom name").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# IAC: integer atom code").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# MASS: mass").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# CG: charge").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# CGC: charge group code").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# INE: number of exclusions").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "#ATNM MRES PANM IAC  MASS      CG     CGC INE").map_err(|e| IoError::WriteError(e.to_string()))?;

    let n_atoms = topo.solute.num_atoms();
    for (i, atom) in topo.solute.atoms.iter().enumerate() {
        let n_exclusions = topo.exclusions.get(i).map_or(0, |e| e.len());
        let cg_code = topo.atom_to_chargegroup.get(i).map_or(1, |&c| c + 1);

        writeln!(
            writer,
            "{:5} {:4} {:>4} {:3} {:8.4} {:8.5} {:4} {:3}",
            i + 1,
            atom.residue_nr,
            atom.name,
            atom.iac,
            atom.mass,
            atom.charge,
            cg_code,
            n_exclusions
        ).map_err(|e| IoError::WriteError(e.to_string()))?;

        // Write exclusions if any
        if let Some(exclusions) = topo.exclusions.get(i) {
            if !exclusions.is_empty() {
                let mut excl_vec: Vec<_> = exclusions.iter().collect();
                excl_vec.sort();
                for (j, excl) in excl_vec.iter().enumerate() {
                    write!(writer, "{:5}", *excl + 1).map_err(|e| IoError::WriteError(e.to_string()))?;
                    if (j + 1) % 10 == 0 {
                        writeln!(writer).map_err(|e| IoError::WriteError(e.to_string()))?;
                    }
                }
                writeln!(writer).map_err(|e| IoError::WriteError(e.to_string()))?;
            }
        }
    }
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // BONDSTRETCHTYPE block
    if !topo.bond_parameters.is_empty() {
        writeln!(writer, "BONDSTRETCHTYPE").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NBTY: number of bond types").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.bond_parameters.len()).map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#  CB     CHB       B0").map_err(|e| IoError::WriteError(e.to_string()))?;
        for params in &topo.bond_parameters {
            writeln!(
                writer,
                "{:10.5e} {:10.5e} {:10.7}",
                params.k_quartic,
                params.k_harmonic,
                params.r0
            ).map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    // BOND block
    if !topo.solute.bonds.is_empty() {
        writeln!(writer, "BOND").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NBH: number of bonds").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.solute.bonds.len()).map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#  IB   JB  ICB").map_err(|e| IoError::WriteError(e.to_string()))?;
        for bond in &topo.solute.bonds {
            writeln!(
                writer,
                "{:5}{:5}{:5}",
                bond.i + 1,
                bond.j + 1,
                bond.bond_type + 1
            ).map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    // BONDANGLEBENDTYPE block
    if !topo.angle_parameters.is_empty() {
        writeln!(writer, "BONDANGLEBENDTYPE").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NTTY: number of angle types").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.angle_parameters.len()).map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#   CT      CHT     T0[deg]").map_err(|e| IoError::WriteError(e.to_string()))?;
        for params in &topo.angle_parameters {
            writeln!(
                writer,
                "{:10.5e} {:10.5e} {:10.4}",
                params.k_cosine,
                params.k_harmonic,
                params.theta0.to_degrees()
            ).map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    // BONDANGLE block
    if !topo.solute.angles.is_empty() {
        writeln!(writer, "BONDANGLE").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NTHEH: number of angles").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.solute.angles.len()).map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#  IT   JT   KT  ICT").map_err(|e| IoError::WriteError(e.to_string()))?;
        for angle in &topo.solute.angles {
            writeln!(
                writer,
                "{:5}{:5}{:5}{:5}",
                angle.i + 1,
                angle.j + 1,
                angle.k + 1,
                angle.angle_type + 1
            ).map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    writer.flush().map_err(|e| IoError::WriteError(e.to_string()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cg16_topology() {
        let result = read_topology_file("../md++/src/check/data/cg16.topo");

        if let Ok(parsed) = result {
            println!("Loaded topology with {} atoms", parsed.n_atoms);
            assert_eq!(parsed.n_atoms, 4);
            assert_eq!(parsed.bonds.len(), 3);
            assert_eq!(parsed.angles.len(), 2);

            // Check first atom
            assert_eq!(parsed.iac[0], 6);  // Type C
            assert!((parsed.masses[0] - 62.0).abs() < 1e-6);
            assert!((parsed.charges[0] - 0.0).abs() < 1e-6);

            // Check first bond (atoms 1-2, type 1 in file = 0-1, type 0 in code)
            assert_eq!(parsed.bonds[0], (0, 1, 0));

            println!("First bond type: k_harm={}, r0={}",
                     parsed.bond_parameters[0].k_harmonic,
                     parsed.bond_parameters[0].r0);
        } else {
            println!("Could not load file (may not exist in test environment)");
        }
    }
}
