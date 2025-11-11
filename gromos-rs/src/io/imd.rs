//! Parser for GROMOS input parameter files (.imd)
//!
//! The IMD format uses block-based keywords to specify simulation parameters.
//!
//! # Format Example
//! ```text
//! TITLE
//!   My simulation
//! END
//! SYSTEM
//!   NPM     1
//!   NSM     0
//! END
//! STEP
//!   NSTLIM  1000
//!   T       0.0
//!   DT      0.002
//! END
//! ```
//!
//! # Major Blocks
//! - TITLE: Simulation description
//! - SYSTEM: System composition (NPM, NSM)
//! - STEP: Integration parameters (NSTLIM, DT, T)
//! - BOUNDCOND: Boundary conditions (NTB, NDFMIN)
//! - MULTIBATH: Temperature coupling groups
//! - PRESSURESCALE: Pressure coupling
//! - FORCE: Force field terms to compute
//! - CONSTRAINT: Constraint algorithm settings
//! - PAIRLIST: Neighbor list parameters
//! - NONBONDED: Cutoffs and long-range electrostatics
//! - INITIALISE: Initial velocities
//! - WRITETRAJ: Trajectory output settings
//! - PRINTOUT: Energy output settings

use crate::io::IoError;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// GROMOS simulation parameters from .imd file
#[derive(Debug, Clone)]
pub struct ImdParameters {
    /// Simulation title
    pub title: String,

    // SYSTEM block
    pub npm: usize,  // Number of (identical) protein molecules
    pub nsm: usize,  // Number of (identical) solvent molecules

    // STEP block
    pub nstlim: usize,  // Number of MD steps
    pub t0: f64,         // Initial time (ps)
    pub dt: f64,         // Time step (ps)

    // BOUNDCOND block
    pub ntb: i32,    // Boundary type (0=vacuum, 1=rectangular, 2=triclinic, -1=truncated octahedron)
    pub ndfmin: i32, // Minimum number of degrees of freedom

    // MULTIBATH block
    pub num_temp_baths: usize,
    pub temp_bath: Vec<TempBathParameters>,

    // PRESSURESCALE block
    pub couple_pressure: bool,
    pub pressure_parameters: Option<PressureParameters>,

    // FORCE block
    pub force_groups: Vec<Vec<(usize, usize)>>,  // Energy group pairs

    // CONSTRAINT block
    pub ntc: i32,  // SHAKE constraints (1=none, 2=H-bonds, 3=all bonds, 4=all)
    pub ntcp: i32, // P-SHAKE (pressure-SHAKE)
    pub ntcs: i32, // Solvent SHAKE/SETTLE

    // PAIRLIST block
    pub algorithm: i32,  // Pairlist algorithm
    pub nsnb: usize,     // Update frequency
    pub rcutp: f64,      // Short-range cutoff (nm)
    pub rcutl: f64,      // Long-range cutoff (nm)
    pub size: f64,       // Grid cell size (nm)
    pub type_: i32,      // Pairlist type

    // NONBONDED block
    pub nlrele: i32,     // Long-range electrostatics (0=cutoff, 1=RF, 2=PME, 3=P3M)
    pub appak: f64,      // Reaction field κ (nm⁻¹)
    pub rcrf: f64,       // Reaction field cutoff (nm)
    pub epsrf: f64,      // Reaction field permittivity
    pub nslfexcl: i32,   // Exclusions

    // PME-specific parameters
    pub grid_x: usize,   // PME grid size X
    pub grid_y: usize,   // PME grid size Y
    pub grid_z: usize,   // PME grid size Z
    pub pme_order: usize,// PME spline order (typically 4)
    pub pme_alpha: f64,  // PME Ewald parameter

    // INITIALISE block
    pub ntivel: i32,     // Initial velocities (0=read, 1=generate)
    pub ntishk: i32,     // SHAKE initial configuration
    pub ntinht: i32,     // Initial temperature assignment
    pub ntinhb: i32,     // Initial bond constraints
    pub ntishi: i32,     // Initial SHAKE iterations
    pub ig: i64,         // Random seed
    pub tempi: f64,      // Initial temperature (K)

    // WRITETRAJ block
    pub ntwx: usize,     // Trajectory write frequency
    pub ntwe: usize,     // Energy write frequency
    pub ntwv: bool,      // Write velocities
    pub ntwf: bool,      // Write forces
    pub ntwe_special: bool, // Special energy format

    // PRINTOUT block
    pub ntpr: usize,     // Print frequency

    /// Raw blocks for custom parsing
    pub raw_blocks: HashMap<String, Vec<String>>,
}

/// Temperature bath parameters (MULTIBATH block)
#[derive(Debug, Clone)]
pub struct TempBathParameters {
    pub algorithm: i32,      // Coupling algorithm
    pub num_bath_groups: usize,
    pub temp0: Vec<f64>,     // Reference temperatures (K)
    pub tau: Vec<f64>,       // Coupling times (ps)
    pub dof: Vec<usize>,     // Degrees of freedom per group
}

/// Pressure bath parameters (PRESSURESCALE block)
#[derive(Debug, Clone)]
pub struct PressureParameters {
    pub algorithm: i32,      // Coupling algorithm (1=Berendsen, 2=Parrinello-Rahman)
    pub pressure0: [[f64; 3]; 3], // Reference pressure tensor (GPa/bar)
    pub compressibility: [[f64; 3]; 3], // Isothermal compressibility
    pub tau_p: f64,          // Coupling time (ps)
    pub virial: i32,         // Virial calculation method
}

impl Default for ImdParameters {
    fn default() -> Self {
        Self {
            title: String::from("GROMOS simulation"),
            npm: 1,
            nsm: 0,
            nstlim: 1000,
            t0: 0.0,
            dt: 0.002,
            ntb: 1,
            ndfmin: 0,
            num_temp_baths: 1,
            temp_bath: vec![TempBathParameters::default()],
            couple_pressure: false,
            pressure_parameters: None,
            force_groups: Vec::new(),
            ntc: 2,
            ntcp: 0,
            ntcs: 1,
            algorithm: 0,
            nsnb: 5,
            rcutp: 0.8,
            rcutl: 1.4,
            size: 0.4,
            type_: 0,
            nlrele: 1,
            appak: 0.0,
            rcrf: 1.4,
            epsrf: 0.0,
            nslfexcl: 0,
            grid_x: 64,
            grid_y: 64,
            grid_z: 64,
            pme_order: 4,
            pme_alpha: 0.0,
            ntivel: 0,
            ntishk: 0,
            ntinht: 0,
            ntinhb: 0,
            ntishi: 1000,
            ig: 12345,
            tempi: 300.0,
            ntwx: 100,
            ntwe: 100,
            ntwv: false,
            ntwf: false,
            ntwe_special: false,
            ntpr: 100,
            raw_blocks: HashMap::new(),
        }
    }
}

impl Default for TempBathParameters {
    fn default() -> Self {
        Self {
            algorithm: 1,
            num_bath_groups: 1,
            temp0: vec![300.0],
            tau: vec![0.1],
            dof: vec![0],
        }
    }
}

/// Parse GROMOS .imd parameter file
pub fn read_imd_file<P: AsRef<Path>>(path: P) -> Result<ImdParameters, IoError> {
    let file = File::open(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    let reader = BufReader::new(file);

    let mut params = ImdParameters::default();
    let mut current_block = String::new();
    let mut block_lines: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Skip comments and empty lines
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        // Check for block start/end
        if trimmed == "END" {
            // Process completed block
            if !current_block.is_empty() {
                parse_block(&mut params, &current_block, &block_lines)?;
                params.raw_blocks.insert(current_block.clone(), block_lines.clone());
                block_lines.clear();
            }
            current_block.clear();
        } else if !current_block.is_empty() {
            // Inside a block - collect lines
            block_lines.push(trimmed.to_string());
        } else {
            // New block starting
            current_block = trimmed.to_string();
        }
    }

    Ok(params)
}

/// Parse a specific IMD block
fn parse_block(params: &mut ImdParameters, block_name: &str, lines: &[String]) -> Result<(), IoError> {
    match block_name {
        "TITLE" => {
            params.title = lines.join(" ");
        },
        "SYSTEM" => {
            for line in lines {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    match parts[0] {
                        "NPM" => params.npm = parts[1].parse().unwrap_or(1),
                        "NSM" => params.nsm = parts[1].parse().unwrap_or(0),
                        _ => {}
                    }
                }
            }
        },
        "STEP" => {
            for line in lines {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    match parts[0] {
                        "NSTLIM" => params.nstlim = parts[1].parse().unwrap_or(1000),
                        "T" => params.t0 = parts[1].parse().unwrap_or(0.0),
                        "DT" => params.dt = parts[1].parse().unwrap_or(0.002),
                        _ => {}
                    }
                }
            }
        },
        "BOUNDCOND" => {
            for line in lines {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    match parts[0] {
                        "NTB" => params.ntb = parts[1].parse().unwrap_or(1),
                        "NDFMIN" => params.ndfmin = parts[1].parse().unwrap_or(0),
                        _ => {}
                    }
                }
            }
        },
        "MULTIBATH" => {
            // Simplified MULTIBATH parsing
            let mut bath = TempBathParameters::default();
            for line in lines {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    match parts[0] {
                        "ALGORITHM" => bath.algorithm = parts[1].parse().unwrap_or(1),
                        "NUM" | "NBATHS" => bath.num_bath_groups = parts[1].parse().unwrap_or(1),
                        "TEMP0" | "TEMP" => {
                            bath.temp0 = parts[1..].iter()
                                .filter_map(|s| s.parse().ok())
                                .collect();
                        },
                        "TAU" | "TAUT" => {
                            bath.tau = parts[1..].iter()
                                .filter_map(|s| s.parse().ok())
                                .collect();
                        },
                        _ => {}
                    }
                }
            }
            params.temp_bath = vec![bath];
        },
        "CONSTRAINT" => {
            for line in lines {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    match parts[0] {
                        "NTC" => params.ntc = parts[1].parse().unwrap_or(2),
                        "NTCP" => params.ntcp = parts[1].parse().unwrap_or(0),
                        "NTCS" => params.ntcs = parts[1].parse().unwrap_or(1),
                        _ => {}
                    }
                }
            }
        },
        "PAIRLIST" => {
            for line in lines {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    match parts[0] {
                        "ALGORITHM" => params.algorithm = parts[1].parse().unwrap_or(0),
                        "NSNB" => params.nsnb = parts[1].parse().unwrap_or(5),
                        "RCUTP" => params.rcutp = parts[1].parse().unwrap_or(0.8),
                        "RCUTL" => params.rcutl = parts[1].parse().unwrap_or(1.4),
                        "SIZE" => params.size = parts[1].parse().unwrap_or(0.4),
                        "TYPE" => params.type_ = parts[1].parse().unwrap_or(0),
                        _ => {}
                    }
                }
            }
        },
        "NONBONDED" => {
            for line in lines {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    match parts[0] {
                        "NLRELE" => params.nlrele = parts[1].parse().unwrap_or(1),
                        "APPAK" => params.appak = parts[1].parse().unwrap_or(0.0),
                        "RCRF" => params.rcrf = parts[1].parse().unwrap_or(1.4),
                        "EPSRF" => params.epsrf = parts[1].parse().unwrap_or(0.0),
                        "NSLFEXCL" => params.nslfexcl = parts[1].parse().unwrap_or(0),
                        // PME parameters
                        "GRIDX" => params.grid_x = parts[1].parse().unwrap_or(64),
                        "GRIDY" => params.grid_y = parts[1].parse().unwrap_or(64),
                        "GRIDZ" => params.grid_z = parts[1].parse().unwrap_or(64),
                        "ORDER" | "PMEORDER" => params.pme_order = parts[1].parse().unwrap_or(4),
                        "ALPHA" | "PMEALPHA" => params.pme_alpha = parts[1].parse().unwrap_or(0.0),
                        _ => {}
                    }
                }
            }
        },
        "INITIALISE" => {
            for line in lines {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    match parts[0] {
                        "NTIVEL" => params.ntivel = parts[1].parse().unwrap_or(0),
                        "NTISHK" => params.ntishk = parts[1].parse().unwrap_or(0),
                        "NTINHT" => params.ntinht = parts[1].parse().unwrap_or(0),
                        "NTINHB" => params.ntinhb = parts[1].parse().unwrap_or(0),
                        "NTISHI" => params.ntishi = parts[1].parse().unwrap_or(1000),
                        "IG" => params.ig = parts[1].parse().unwrap_or(12345),
                        "TEMPI" => params.tempi = parts[1].parse().unwrap_or(300.0),
                        _ => {}
                    }
                }
            }
        },
        "WRITETRAJ" => {
            for line in lines {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    match parts[0] {
                        "NTWX" => params.ntwx = parts[1].parse().unwrap_or(100),
                        "NTWE" => params.ntwe = parts[1].parse().unwrap_or(100),
                        "NTWV" => params.ntwv = parts[1].parse::<i32>().unwrap_or(0) != 0,
                        "NTWF" => params.ntwf = parts[1].parse::<i32>().unwrap_or(0) != 0,
                        _ => {}
                    }
                }
            }
        },
        "PRINTOUT" => {
            for line in lines {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    match parts[0] {
                        "NTPR" => params.ntpr = parts[1].parse().unwrap_or(100),
                        _ => {}
                    }
                }
            }
        },
        _ => {
            // Unknown block - store in raw_blocks for later use
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_parameters() {
        let params = ImdParameters::default();
        assert_eq!(params.nstlim, 1000);
        assert_eq!(params.dt, 0.002);
        assert_eq!(params.ntc, 2);  // SHAKE on H-bonds
    }
}
