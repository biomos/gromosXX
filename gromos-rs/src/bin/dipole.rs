//! dipole - Calculate Dipole Moment
//!
//! Usage: dipole @topo <topology> @traj <trajectory> @atoms <selection> [@cog]
//!
//! Calculates the dipole moment for selected atoms over a trajectory.
//! The dipole moment measures the polarity and charge distribution in a molecule.
//!
//! ## Formula
//!
//! μ = Σ(qi * ri)
//!
//! Where:
//! - μ = dipole moment vector (e·nm)
//! - qi = partial charge of atom i (e)
//! - ri = position of atom i (nm)
//!
//! |μ| = magnitude of dipole (reported in Debye)
//! 1 Debye = 0.20819434 e·Å = 0.020819434 e·nm
//!
//! ## Important Notes
//! - For molecules with net charge, dipole depends on the origin
//! - Use @cog to center at center-of-geometry for consistent results

use gromos_rs::io::topology::{read_topology_file, build_topology};
use gromos_rs::io::trajectory::TrajectoryReader;
use gromos_rs::math::Vec3;
use std::env;
use std::process;

// Conversion factor: e·nm to Debye
const ENM_TO_DEBYE: f64 = 48.0321;

fn print_usage() {
    eprintln!("dipole - Dipole Moment calculation");
    eprintln!();
    eprintln!("Usage: dipole @topo <topology> @traj <trajectory> @atoms <selection> [@cog]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo    Molecular topology file");
    eprintln!("  @traj    Input trajectory file (.trc)");
    eprintln!("  @atoms   Atom selection (default: all)");
    eprintln!("           Format: 'all', '1-10', '1,5,10-20'");
    eprintln!("  @cog     Center at center-of-geometry (recommended for charged molecules)");
    eprintln!();
    eprintln!("Description:");
    eprintln!("  Calculates dipole moment: μ = Σ(qi * ri)");
    eprintln!("  - μ = dipole moment vector (e·nm)");
    eprintln!("  - qi = atomic charge (e)");
    eprintln!("  - ri = atomic position (nm)");
    eprintln!();
    eprintln!("  Output in Debye (D): 1 D = 0.020819434 e·nm");
    eprintln!();
    eprintln!("  ⚠️  For charged molecules: Use @cog to avoid origin-dependent results");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  Time (ps)  |μ| (D)  μx (D)  μy (D)  μz (D)  Net_charge (e)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  # Calculate dipole for all atoms");
    eprintln!("  dipole @topo system.top @traj output.trc @atoms all");
    eprintln!();
    eprintln!("  # Protein molecule with center-of-geometry centering");
    eprintln!("  dipole @topo system.top @traj output.trc @atoms 1-1000 @cog");
}

#[derive(Debug)]
struct DipoleArgs {
    topo_file: String,
    traj_file: String,
    atom_selection: Vec<usize>,
    use_cog: bool,
}

fn parse_atom_selection(spec: &str, n_atoms: usize) -> Result<Vec<usize>, String> {
    if spec == "all" {
        return Ok((0..n_atoms).collect());
    }

    let mut atoms = Vec::new();

    for part in spec.split(',') {
        if part.contains('-') {
            let range: Vec<&str> = part.split('-').collect();
            if range.len() != 2 {
                return Err(format!("Invalid range: {}", part));
            }
            let start: usize = range[0]
                .parse()
                .map_err(|_| format!("Invalid start: {}", range[0]))?;
            let end: usize = range[1]
                .parse()
                .map_err(|_| format!("Invalid end: {}", range[1]))?;

            if start < 1 || end < 1 || start > end {
                return Err(format!("Invalid range: {}-{}", start, end));
            }

            // Convert from 1-based to 0-based
            for i in (start - 1)..end {
                if i < n_atoms {
                    atoms.push(i);
                }
            }
        } else {
            let atom: usize = part
                .parse()
                .map_err(|_| format!("Invalid atom: {}", part))?;
            if atom < 1 {
                return Err("Atom numbers must be >= 1".to_string());
            }
            // Convert from 1-based to 0-based
            if atom - 1 < n_atoms {
                atoms.push(atom - 1);
            }
        }
    }

    if atoms.is_empty() {
        return Err("No valid atoms in selection".to_string());
    }

    atoms.sort_unstable();
    atoms.dedup();
    Ok(atoms)
}

fn parse_args(args: Vec<String>) -> Result<DipoleArgs, String> {
    let mut topo_file = None;
    let mut traj_file = None;
    let mut atom_spec = None;
    let mut use_cog = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing value for @topo".to_string());
                }
                topo_file = Some(args[i].clone());
            }
            "@traj" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing value for @traj".to_string());
                }
                traj_file = Some(args[i].clone());
            }
            "@atoms" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing value for @atoms".to_string());
                }
                atom_spec = Some(args[i].clone());
            }
            "@cog" => {
                use_cog = true;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    let topo_file = topo_file.ok_or("Missing required argument: @topo")?;
    let traj_file = traj_file.ok_or("Missing required argument: @traj")?;

    Ok(DipoleArgs {
        topo_file,
        traj_file,
        atom_selection: vec![],
        use_cog,
    })
}

/// Calculate center of geometry (unweighted center)
fn calc_center_of_geometry(positions: &[Vec3], atom_selection: &[usize]) -> Vec3 {
    let mut cog = Vec3::ZERO;
    for &idx in atom_selection {
        cog = cog + positions[idx];
    }
    cog / (atom_selection.len() as f32)
}

/// Calculate dipole moment
fn calc_dipole(
    positions: &[Vec3],
    charges: &[f64],
    atom_selection: &[usize],
    use_cog: bool,
) -> (Vec3, f64) {
    // Calculate center of geometry if requested
    let origin = if use_cog {
        calc_center_of_geometry(positions, atom_selection)
    } else {
        Vec3::ZERO
    };

    // Calculate dipole moment: μ = Σ(qi * ri)
    let mut dipole = Vec3::ZERO;
    let mut net_charge = 0.0;

    for &idx in atom_selection {
        let r = positions[idx] - origin;
        let q = charges[idx] as f32;
        dipole = dipole + r * q;
        net_charge += charges[idx];
    }

    // Convert to f64 for magnitude calculation
    let dipole_x = dipole.x as f64;
    let dipole_y = dipole.y as f64;
    let dipole_z = dipole.z as f64;

    // Calculate magnitude in e·nm
    let magnitude_enm = (dipole_x * dipole_x + dipole_y * dipole_y + dipole_z * dipole_z).sqrt();

    // Convert to Debye
    let magnitude_debye = magnitude_enm * ENM_TO_DEBYE;

    (dipole, magnitude_debye)
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args[1] == "-h" || args[1] == "--help" {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let mut dipole_args = match parse_args(args) {
        Ok(args) => args,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    // Read topology
    eprintln!("Reading topology: {}", dipole_args.topo_file);
    let blocks = match read_topology_file(&dipole_args.topo_file) {
        Ok(blocks) => blocks,
        Err(e) => {
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        }
    };

    let topo = build_topology(blocks);
    let n_atoms = topo.num_atoms();
    eprintln!("  Total atoms: {}", n_atoms);

    // Parse atom selection
    dipole_args.atom_selection = match parse_atom_selection(&"all", n_atoms) {
        Ok(sel) => sel,
        Err(e) => {
            eprintln!("Error parsing atom selection: {}", e);
            process::exit(1);
        }
    };

    eprintln!("  Selected atoms: {}", dipole_args.atom_selection.len());

    // Calculate net charge
    let mut net_charge = 0.0;
    for &idx in &dipole_args.atom_selection {
        net_charge += topo.charge[idx];
    }
    eprintln!("  Net charge: {:.6} e", net_charge);

    if net_charge.abs() > 0.01 && !dipole_args.use_cog {
        eprintln!("  ⚠️  WARNING: Molecule has net charge. Consider using @cog flag.");
    }

    eprintln!("  Center: {}", if dipole_args.use_cog { "center-of-geometry" } else { "origin" });
    eprintln!();

    // Open trajectory
    eprintln!("Reading trajectory: {}", dipole_args.traj_file);
    let mut traj = match TrajectoryReader::new(&dipole_args.traj_file) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("Error opening trajectory: {}", e);
            process::exit(1);
        }
    };

    // Output header
    println!("# Time (ps)    |μ| (D)      μx (D)      μy (D)      μz (D)  Net_Q (e)");

    let mut frame_count = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                frame_count += 1;

                // Calculate dipole
                let (dipole_vec, magnitude) = calc_dipole(
                    &frame.positions,
                    &topo.charge,
                    &dipole_args.atom_selection,
                    dipole_args.use_cog,
                );

                // Convert components to Debye
                let mu_x = dipole_vec.x as f64 * ENM_TO_DEBYE;
                let mu_y = dipole_vec.y as f64 * ENM_TO_DEBYE;
                let mu_z = dipole_vec.z as f64 * ENM_TO_DEBYE;

                println!(
                    "{:12.6}{:12.6}{:12.6}{:12.6}{:12.6}{:10.6}",
                    frame.time, magnitude, mu_x, mu_y, mu_z, net_charge
                );
            }
            Ok(None) => break,
            Err(e) => {
                eprintln!("Error reading frame {}: {}", frame_count, e);
                process::exit(1);
            }
        }
    }

    eprintln!();
    eprintln!("Processed {} frames", frame_count);
}
