//! rgyr - Calculate Radius of Gyration
//!
//! Usage: rgyr @topo <topology> @traj <trajectory> @atoms <selection> [@massweighted]
//!
//! Calculates the radius of gyration (Rg) for selected atoms over a trajectory.
//! Rg measures the compactness of a molecular structure and is commonly used
//! to track protein folding/unfolding and conformational changes.
//!
//! ## Formulas
//!
//! Unweighted: Rg = sqrt((1/N) * Σ(ri - rcom)²)
//! Mass-weighted: Rg = sqrt((1/M) * Σ(mi * (ri - rcom)²))
//!
//! Where:
//! - N = number of atoms
//! - ri = position of atom i
//! - rcom = center of mass
//! - M = total mass
//! - mi = mass of atom i

use gromos_rs::io::topology::{read_topology_file, build_topology};
use gromos_rs::io::trajectory::TrajectoryReader;
use gromos_rs::selection::AtomSelection;
use gromos_rs::math::Vec3;
use std::env;
use std::process;

fn print_usage() {
    eprintln!("rgyr - Radius of Gyration calculation");
    eprintln!();
    eprintln!("Usage: rgyr @topo <topology> @traj <trajectory> @atoms <selection> [@massweighted]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo          Molecular topology file");
    eprintln!("  @traj          Input trajectory file (.trc)");
    eprintln!("  @atoms         Atom selection (default: all)");
    eprintln!("                 Formats: 'all', '1-10', '1,5,10-20'");
    eprintln!("                          '1:1-10' (mol 1, atoms 1-10)");
    eprintln!("                          'r:1-5' (residues 1-5)");
    eprintln!("                          'a:CA' (all CA atoms)");
    eprintln!("  @massweighted  Use mass-weighted formula (default: unweighted)");
    eprintln!();
    eprintln!("Description:");
    eprintln!("  Calculates radius of gyration over trajectory:");
    eprintln!("  - Unweighted: Rg = sqrt((1/N) * Σ(ri - rcom)²)");
    eprintln!("  - Mass-weighted: Rg = sqrt((1/M) * Σ(mi * (ri - rcom)²))");
    eprintln!();
    eprintln!("  Rg measures molecular compactness:");
    eprintln!("  - Smaller Rg → more compact structure");
    eprintln!("  - Larger Rg → more extended structure");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  Time (ps)  Rg (nm)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  # Calculate Rg for all atoms");
    eprintln!("  rgyr @topo system.top @traj output.trc @atoms all");
    eprintln!();
    eprintln!("  # Mass-weighted Rg for backbone atoms");
    eprintln!("  rgyr @topo system.top @traj output.trc @atoms 1-100 @massweighted");
    eprintln!();
    eprintln!("  # Protein atoms only");
    eprintln!("  rgyr @topo system.top @traj output.trc @atoms 1-1000");
}

#[derive(Debug)]
struct RgyrArgs {
    topo_file: String,
    traj_file: String,
    atom_spec: String,  // Store string, parse with AtomSelection later
    mass_weighted: bool,
}

fn parse_args(args: Vec<String>) -> Result<RgyrArgs, String> {
    let mut topo_file = None;
    let mut traj_file = None;
    let mut atom_spec = None;
    let mut mass_weighted = false;

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
            "@massweighted" => {
                mass_weighted = true;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    let topo_file = topo_file.ok_or("Missing required argument: @topo")?;
    let traj_file = traj_file.ok_or("Missing required argument: @traj")?;
    let atom_spec = atom_spec.unwrap_or_else(|| "all".to_string());

    Ok(RgyrArgs {
        topo_file,
        traj_file,
        atom_spec,  // Store as string, parse later with AtomSelection
        mass_weighted,
    })
}

/// Calculate center of mass
fn calc_center_of_mass(positions: &[Vec3], masses: &[f64], atom_selection: &[usize]) -> Vec3 {
    let mut com = Vec3::ZERO;
    let mut total_mass = 0.0;

    for &idx in atom_selection {
        com = com + positions[idx] * (masses[idx] as f32);
        total_mass += masses[idx];
    }

    if total_mass > 0.0 {
        com = com / (total_mass as f32);
    }

    com
}

/// Calculate unweighted radius of gyration
fn calc_rgyr_unweighted(positions: &[Vec3], atom_selection: &[usize]) -> f32 {
    // Calculate geometric center
    let mut center = Vec3::ZERO;
    for &idx in atom_selection {
        center = center + positions[idx];
    }
    center = center / atom_selection.len() as f32;

    // Calculate Rg²
    let mut rg_squared = 0.0;
    for &idx in atom_selection {
        let diff = positions[idx] - center;
        rg_squared += diff.dot(diff);
    }

    rg_squared /= atom_selection.len() as f32;
    rg_squared.sqrt()
}

/// Calculate mass-weighted radius of gyration
fn calc_rgyr_mass_weighted(
    positions: &[Vec3],
    masses: &[f64],
    atom_selection: &[usize],
) -> f32 {
    // Calculate center of mass
    let com = calc_center_of_mass(positions, masses, atom_selection);

    // Calculate mass-weighted Rg²
    let mut rg_squared = 0.0;
    let mut total_mass = 0.0;

    for &idx in atom_selection {
        let diff = positions[idx] - com;
        rg_squared += masses[idx] * (diff.dot(diff) as f64);
        total_mass += masses[idx];
    }

    if total_mass > 0.0 {
        rg_squared /= total_mass;
    }

    (rg_squared.sqrt() as f32)
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args[1] == "-h" || args[1] == "--help" {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let mut rgyr_args = match parse_args(args) {
        Ok(args) => args,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    // Read topology
    eprintln!("Reading topology: {}", rgyr_args.topo_file);
    let blocks = match read_topology_file(&rgyr_args.topo_file) {
        Ok(blocks) => blocks,
        Err(e) => {
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        }
    };

    let topo = build_topology(blocks);
    let n_atoms = topo.num_atoms();
    eprintln!("  Total atoms: {}", n_atoms);

    // Parse atom selection using AtomSelection (GROMOS++ compatible)
    let atom_selection = match AtomSelection::from_string(&rgyr_args.atom_spec, &topo) {
        Ok(sel) => sel,
        Err(e) => {
            eprintln!("Error parsing atom selection '{}': {}", rgyr_args.atom_spec, e);
            process::exit(1);
        }
    };

    eprintln!("  Selected atoms: {}", atom_selection.len());
    eprintln!(
        "  Mode: {}",
        if rgyr_args.mass_weighted {
            "mass-weighted"
        } else {
            "unweighted"
        }
    );
    eprintln!();

    // Open trajectory
    eprintln!("Reading trajectory: {}", rgyr_args.traj_file);
    let mut traj = match TrajectoryReader::new(&rgyr_args.traj_file) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("Error opening trajectory: {}", e);
            process::exit(1);
        }
    };

    // Output header
    println!("# Time (ps)     Rg (nm)");

    let mut frame_count = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                frame_count += 1;

                // Calculate Rg
                let rg = if rgyr_args.mass_weighted {
                    calc_rgyr_mass_weighted(
                        &frame.positions,
                        &topo.mass,
                        atom_selection.indices(),
                    )
                } else {
                    calc_rgyr_unweighted(&frame.positions, atom_selection.indices())
                };

                println!("{:12.6}{:12.6}", frame.time, rg);
            }
            Ok(None) => break, // End of trajectory
            Err(e) => {
                eprintln!("Error reading frame {}: {}", frame_count, e);
                process::exit(1);
            }
        }
    }

    eprintln!();
    eprintln!("Processed {} frames", frame_count);
}
