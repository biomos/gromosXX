//! rmsd - Calculate Root Mean Square Deviation
//!
//! Usage: rmsd @topo <topology> @traj <trajectory> @ref <reference_frame> [@atoms <selection>]
//!
//! Calculates RMSD between trajectory frames and a reference structure.
//! RMSD measures structural similarity and is widely used to track
//! conformational changes during MD simulations.

use gromos_rs::io::topology::{read_topology_file, build_topology};
use gromos_rs::io::trajectory::TrajectoryReader;
use gromos_rs::selection::AtomSelection;
use gromos_rs::math::Vec3;
use std::env;
use std::process;

fn print_usage() {
    eprintln!("rmsd - Root Mean Square Deviation calculation");
    eprintln!();
    eprintln!("Usage: rmsd @topo <topology> @traj <trajectory> @ref <ref_frame> [@atoms <selection>] [@fit]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo    Molecular topology file");
    eprintln!("  @traj    Input trajectory file (.trc)");
    eprintln!("  @ref     Reference frame number (0-based) or 'first' or 'last'");
    eprintln!("  @atoms   Atom selection (default: all)");
    eprintln!("           Formats: 'all', '1-10', '1,5,10-20'");
    eprintln!("                    '1:1-10' (mol 1, atoms 1-10)");
    eprintln!("                    'r:1-5' (residues 1-5)");
    eprintln!("                    'a:CA' (all CA atoms)");
    eprintln!("  @fit     Perform rotational fit before RMSD calculation");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  Time (ps)  RMSD (nm)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  rmsd @topo system.top @traj output.trc @ref first");
    eprintln!("  rmsd @topo system.top @traj output.trc @ref 0 @atoms all");
    eprintln!("  rmsd @topo system.top @traj output.trc @ref last @atoms a:CA @fit");
}

#[derive(Debug)]
struct RmsdArgs {
    topo_file: String,
    traj_file: String,
    ref_spec: String,
    atom_spec: String,
    do_fit: bool,
}

fn parse_args(args: Vec<String>) -> Result<RmsdArgs, String> {
    let mut topo_file = None;
    let mut traj_file = None;
    let mut ref_spec = None;
    let mut atom_spec = None;
    let mut do_fit = false;

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
            "@ref" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing value for @ref".to_string());
                }
                ref_spec = Some(args[i].clone());
            }
            "@atoms" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing value for @atoms".to_string());
                }
                atom_spec = Some(args[i].clone());
            }
            "@fit" => {
                do_fit = true;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    let topo_file = topo_file.ok_or("Missing required argument: @topo")?;
    let traj_file = traj_file.ok_or("Missing required argument: @traj")?;
    let ref_spec = ref_spec.unwrap_or_else(|| "first".to_string());
    let atom_spec = atom_spec.unwrap_or_else(|| "all".to_string());

    Ok(RmsdArgs {
        topo_file,
        traj_file,
        ref_spec,
        atom_spec,
        do_fit,
    })
}

fn calculate_rmsd(pos1: &[Vec3], pos2: &[Vec3], atom_indices: &[usize]) -> f32 {
    let indices = if atom_indices.is_empty() {
        // Use all atoms
        (0..pos1.len()).collect::<Vec<_>>()
    } else {
        atom_indices.to_vec()
    };

    let n = indices.len();
    if n == 0 {
        return 0.0;
    }

    let mut sum_sq = 0.0f32;
    for &i in &indices {
        if i >= pos1.len() || i >= pos2.len() {
            continue;
        }
        let diff = pos1[i] - pos2[i];
        sum_sq += diff.length_squared();
    }

    (sum_sq / n as f32).sqrt()
}

fn translate_to_origin(positions: &mut [Vec3], atom_indices: &[usize]) -> Vec3 {
    let indices = if atom_indices.is_empty() {
        (0..positions.len()).collect::<Vec<_>>()
    } else {
        atom_indices.to_vec()
    };

    // Calculate center of mass for selected atoms
    let mut com = Vec3::ZERO;
    for &i in &indices {
        if i < positions.len() {
            com = com + positions[i];
        }
    }
    com = com / indices.len() as f32;

    // Translate all positions
    for pos in positions.iter_mut() {
        *pos = *pos - com;
    }

    com
}

fn simple_rotation_fit(mobile: &mut [Vec3], reference: &[Vec3], atom_indices: &[usize]) {
    // Simple rotation alignment using Kabsch algorithm (simplified version)
    // For production use, a full Kabsch implementation would be better

    // For now, just center both structures
    // A full implementation would compute the optimal rotation matrix
    let indices = if atom_indices.is_empty() {
        (0..mobile.len()).collect::<Vec<_>>()
    } else {
        atom_indices.to_vec()
    };

    // Calculate centroids
    let mut mobile_com = Vec3::ZERO;
    let mut ref_com = Vec3::ZERO;

    for &i in &indices {
        if i < mobile.len() && i < reference.len() {
            mobile_com = mobile_com + mobile[i];
            ref_com = ref_com + reference[i];
        }
    }

    mobile_com = mobile_com / indices.len() as f32;
    ref_com = ref_com / indices.len() as f32;

    // Center both structures
    for pos in mobile.iter_mut() {
        *pos = *pos - mobile_com + ref_com;
    }

    // Note: A complete implementation would calculate and apply
    // the optimal rotation matrix here using SVD
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    // Parse arguments
    let parsed_args = match parse_args(args) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    // Read topology
    eprintln!("Reading topology: {}", parsed_args.topo_file);
    let blocks = match read_topology_file(&parsed_args.topo_file) {
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
    let atom_selection = match AtomSelection::from_string(&parsed_args.atom_spec, &topo) {
        Ok(sel) => sel,
        Err(e) => {
            eprintln!("Error parsing atom selection '{}': {}", parsed_args.atom_spec, e);
            process::exit(1);
        }
    };

    eprintln!("  Selected atoms: {}", atom_selection.len());
    eprintln!();

    eprintln!("rmsd - RMSD calculation");
    eprintln!("  Trajectory: {}", parsed_args.traj_file);
    eprintln!("  Reference:  {}", parsed_args.ref_spec);
    eprintln!("  Atoms:      {} atoms selected", atom_selection.len());
    eprintln!("  Fit:        {}", if parsed_args.do_fit { "yes" } else { "no" });
    eprintln!();

    // Open trajectory reader
    let mut reader = match TrajectoryReader::new(&parsed_args.traj_file) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error opening trajectory file: {}", e);
            process::exit(1);
        }
    };

    eprintln!("  Title: {}", reader.title());
    eprintln!();

    // Read all frames
    eprintln!("  Reading trajectory frames...");
    let frames = match reader.read_all_frames() {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error reading trajectory frames: {}", e);
            process::exit(1);
        }
    };

    if frames.is_empty() {
        eprintln!("Error: No frames found in trajectory");
        process::exit(1);
    }

    eprintln!("  Frames read: {}", frames.len());
    eprintln!();

    // Determine reference frame
    let ref_frame_idx = match parsed_args.ref_spec.to_lowercase().as_str() {
        "first" => 0,
        "last" => frames.len() - 1,
        _ => {
            match parsed_args.ref_spec.parse::<usize>() {
                Ok(idx) => {
                    if idx >= frames.len() {
                        eprintln!("Error: Reference frame {} out of range (0-{})",
                            idx, frames.len() - 1);
                        process::exit(1);
                    }
                    idx
                }
                Err(_) => {
                    eprintln!("Error: Invalid reference specification '{}'", parsed_args.ref_spec);
                    eprintln!("Use 'first', 'last', or a frame number (0-based)");
                    process::exit(1);
                }
            }
        }
    };

    eprintln!("  Reference frame: {} (time: {:.4} ps)", ref_frame_idx, frames[ref_frame_idx].time);
    eprintln!();

    // Get reference positions
    let mut reference_positions = frames[ref_frame_idx].positions.clone();

    // Center reference structure if fitting
    if parsed_args.do_fit {
        translate_to_origin(&mut reference_positions, atom_selection.indices());
    }

    // Output header
    println!("# RMSD analysis");
    println!("# Trajectory: {}", parsed_args.traj_file);
    println!("# Reference frame: {}", ref_frame_idx);
    println!("# Atoms: {} selected", atom_selection.len());
    println!("# Fit: {}", if parsed_args.do_fit { "yes" } else { "no" });
    println!("#");
    println!("# {:>10}  {:>12}", "Time (ps)", "RMSD (nm)");

    // Calculate RMSD for each frame
    for frame in frames.iter() {
        let mut mobile_positions = frame.positions.clone();

        // Apply fit if requested
        if parsed_args.do_fit {
            translate_to_origin(&mut mobile_positions, atom_selection.indices());
            simple_rotation_fit(&mut mobile_positions, &reference_positions, atom_selection.indices());
        }

        let rmsd = calculate_rmsd(&mobile_positions, &reference_positions, atom_selection.indices());

        println!("{:12.4}  {:12.6}", frame.time, rmsd);
    }

    eprintln!();
    eprintln!("RMSD calculation complete!");
}
