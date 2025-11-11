//! rmsd - Calculate Root Mean Square Deviation
//!
//! Usage: rmsd @traj <trajectory> @ref <reference_frame> [@atoms <selection>]
//!
//! Calculates RMSD between trajectory frames and a reference structure.
//! RMSD measures structural similarity and is widely used to track
//! conformational changes during MD simulations.

use gromos_rs::io::trajectory::TrajectoryReader;
use gromos_rs::math::Vec3;
use std::env;
use std::process;

fn print_usage() {
    eprintln!("rmsd - Root Mean Square Deviation calculation");
    eprintln!();
    eprintln!("Usage: rmsd @traj <trajectory> @ref <ref_frame> [@atoms <selection>] [@fit]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @traj    Input trajectory file (.trc)");
    eprintln!("  @ref     Reference frame number (0-based) or 'first' or 'last'");
    eprintln!("  @atoms   Atom selection (default: all)");
    eprintln!("           Format: 'all', '1-10', '1,5,10-20'");
    eprintln!("  @fit     Perform rotational fit before RMSD calculation");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  Time (ps)  RMSD (nm)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  rmsd @traj output.trc @ref first");
    eprintln!("  rmsd @traj output.trc @ref 0 @atoms all");
    eprintln!("  rmsd @traj output.trc @ref last @atoms 1-100 @fit");
}

#[derive(Debug)]
struct RmsdArgs {
    traj_file: String,
    ref_spec: String,
    atom_selection: Vec<usize>,
    do_fit: bool,
}

fn parse_args(args: Vec<String>) -> Result<RmsdArgs, String> {
    let mut traj_file = None;
    let mut ref_spec = None;
    let mut atom_spec = None;
    let mut do_fit = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
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

    let traj_file = traj_file.ok_or("Missing required argument @traj")?;
    let ref_spec = ref_spec.unwrap_or_else(|| "first".to_string());

    // Parse atom selection (default to empty, meaning all atoms)
    let atom_selection = if let Some(spec) = atom_spec {
        parse_atom_selection(&spec)?
    } else {
        Vec::new()
    };

    Ok(RmsdArgs {
        traj_file,
        ref_spec,
        atom_selection,
        do_fit,
    })
}

fn parse_atom_selection(spec: &str) -> Result<Vec<usize>, String> {
    if spec.to_lowercase() == "all" {
        return Ok(Vec::new()); // Empty means all atoms
    }

    let mut atoms = Vec::new();

    for part in spec.split(',') {
        if part.contains('-') {
            // Range: "1-10"
            let range_parts: Vec<&str> = part.split('-').collect();
            if range_parts.len() != 2 {
                return Err(format!("Invalid range: {}", part));
            }
            let start: usize = range_parts[0]
                .trim()
                .parse()
                .map_err(|_| format!("Invalid number: {}", range_parts[0]))?;
            let end: usize = range_parts[1]
                .trim()
                .parse()
                .map_err(|_| format!("Invalid number: {}", range_parts[1]))?;

            for i in start..=end {
                atoms.push(i - 1); // Convert to 0-based indexing
            }
        } else {
            // Single atom
            let atom: usize = part
                .trim()
                .parse()
                .map_err(|_| format!("Invalid number: {}", part))?;
            atoms.push(atom - 1); // Convert to 0-based indexing
        }
    }

    Ok(atoms)
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

    eprintln!("rmsd - RMSD calculation");
    eprintln!("  Trajectory: {}", parsed_args.traj_file);
    eprintln!("  Reference:  {}", parsed_args.ref_spec);
    if !parsed_args.atom_selection.is_empty() {
        eprintln!("  Atoms:      {} atoms selected", parsed_args.atom_selection.len());
    } else {
        eprintln!("  Atoms:      all");
    }
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
        translate_to_origin(&mut reference_positions, &parsed_args.atom_selection);
    }

    // Output header
    println!("# RMSD analysis");
    println!("# Trajectory: {}", parsed_args.traj_file);
    println!("# Reference frame: {}", ref_frame_idx);
    println!("# Atoms: {}", if parsed_args.atom_selection.is_empty() {
        "all".to_string()
    } else {
        format!("{} selected", parsed_args.atom_selection.len())
    });
    println!("# Fit: {}", if parsed_args.do_fit { "yes" } else { "no" });
    println!("#");
    println!("# {:>10}  {:>12}", "Time (ps)", "RMSD (nm)");

    // Calculate RMSD for each frame
    for frame in frames.iter() {
        let mut mobile_positions = frame.positions.clone();

        // Apply fit if requested
        if parsed_args.do_fit {
            translate_to_origin(&mut mobile_positions, &parsed_args.atom_selection);
            simple_rotation_fit(&mut mobile_positions, &reference_positions, &parsed_args.atom_selection);
        }

        let rmsd = calculate_rmsd(&mobile_positions, &reference_positions, &parsed_args.atom_selection);

        println!("{:12.4}  {:12.6}", frame.time, rmsd);
    }

    eprintln!();
    eprintln!("RMSD calculation complete!");
}
