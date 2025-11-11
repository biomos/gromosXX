//! hbond - Hydrogen bond analysis
//!
//! Usage: hbond @traj <trajectory> [@hbdist <cutoff>] [@hbangle <angle>]
//!
//! Monitors hydrogen bonds in molecular trajectories using geometric criteria:
//! - H-A distance < cutoff (default: 0.25 nm)
//! - D-H-A angle > min_angle (default: 135°)

use gromos_rs::{
    io::trajectory::TrajectoryReader,
    math::Vec3,
};
use std::env;
use std::process;
use std::collections::HashMap;

fn print_usage() {
    eprintln!("hbond - Hydrogen bond analysis");
    eprintln!();
    eprintln!("Usage: hbond @traj <trajectory> [@hbdist <cutoff>] [@hbangle <angle>] [@output <file>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @traj      Input trajectory file (.trc)");
    eprintln!("  @hbdist    H-A distance cutoff in nm (default: 0.25)");
    eprintln!("  @hbangle   D-H-A minimum angle in degrees (default: 135)");
    eprintln!("  @output    Output file (default: hbond.out)");
    eprintln!();
    eprintln!("Geometric criteria:");
    eprintln!("  H-bond exists if:");
    eprintln!("    1. Distance H-A < cutoff");
    eprintln!("    2. Angle D-H-A > min_angle");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  hbond @traj md.trc");
    eprintln!("  hbond @traj md.trc @hbdist 0.27 @hbangle 120");
}

#[derive(Debug)]
struct HbondArgs {
    traj_file: String,
    distance_cutoff: f32,  // nm
    angle_cutoff: f32,     // degrees
    output_file: String,
}

impl Default for HbondArgs {
    fn default() -> Self {
        Self {
            traj_file: String::new(),
            distance_cutoff: 0.25,  // 2.5 Angstrom
            angle_cutoff: 135.0,    // degrees
            output_file: "hbond.out".to_string(),
        }
    }
}

fn parse_args(args: Vec<String>) -> Result<HbondArgs, String> {
    let mut hb_args = HbondArgs::default();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @traj".to_string()); }
                hb_args.traj_file = args[i].clone();
            }
            "@hbdist" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @hbdist".to_string()); }
                hb_args.distance_cutoff = args[i].parse()
                    .map_err(|_| format!("Invalid value for @hbdist: {}", args[i]))?;
            }
            "@hbangle" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @hbangle".to_string()); }
                hb_args.angle_cutoff = args[i].parse()
                    .map_err(|_| format!("Invalid value for @hbangle: {}", args[i]))?;
            }
            "@output" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @output".to_string()); }
                hb_args.output_file = args[i].clone();
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    if hb_args.traj_file.is_empty() {
        return Err("Missing required argument @traj".to_string());
    }

    Ok(hb_args)
}

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct HydrogenBond {
    donor: usize,
    hydrogen: usize,
    acceptor: usize,
}

#[derive(Debug)]
struct HbondStats {
    count: usize,
    avg_distance: f32,
    avg_angle: f32,
    occupancy: f32,  // fraction of frames where H-bond present
}

/// Calculate angle between three points (in degrees)
fn calculate_angle(p1: Vec3, p2: Vec3, p3: Vec3) -> f32 {
    let v1 = (p1 - p2).normalize();
    let v2 = (p3 - p2).normalize();
    let dot = v1.dot(v2).clamp(-1.0, 1.0);
    dot.acos().to_degrees()
}

/// Identify potential hydrogen bonds in a frame using simple heuristics
/// This is a simplified version - a full implementation would use topology info
fn find_hydrogen_bonds(
    positions: &[Vec3],
    distance_cutoff: f32,
    angle_cutoff: f32,
) -> Vec<(HydrogenBond, f32, f32)> {
    let mut hbonds = Vec::new();
    let n_atoms = positions.len();

    // Simple heuristic: assume atoms are ordered such that hydrogens
    // are close to their donor atoms in the index
    // This is a simplified approach - real implementation needs topology

    for h in 0..n_atoms {
        // Try atoms near H as potential donors
        for d in 0..n_atoms {
            if d == h { continue; }

            let dh_dist = (positions[h] - positions[d]).length();

            // Typical D-H bond length is ~0.1 nm
            if dh_dist > 0.15 { continue; }

            // Look for acceptors
            for a in 0..n_atoms {
                if a == h || a == d { continue; }

                let ha_dist = (positions[a] - positions[h]).length();

                if ha_dist < distance_cutoff {
                    let angle = calculate_angle(positions[d], positions[h], positions[a]);

                    if angle > angle_cutoff {
                        let hbond = HydrogenBond {
                            donor: d,
                            hydrogen: h,
                            acceptor: a,
                        };
                        hbonds.push((hbond, ha_dist, angle));
                    }
                }
            }
        }
    }

    hbonds
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    // Parse arguments
    let hb_args = match parse_args(args) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    println!("hbond - Hydrogen bond analysis");
    println!("  Trajectory: {}", hb_args.traj_file);
    println!("  Distance cutoff: {:.3} nm", hb_args.distance_cutoff);
    println!("  Angle cutoff: {:.1}°", hb_args.angle_cutoff);
    println!("  Output: {}", hb_args.output_file);
    println!();

    // Open trajectory
    let mut reader = match TrajectoryReader::new(&hb_args.traj_file) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error opening trajectory: {}", e);
            process::exit(1);
        }
    };

    println!("  Trajectory title: {}", reader.title());
    println!();

    // Read all frames and analyze
    println!("  Analyzing hydrogen bonds...");

    let mut frame_count = 0;
    let mut hbond_stats: HashMap<HydrogenBond, Vec<(f32, f32)>> = HashMap::new();
    let mut total_hbonds_per_frame = Vec::new();

    loop {
        match reader.read_frame() {
            Ok(Some(frame)) => {
                frame_count += 1;

                let hbonds = find_hydrogen_bonds(
                    &frame.positions,
                    hb_args.distance_cutoff,
                    hb_args.angle_cutoff,
                );

                total_hbonds_per_frame.push(hbonds.len());

                for (hbond, dist, angle) in hbonds {
                    hbond_stats.entry(hbond)
                        .or_insert_with(Vec::new)
                        .push((dist, angle));
                }

                if frame_count % 100 == 0 {
                    println!("    Processed {} frames...", frame_count);
                }
            }
            Ok(None) => break,
            Err(e) => {
                eprintln!("Error reading frame {}: {}", frame_count, e);
                break;
            }
        }
    }

    println!("  Total frames analyzed: {}", frame_count);
    println!("  Unique H-bonds found: {}", hbond_stats.len());
    println!();

    // Calculate statistics
    let mut stats = Vec::new();
    for (hbond, data) in &hbond_stats {
        let count = data.len();
        let avg_dist = data.iter().map(|(d, _)| d).sum::<f32>() / count as f32;
        let avg_angle = data.iter().map(|(_, a)| a).sum::<f32>() / count as f32;
        let occupancy = count as f32 / frame_count as f32;

        stats.push((
            hbond.clone(),
            HbondStats {
                count,
                avg_distance: avg_dist,
                avg_angle,
                occupancy,
            },
        ));
    }

    // Sort by occupancy (most persistent bonds first)
    stats.sort_by(|a, b| b.1.occupancy.partial_cmp(&a.1.occupancy).unwrap());

    // Write output
    use std::fs::File;
    use std::io::Write;

    let mut output = File::create(&hb_args.output_file)
        .expect("Failed to create output file");

    writeln!(output, "# Hydrogen Bond Analysis").unwrap();
    writeln!(output, "# Trajectory: {}", hb_args.traj_file).unwrap();
    writeln!(output, "# Frames: {}", frame_count).unwrap();
    writeln!(output, "# Distance cutoff: {:.3} nm", hb_args.distance_cutoff).unwrap();
    writeln!(output, "# Angle cutoff: {:.1}°", hb_args.angle_cutoff).unwrap();
    writeln!(output, "#").unwrap();
    writeln!(output, "# {:>6} {:>6} {:>6} {:>8} {:>8} {:>8} {:>8}",
        "Donor", "H", "Accept", "Count", "Occupancy", "Avg_Dist", "Avg_Angle").unwrap();

    for (hbond, stat) in &stats {
        writeln!(output, "  {:6} {:6} {:6} {:8} {:8.4} {:8.4} {:8.2}",
            hbond.donor + 1,  // 1-indexed
            hbond.hydrogen + 1,
            hbond.acceptor + 1,
            stat.count,
            stat.occupancy,
            stat.avg_distance,
            stat.avg_angle
        ).unwrap();
    }

    writeln!(output, "#").unwrap();
    writeln!(output, "# Average H-bonds per frame: {:.2}",
        total_hbonds_per_frame.iter().sum::<usize>() as f32 / frame_count as f32).unwrap();

    println!("Results:");
    println!("  Total unique H-bonds: {}", stats.len());
    println!("  Average H-bonds/frame: {:.2}",
        total_hbonds_per_frame.iter().sum::<usize>() as f32 / frame_count as f32);

    if !stats.is_empty() {
        println!("\nTop 10 most persistent H-bonds:");
        println!("  {:>6} {:>6} {:>6} {:>10} {:>10}",
            "Donor", "H", "Accept", "Occupancy", "Avg_Dist");
        for (hbond, stat) in stats.iter().take(10) {
            println!("  {:6} {:6} {:6} {:9.1}% {:8.3} nm",
                hbond.donor + 1,
                hbond.hydrogen + 1,
                hbond.acceptor + 1,
                stat.occupancy * 100.0,
                stat.avg_distance
            );
        }
    }

    println!("\nOutput written to: {}", hb_args.output_file);
    println!("Done!");
}
