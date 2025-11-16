//! trs_ana - Trajectory statistics analysis
//!
//! Usage: trs_ana @traj <trajectory> [@stats]
//!
//! Analyzes GROMOS trajectory files and provides statistics including:
//! - Number of atoms and frames
//! - Coordinate statistics (mean, std dev, min, max)
//! - Box dimensions analysis
//! - Trajectory quality metrics

use gromos_rs::io::trajectory::TrajectoryReader;
use gromos_rs::math::Vec3;
use std::env;
use std::process;

fn print_usage() {
    eprintln!("trs_ana - Trajectory statistics analysis");
    eprintln!();
    eprintln!("Usage: trs_ana @traj <trajectory> [@verbose]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @traj      Input trajectory file (.trc)");
    eprintln!("  @verbose   Show detailed per-frame statistics");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  trs_ana @traj output.trc");
    eprintln!("  trs_ana @traj output.trc @verbose");
}

#[derive(Debug)]
struct TrajectoryArgs {
    traj_file: String,
    verbose: bool,
}

fn parse_args(args: Vec<String>) -> Result<TrajectoryArgs, String> {
    let mut traj_file = None;
    let mut verbose = false;

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
            "@verbose" => {
                verbose = true;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    let traj_file = traj_file.ok_or("Missing required argument @traj")?;

    Ok(TrajectoryArgs {
        traj_file,
        verbose,
    })
}

#[derive(Debug)]
struct CoordinateStats {
    mean: Vec3,
    std_dev: Vec3,
    min: Vec3,
    max: Vec3,
}

impl CoordinateStats {
    fn from_positions(positions: &[Vec<Vec3>]) -> Self {
        if positions.is_empty() {
            return Self {
                mean: Vec3::ZERO,
                std_dev: Vec3::ZERO,
                min: Vec3::ZERO,
                max: Vec3::ZERO,
            };
        }

        let n_frames = positions.len();
        let n_atoms = positions[0].len();

        // Calculate mean position for each atom
        let mut mean_positions = vec![Vec3::ZERO; n_atoms];
        for frame_positions in positions {
            for (i, pos) in frame_positions.iter().enumerate() {
                mean_positions[i] = mean_positions[i] + *pos;
            }
        }
        for mean_pos in &mut mean_positions {
            *mean_pos = *mean_pos / n_frames as f32;
        }

        // Calculate overall mean
        let overall_mean = mean_positions.iter().fold(Vec3::ZERO, |acc, p| acc + *p) / n_atoms as f32;

        // Calculate standard deviation
        let mut variance = Vec3::ZERO;
        for frame_positions in positions {
            for pos in frame_positions {
                let diff = *pos - overall_mean;
                variance.x += diff.x * diff.x;
                variance.y += diff.y * diff.y;
                variance.z += diff.z * diff.z;
            }
        }
        variance = variance / (n_frames * n_atoms) as f32;
        let std_dev = Vec3::new(variance.x.sqrt(), variance.y.sqrt(), variance.z.sqrt());

        // Find min and max
        let mut min = Vec3::new(f32::INFINITY, f32::INFINITY, f32::INFINITY);
        let mut max = Vec3::new(f32::NEG_INFINITY, f32::NEG_INFINITY, f32::NEG_INFINITY);

        for frame_positions in positions {
            for pos in frame_positions {
                min.x = min.x.min(pos.x);
                min.y = min.y.min(pos.y);
                min.z = min.z.min(pos.z);
                max.x = max.x.max(pos.x);
                max.y = max.y.max(pos.y);
                max.z = max.z.max(pos.z);
            }
        }

        Self {
            mean: overall_mean,
            std_dev,
            min,
            max,
        }
    }
}

fn calculate_center_of_mass(positions: &[Vec3]) -> Vec3 {
    if positions.is_empty() {
        return Vec3::ZERO;
    }
    positions.iter().fold(Vec3::ZERO, |acc, p| acc + *p) / positions.len() as f32
}

fn calculate_radius_of_gyration(positions: &[Vec3], com: Vec3) -> f32 {
    if positions.is_empty() {
        return 0.0;
    }

    let sum_sq_dist: f32 = positions.iter()
        .map(|p| (*p - com).length_squared())
        .sum();

    (sum_sq_dist / positions.len() as f32).sqrt()
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

    println!("trs_ana - Trajectory statistics");
    println!("  Trajectory: {}", parsed_args.traj_file);
    println!();

    // Open trajectory reader
    let mut reader = match TrajectoryReader::new(&parsed_args.traj_file) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error opening trajectory file: {}", e);
            process::exit(1);
        }
    };

    println!("  Title: {}", reader.title());
    println!();

    // Read all frames
    println!("  Reading trajectory frames...");
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

    let n_frames = frames.len();
    let n_atoms = frames[0].positions.len();

    println!("  Number of frames: {}", n_frames);
    println!("  Number of atoms:  {}", n_atoms);
    println!();

    // Time range
    let times: Vec<f64> = frames.iter().map(|f| f.time).collect();
    let time_start = times.first().unwrap();
    let time_end = times.last().unwrap();
    let time_step = if n_frames > 1 {
        (time_end - time_start) / (n_frames - 1) as f64
    } else {
        0.0
    };

    println!("Time information:");
    println!("  Start time: {:.4} ps", time_start);
    println!("  End time:   {:.4} ps", time_end);
    println!("  Time step:  {:.4} ps", time_step);
    println!("  Duration:   {:.4} ps", time_end - time_start);
    println!();

    // Box dimensions statistics
    let box_dims: Vec<Vec3> = frames.iter().map(|f| f.box_dims).collect();
    let box_mean = box_dims.iter().fold(Vec3::ZERO, |acc, b| acc + *b) / n_frames as f32;

    let box_variance = box_dims.iter()
        .map(|b| {
            let diff = *b - box_mean;
            Vec3::new(diff.x * diff.x, diff.y * diff.y, diff.z * diff.z)
        })
        .fold(Vec3::ZERO, |acc, v| acc + v) / n_frames as f32;

    let box_std = Vec3::new(box_variance.x.sqrt(), box_variance.y.sqrt(), box_variance.z.sqrt());

    println!("Box dimensions:");
    println!("  Mean:   ({:.4}, {:.4}, {:.4}) nm", box_mean.x, box_mean.y, box_mean.z);
    println!("  Std:    ({:.4}, {:.4}, {:.4}) nm", box_std.x, box_std.y, box_std.z);
    println!("  Volume: {:.4} ± {:.4} nm³",
        box_mean.x * box_mean.y * box_mean.z,
        (box_std.x + box_std.y + box_std.z) / 3.0 * box_mean.x * box_mean.y * box_mean.z
    );
    println!();

    // Coordinate statistics
    let all_positions: Vec<Vec<Vec3>> = frames.iter().map(|f| f.positions.clone()).collect();
    let coord_stats = CoordinateStats::from_positions(&all_positions);

    println!("Coordinate statistics:");
    println!("  Mean:   ({:.4}, {:.4}, {:.4}) nm", coord_stats.mean.x, coord_stats.mean.y, coord_stats.mean.z);
    println!("  Std:    ({:.4}, {:.4}, {:.4}) nm", coord_stats.std_dev.x, coord_stats.std_dev.y, coord_stats.std_dev.z);
    println!("  Min:    ({:.4}, {:.4}, {:.4}) nm", coord_stats.min.x, coord_stats.min.y, coord_stats.min.z);
    println!("  Max:    ({:.4}, {:.4}, {:.4}) nm", coord_stats.max.x, coord_stats.max.y, coord_stats.max.z);
    println!("  Range:  ({:.4}, {:.4}, {:.4}) nm",
        coord_stats.max.x - coord_stats.min.x,
        coord_stats.max.y - coord_stats.min.y,
        coord_stats.max.z - coord_stats.min.z
    );
    println!();

    // Center of mass and radius of gyration analysis
    println!("Structural analysis:");

    let mut coms = Vec::new();
    let mut rgys = Vec::new();

    for frame in &frames {
        let com = calculate_center_of_mass(&frame.positions);
        let rgy = calculate_radius_of_gyration(&frame.positions, com);
        coms.push(com);
        rgys.push(rgy);
    }

    let rgy_mean = rgys.iter().sum::<f32>() / rgys.len() as f32;
    let rgy_variance = rgys.iter().map(|r| (r - rgy_mean).powi(2)).sum::<f32>() / rgys.len() as f32;
    let rgy_std = rgy_variance.sqrt();
    let rgy_min = rgys.iter().cloned().fold(f32::INFINITY, f32::min);
    let rgy_max = rgys.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

    println!("  Radius of gyration:");
    println!("    Mean:  {:.4} ± {:.4} nm", rgy_mean, rgy_std);
    println!("    Range: {:.4} - {:.4} nm", rgy_min, rgy_max);
    println!();

    // Velocity statistics (if available)
    let has_velocities = frames[0].velocities.is_some();
    if has_velocities {
        println!("Velocity information:");
        println!("  Velocities: Present in trajectory");

        let mut max_vel = 0.0f32;
        for frame in &frames {
            if let Some(ref vels) = frame.velocities {
                for vel in vels {
                    let speed = vel.length();
                    max_vel = max_vel.max(speed);
                }
            }
        }
        println!("  Maximum speed: {:.4} nm/ps", max_vel);
        println!();
    } else {
        println!("Velocity information:");
        println!("  Velocities: Not present in trajectory");
        println!();
    }

    // Force statistics (if available)
    let has_forces = frames[0].forces.is_some();
    if has_forces {
        println!("Force information:");
        println!("  Forces: Present in trajectory");

        let mut max_force = 0.0f32;
        for frame in &frames {
            if let Some(ref forces) = frame.forces {
                for force in forces {
                    let force_mag = force.length();
                    max_force = max_force.max(force_mag);
                }
            }
        }
        println!("  Maximum force: {:.4} kJ/(mol·nm)", max_force);
        println!();
    } else {
        println!("Force information:");
        println!("  Forces: Not present in trajectory");
        println!();
    }

    // Verbose mode: per-frame statistics
    if parsed_args.verbose && n_frames <= 100 {
        println!("Per-frame statistics:");
        println!("  {:>6} {:>10} {:>10} {:>10} {:>10}",
            "Frame", "Time (ps)", "Rg (nm)", "Box X", "Box Y");
        println!("  {}", "-".repeat(56));

        for (i, frame) in frames.iter().enumerate() {
            println!("  {:6} {:10.4} {:10.4} {:10.4} {:10.4}",
                i, frame.time, rgys[i], frame.box_dims.x, frame.box_dims.y);
        }
        println!();
    } else if parsed_args.verbose {
        println!("(Verbose output suppressed for trajectories > 100 frames)");
        println!();
    }

    println!("Analysis complete!");
}
