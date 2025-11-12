//! rmsf - Root Mean Square Fluctuation Calculator
//!
//! Calculates RMSF to measure atomic flexibility and dynamics.
//! Shows which atoms/residues are most mobile during the simulation.

use gromos_rs::{
    io::topology::{read_topology_file, build_topology},
    io::trajectory::TrajectoryReader,
    selection::AtomSelection,
    math::Vec3,
    logging::{set_log_level, LogLevel, ProgressBar},
    log_info, log_debug,
};
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process;

fn print_usage() {
    eprintln!("rmsf - Root Mean Square Fluctuation Calculator");
    eprintln!();
    eprintln!("Usage:");
    eprintln!("  rmsf @topo <topology> @traj <trajectory> [@options]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo     Molecular topology file");
    eprintln!("  @traj     Trajectory file (.trc)");
    eprintln!("  @atoms    Atom selection (default: all)");
    eprintln!("            Formats: 'all', '1-10', '1,5,10-20'");
    eprintln!("                     '1:1-10' (mol 1, atoms 1-10)");
    eprintln!("                     'r:1-5' (residues 1-5)");
    eprintln!("                     'a:CA' (all CA atoms)");
    eprintln!("  @out      Output file (default: rmsf.dat)");
    eprintln!("  @skip     Skip first N frames (default: 0)");
    eprintln!("  @every    Use every Nth frame (default: 1)");
    eprintln!("  @verbose  Verbose output (0=normal, 1=verbose)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  rmsf @topo system.top @traj md.trc");
    eprintln!("  rmsf @topo system.top @traj md.trc @atoms a:CA @skip 100");
    eprintln!("  rmsf @topo system.top @traj md.trc @atoms r:1-10 @every 10");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  Column 1: Atom index");
    eprintln!("  Column 2: RMSF (nm)");
}

#[derive(Debug)]
struct RmsfArgs {
    topo_file: String,
    traj_file: String,
    atom_spec: String,
    output_file: String,
    skip: usize,
    every: usize,
    verbose: usize,
}

impl Default for RmsfArgs {
    fn default() -> Self {
        Self {
            topo_file: String::new(),
            traj_file: String::new(),
            atom_spec: "all".to_string(),
            output_file: "rmsf.dat".to_string(),
            skip: 0,
            every: 1,
            verbose: 0,
        }
    }
}

fn parse_args(args: Vec<String>) -> Result<RmsfArgs, String> {
    let mut rmsf_args = RmsfArgs::default();

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" | "@topology" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @topo".to_string()); }
                rmsf_args.topo_file = args[i].clone();
            }
            "@traj" | "@trajectory" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @traj".to_string()); }
                rmsf_args.traj_file = args[i].clone();
            }
            "@atoms" | "@atom" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @atoms".to_string()); }
                rmsf_args.atom_spec = args[i].clone();
            }
            "@out" | "@output" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @out".to_string()); }
                rmsf_args.output_file = args[i].clone();
            }
            "@skip" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @skip".to_string()); }
                rmsf_args.skip = args[i].parse()
                    .map_err(|_| format!("Invalid value for @skip: {}", args[i]))?;
            }
            "@every" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @every".to_string()); }
                rmsf_args.every = args[i].parse()
                    .map_err(|_| format!("Invalid value for @every: {}", args[i]))?;
            }
            "@verbose" | "@v" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @verbose".to_string()); }
                rmsf_args.verbose = args[i].parse()
                    .map_err(|_| format!("Invalid value for @verbose: {}", args[i]))?;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    if rmsf_args.topo_file.is_empty() {
        return Err("Missing required argument: @topo".to_string());
    }
    if rmsf_args.traj_file.is_empty() {
        return Err("Missing required argument: @traj".to_string());
    }

    Ok(rmsf_args)
}

fn calculate_rmsf(
    traj_reader: &mut TrajectoryReader,
    atom_selection: Option<&[usize]>,
    skip: usize,
    every: usize,
) -> Result<(Vec<usize>, Vec<f32>), String> {
    log_info!("Calculating RMSF");

    // Read all frames
    log_info!("Reading all frames from trajectory");
    let all_frames = traj_reader.read_all_frames()
        .map_err(|e| format!("Failed to read frames: {}", e))?;

    if all_frames.is_empty() {
        return Err("Empty trajectory".to_string());
    }

    let total_frames = all_frames.len();
    let n_atoms = all_frames[0].positions.len();
    log_info!("Total frames in trajectory: {}", total_frames);
    log_debug!("Number of atoms: {}", n_atoms);

    // Determine which atoms to analyze
    let atoms: Vec<usize> = match atom_selection {
        Some(sel) => {
            // Validate selection
            for &idx in sel {
                if idx >= n_atoms {
                    return Err(format!("Atom index {} out of range (max {})",
                        idx + 1, n_atoms));
                }
            }
            sel.to_vec()
        }
        None => (0..n_atoms).collect(),
    };

    log_info!("Analyzing {} atoms", atoms.len());
    log_info!("Skipping first {} frames, using every {} frame", skip, every);

    // First pass: Calculate average positions
    let mut avg_pos = vec![Vec3::ZERO; atoms.len()];
    let mut frames_used = 0;

    let mut progress = ProgressBar::new(total_frames);
    log_debug!("Pass 1: Calculating average positions");

    for (frame_idx, frame) in all_frames.iter().enumerate() {
        progress.update(frame_idx + 1);

        if frame_idx < skip {
            continue;
        }

        if (frame_idx - skip) % every != 0 {
            continue;
        }

        for (i, &atom_idx) in atoms.iter().enumerate() {
            avg_pos[i] = avg_pos[i] + frame.positions[atom_idx];
        }

        frames_used += 1;
    }

    if frames_used == 0 {
        return Err("No frames were processed".to_string());
    }

    log_info!("Processed {} frames", frames_used);

    // Normalize averages
    for pos in &mut avg_pos {
        *pos = *pos / frames_used as f32;
    }

    // Second pass: Calculate fluctuations
    let mut rmsf = vec![0.0f32; atoms.len()];

    progress = ProgressBar::new(total_frames);
    log_debug!("Pass 2: Calculating fluctuations");

    for (frame_idx, frame) in all_frames.iter().enumerate() {
        progress.update(frame_idx + 1);

        if frame_idx < skip {
            continue;
        }

        if (frame_idx - skip) % every != 0 {
            continue;
        }

        for (i, &atom_idx) in atoms.iter().enumerate() {
            let deviation = frame.positions[atom_idx] - avg_pos[i];
            rmsf[i] += deviation.length_squared();
        }
    }

    // Calculate RMSF = sqrt(<deviation^2>)
    for fluct in &mut rmsf {
        *fluct = (*fluct / frames_used as f32).sqrt();
    }

    Ok((atoms, rmsf))
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();

    if args.is_empty() || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        return;
    }

    let rmsf_args = match parse_args(args) {
        Ok(args) => args,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    let log_level = if rmsf_args.verbose > 0 {
        LogLevel::Debug
    } else {
        LogLevel::Info
    };
    set_log_level(log_level);

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║              RMSF - Root Mean Square Fluctuation             ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    // Read topology
    log_info!("Reading topology: {}", rmsf_args.topo_file);
    let blocks = match read_topology_file(&rmsf_args.topo_file) {
        Ok(blocks) => blocks,
        Err(e) => {
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        }
    };

    let topo = build_topology(blocks);
    log_info!("  Total atoms: {}", topo.num_atoms());

    // Parse atom selection using AtomSelection (GROMOS++ compatible)
    let atom_selection = match AtomSelection::from_string(&rmsf_args.atom_spec, &topo) {
        Ok(sel) => sel,
        Err(e) => {
            eprintln!("Error parsing atom selection '{}': {}", rmsf_args.atom_spec, e);
            process::exit(1);
        }
    };

    log_info!("  Selected atoms: {}", atom_selection.len());

    log_info!("Reading trajectory: {}", rmsf_args.traj_file);

    let mut traj_reader = match TrajectoryReader::new(&rmsf_args.traj_file) {
        Ok(reader) => reader,
        Err(e) => {
            eprintln!("Error opening trajectory: {}", e);
            process::exit(1);
        }
    };

    println!("Analyzing: {} selected atoms", atom_selection.len());
    println!("Output:    {}", rmsf_args.output_file);
    println!();

    log_info!("Calculating RMSF...");

    let (atoms, rmsf) = match calculate_rmsf(
        &mut traj_reader,
        Some(atom_selection.indices()),
        rmsf_args.skip,
        rmsf_args.every,
    ) {
        Ok(result) => result,
        Err(e) => {
            eprintln!("Error calculating RMSF: {}", e);
            process::exit(1);
        }
    };

    // Write output
    log_info!("Writing output to: {}", rmsf_args.output_file);

    let file = match File::create(&rmsf_args.output_file) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error creating output file: {}", e);
            process::exit(1);
        }
    };

    let mut writer = BufWriter::new(file);

    writeln!(writer, "# Root Mean Square Fluctuation (RMSF)").unwrap();
    writeln!(writer, "# Trajectory: {}", rmsf_args.traj_file).unwrap();
    writeln!(writer, "# Atoms analyzed: {}", atoms.len()).unwrap();
    writeln!(writer, "#").unwrap();
    writeln!(writer, "# Column 1: Atom index").unwrap();
    writeln!(writer, "# Column 2: RMSF (nm)").unwrap();
    writeln!(writer, "#").unwrap();

    for (atom_idx, fluct) in atoms.iter().zip(rmsf.iter()) {
        writeln!(writer, "{:6} {:12.6}", atom_idx + 1, fluct).unwrap();
    }

    writer.flush().unwrap();

    // Calculate statistics
    let avg_rmsf = rmsf.iter().sum::<f32>() / rmsf.len() as f32;
    let max_rmsf = rmsf.iter().cloned().fold(0.0f32, f32::max);
    let min_rmsf = rmsf.iter().cloned().fold(f32::MAX, f32::min);

    let max_idx = rmsf.iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(idx, _)| atoms[idx])
        .unwrap();

    println!();
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                     RMSF Calculation Complete                ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();
    println!("Statistics:");
    println!("  Average RMSF: {:.4} nm", avg_rmsf);
    println!("  Min RMSF:     {:.4} nm", min_rmsf);
    println!("  Max RMSF:     {:.4} nm (atom {})", max_rmsf, max_idx + 1);
    println!();
    println!("Output written to: {}", rmsf_args.output_file);
    println!();

    log_info!("RMSF calculation complete");
}
