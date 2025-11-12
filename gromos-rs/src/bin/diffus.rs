//! diffus - Diffusion Coefficient Calculator
//!
//! Calculates diffusion coefficients from mean square displacement (MSD).
//! Uses the Einstein relation: D = lim(t->∞) <|r(t)-r(0)|²> / (6t)

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
    eprintln!("diffus - Diffusion Coefficient Calculator");
    eprintln!();
    eprintln!("Usage:");
    eprintln!("  diffus @topo <topology> @traj <trajectory> [@options]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo     Molecular topology file");
    eprintln!("  @traj     Trajectory file (.trc)");
    eprintln!("  @atoms    Atom selection (default: all)");
    eprintln!("            Formats: 'all', '1-10', '1,5,10-20'");
    eprintln!("                     '1:1-10' (mol 1, atoms 1-10)");
    eprintln!("                     'r:1-5' (residues 1-5)");
    eprintln!("                     'a:OW' (all OW atoms)");
    eprintln!("  @out      Output file for MSD vs time (default: msd.dat)");
    eprintln!("  @skip     Skip first N frames (default: 0)");
    eprintln!("  @every    Use every Nth frame (default: 1)");
    eprintln!("  @fitstart Start of linear fit region in ps (default: auto)");
    eprintln!("  @fitend   End of linear fit region in ps (default: auto)");
    eprintln!("  @verbose  Verbose output (0=normal, 1=verbose)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  diffus @topo system.top @traj md.trc");
    eprintln!("  diffus @topo system.top @traj md.trc @atoms a:OW @skip 1000");
    eprintln!("  diffus @topo system.top @traj md.trc @fitstart 50 @fitend 200");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  MSD file - Column 1: Time (ps), Column 2: MSD (nm²)");
    eprintln!("  Diffusion coefficient printed to stdout");
}

#[derive(Debug)]
struct DiffusArgs {
    topo_file: String,
    traj_file: String,
    atom_spec: String,
    output_file: String,
    skip: usize,
    every: usize,
    fit_start: Option<f64>,
    fit_end: Option<f64>,
    verbose: usize,
}

impl Default for DiffusArgs {
    fn default() -> Self {
        Self {
            topo_file: String::new(),
            traj_file: String::new(),
            atom_spec: "all".to_string(),
            output_file: "msd.dat".to_string(),
            skip: 0,
            every: 1,
            fit_start: None,
            fit_end: None,
            verbose: 0,
        }
    }
}

fn parse_args(args: Vec<String>) -> Result<DiffusArgs, String> {
    let mut diffus_args = DiffusArgs::default();

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" | "@topology" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @topo".to_string()); }
                diffus_args.topo_file = args[i].clone();
            }
            "@traj" | "@trajectory" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @traj".to_string()); }
                diffus_args.traj_file = args[i].clone();
            }
            "@atoms" | "@atom" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @atoms".to_string()); }
                diffus_args.atom_spec = args[i].clone();
            }
            "@out" | "@output" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @out".to_string()); }
                diffus_args.output_file = args[i].clone();
            }
            "@skip" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @skip".to_string()); }
                diffus_args.skip = args[i].parse()
                    .map_err(|_| format!("Invalid value for @skip: {}", args[i]))?;
            }
            "@every" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @every".to_string()); }
                diffus_args.every = args[i].parse()
                    .map_err(|_| format!("Invalid value for @every: {}", args[i]))?;
            }
            "@fitstart" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @fitstart".to_string()); }
                diffus_args.fit_start = Some(args[i].parse()
                    .map_err(|_| format!("Invalid value for @fitstart: {}", args[i]))?);
            }
            "@fitend" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @fitend".to_string()); }
                diffus_args.fit_end = Some(args[i].parse()
                    .map_err(|_| format!("Invalid value for @fitend: {}", args[i]))?);
            }
            "@verbose" | "@v" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @verbose".to_string()); }
                diffus_args.verbose = args[i].parse()
                    .map_err(|_| format!("Invalid value for @verbose: {}", args[i]))?;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    if diffus_args.topo_file.is_empty() {
        return Err("Missing required argument: @topo".to_string());
    }
    if diffus_args.traj_file.is_empty() {
        return Err("Missing required argument: @traj".to_string());
    }

    Ok(diffus_args)
}

fn calculate_msd(
    traj_reader: &mut TrajectoryReader,
    atom_selection: Option<&[usize]>,
    skip: usize,
    every: usize,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    log_info!("Calculating Mean Square Displacement (MSD)");

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

    let atoms: Vec<usize> = match atom_selection {
        Some(sel) => {
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

    // Store selected positions and times
    let mut all_positions: Vec<Vec<Vec3>> = Vec::new();
    let mut times: Vec<f64> = Vec::new();

    let mut progress = ProgressBar::new(total_frames);
    log_debug!("Extracting selected atom positions");

    for (frame_idx, frame) in all_frames.iter().enumerate() {
        progress.update(frame_idx + 1);

        if frame_idx < skip {
            continue;
        }

        if (frame_idx - skip) % every != 0 {
            continue;
        }

        let positions: Vec<Vec3> = atoms.iter()
            .map(|&idx| frame.positions[idx])
            .collect();

        all_positions.push(positions);
        times.push(frame.time);
    }

    let n_frames = all_positions.len();
    if n_frames < 2 {
        return Err("Need at least 2 frames to calculate MSD".to_string());
    }

    log_info!("Loaded {} frames", n_frames);

    // Calculate MSD for each time interval
    let mut msd_times: Vec<f64> = Vec::new();
    let mut msd_values: Vec<f64> = Vec::new();

    log_debug!("Calculating MSD for all time intervals");
    let mut progress = ProgressBar::new(n_frames);

    for dt_idx in 0..n_frames {
        progress.update(dt_idx + 1);

        let dt = times[dt_idx] - times[0];
        let mut msd = 0.0f64;
        let mut n_samples = 0;

        // Average over all time origins
        for t0 in 0..(n_frames - dt_idx) {
            let t1 = t0 + dt_idx;

            // Average over all selected atoms
            for atom_local_idx in 0..atoms.len() {
                let r0 = all_positions[t0][atom_local_idx];
                let r1 = all_positions[t1][atom_local_idx];
                let dr = r1 - r0;
                msd += dr.length_squared() as f64;
                n_samples += 1;
            }
        }

        msd /= n_samples as f64;

        msd_times.push(dt);
        msd_values.push(msd);

        log_debug!("MSD at t={:.3} ps: {:.6} nm²", dt, msd);
    }

    Ok((msd_times, msd_values))
}

fn linear_fit(x: &[f64], y: &[f64], start_idx: usize, end_idx: usize) -> (f64, f64) {
    let n = end_idx - start_idx + 1;
    if n < 2 {
        return (0.0, 0.0);
    }

    let x_slice = &x[start_idx..=end_idx];
    let y_slice = &y[start_idx..=end_idx];

    let sum_x: f64 = x_slice.iter().sum();
    let sum_y: f64 = y_slice.iter().sum();
    let sum_xx: f64 = x_slice.iter().map(|&xi| xi * xi).sum();
    let sum_xy: f64 = x_slice.iter().zip(y_slice.iter())
        .map(|(&xi, &yi)| xi * yi).sum();

    let n_f = n as f64;
    let slope = (n_f * sum_xy - sum_x * sum_y) / (n_f * sum_xx - sum_x * sum_x);
    let intercept = (sum_y - slope * sum_x) / n_f;

    (slope, intercept)
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();

    if args.is_empty() || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        return;
    }

    let diffus_args = match parse_args(args) {
        Ok(args) => args,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    let log_level = if diffus_args.verbose > 0 {
        LogLevel::Debug
    } else {
        LogLevel::Info
    };
    set_log_level(log_level);

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║              Diffusion Coefficient Calculator                ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    // Read topology
    log_info!("Reading topology: {}", diffus_args.topo_file);
    let blocks = match read_topology_file(&diffus_args.topo_file) {
        Ok(blocks) => blocks,
        Err(e) => {
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        }
    };

    let topo = build_topology(blocks);
    log_info!("  Total atoms: {}", topo.num_atoms());

    // Parse atom selection using AtomSelection (GROMOS++ compatible)
    let atom_selection = match AtomSelection::from_string(&diffus_args.atom_spec, &topo) {
        Ok(sel) => sel,
        Err(e) => {
            eprintln!("Error parsing atom selection '{}': {}", diffus_args.atom_spec, e);
            process::exit(1);
        }
    };

    log_info!("  Selected atoms: {}", atom_selection.len());

    log_info!("Reading trajectory: {}", diffus_args.traj_file);

    let mut traj_reader = match TrajectoryReader::new(&diffus_args.traj_file) {
        Ok(reader) => reader,
        Err(e) => {
            eprintln!("Error opening trajectory: {}", e);
            process::exit(1);
        }
    };

    println!("Analyzing: {} selected atoms", atom_selection.len());
    println!("Output:    {}", diffus_args.output_file);
    println!();

    log_info!("Calculating MSD...");

    let (times, msd) = match calculate_msd(
        &mut traj_reader,
        Some(atom_selection.indices()),
        diffus_args.skip,
        diffus_args.every,
    ) {
        Ok(result) => result,
        Err(e) => {
            eprintln!("Error calculating MSD: {}", e);
            process::exit(1);
        }
    };

    // Write MSD output
    log_info!("Writing MSD to: {}", diffus_args.output_file);

    let file = match File::create(&diffus_args.output_file) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error creating output file: {}", e);
            process::exit(1);
        }
    };

    let mut writer = BufWriter::new(file);

    writeln!(writer, "# Mean Square Displacement (MSD)").unwrap();
    writeln!(writer, "# Trajectory: {}", diffus_args.traj_file).unwrap();
    writeln!(writer, "#").unwrap();
    writeln!(writer, "# Column 1: Time (ps)").unwrap();
    writeln!(writer, "# Column 2: MSD (nm²)").unwrap();
    writeln!(writer, "#").unwrap();

    for (t, m) in times.iter().zip(msd.iter()) {
        writeln!(writer, "{:12.4} {:12.6}", t, m).unwrap();
    }

    writer.flush().unwrap();

    // Determine fit region
    let max_time = times[times.len() - 1];
    let fit_start_time = diffus_args.fit_start.unwrap_or(max_time * 0.5);
    let fit_end_time = diffus_args.fit_end.unwrap_or(max_time * 0.9);

    let start_idx = times.iter()
        .position(|&t| t >= fit_start_time)
        .unwrap_or(times.len() / 2);
    let end_idx = times.iter()
        .position(|&t| t >= fit_end_time)
        .unwrap_or(times.len() - 1);

    // Linear fit: MSD = 6Dt => D = slope/6
    let (slope, intercept) = linear_fit(&times, &msd, start_idx, end_idx);
    let diffusion_coeff = slope / 6.0;

    // Convert to commonly used units
    let d_nm2_ps = diffusion_coeff;  // nm²/ps
    let d_cm2_s = diffusion_coeff * 1e-3;  // cm²/s (multiply by 10^-3 to convert nm²/ps to cm²/s)
    let d_m2_s = diffusion_coeff * 1e-15;  // m²/s

    println!();
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                 Diffusion Coefficient Results                ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();
    println!("Linear fit region: {:.2} - {:.2} ps", times[start_idx], times[end_idx]);
    println!("MSD = {:.6} + {:.6} * t", intercept, slope);
    println!();
    println!("Diffusion coefficient D:");
    println!("  {:.6e} nm²/ps", d_nm2_ps);
    println!("  {:.6e} cm²/s", d_cm2_s);
    println!("  {:.6e} m²/s", d_m2_s);
    println!();
    println!("MSD data written to: {}", diffus_args.output_file);
    println!();

    log_info!("Diffusion calculation complete");
}
