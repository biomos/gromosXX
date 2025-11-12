//! rdf - Radial Distribution Function Calculator
//!
//! Calculates the radial distribution function g(r) between atom groups.
//! Essential for analyzing liquid structure, solvation shells, coordination.

use gromos_rs::{
    io::topology::{read_topology_file, build_topology},
    io::trajectory::TrajectoryReader,
    selection::AtomSelection,
    math::Vec3,
    logging::{set_log_level, LogLevel, ProgressBar},
    log_info, log_debug,
};
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process;

fn print_usage() {
    eprintln!("rdf - Radial Distribution Function Calculator");
    eprintln!();
    eprintln!("Usage:");
    eprintln!("  rdf @topo <topology> @traj <trajectory> @group1 <atoms> @group2 <atoms> [@options]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo     Molecular topology file");
    eprintln!("  @traj     Trajectory file (.trc)");
    eprintln!("  @group1   First atom group selection");
    eprintln!("            Formats: 'all', '1-10', '1,5,10-20'");
    eprintln!("                     '1:1-10' (mol 1, atoms 1-10)");
    eprintln!("                     'r:1-5' (residues 1-5)");
    eprintln!("                     'a:OW' (all OW atoms)");
    eprintln!("  @group2   Second atom group selection (same formats)");
    eprintln!("  @out      Output file (default: rdf.dat)");
    eprintln!("  @rmax     Maximum distance in nm (default: 1.5)");
    eprintln!("  @bins     Number of bins (default: 300)");
    eprintln!("  @skip     Skip first N frames (default: 0)");
    eprintln!("  @every    Use every Nth frame (default: 1)");
    eprintln!("  @verbose  Verbose output (0=normal, 1=verbose)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  rdf @topo system.top @traj md.trc @group1 a:OW @group2 a:CA");
    eprintln!("  rdf @topo system.top @traj md.trc @group1 r:1-5 @group2 s:1 @rmax 2.0");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  Column 1: Distance r (nm)");
    eprintln!("  Column 2: g(r) (dimensionless)");
}

#[derive(Debug)]
struct RdfArgs {
    topo_file: String,
    traj_file: String,
    group1_spec: String,
    group2_spec: String,
    output_file: String,
    rmax: f32,
    bins: usize,
    skip: usize,
    every: usize,
    verbose: usize,
}

impl Default for RdfArgs {
    fn default() -> Self {
        Self {
            topo_file: String::new(),
            traj_file: String::new(),
            group1_spec: String::new(),
            group2_spec: String::new(),
            output_file: "rdf.dat".to_string(),
            rmax: 1.5,
            bins: 300,
            skip: 0,
            every: 1,
            verbose: 0,
        }
    }
}

fn parse_args(args: Vec<String>) -> Result<RdfArgs, String> {
    let mut rdf_args = RdfArgs::default();

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" | "@topology" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @topo".to_string()); }
                rdf_args.topo_file = args[i].clone();
            }
            "@traj" | "@trajectory" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @traj".to_string()); }
                rdf_args.traj_file = args[i].clone();
            }
            "@group1" | "@g1" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @group1".to_string()); }
                rdf_args.group1_spec = args[i].clone();
            }
            "@group2" | "@g2" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @group2".to_string()); }
                rdf_args.group2_spec = args[i].clone();
            }
            "@out" | "@output" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @out".to_string()); }
                rdf_args.output_file = args[i].clone();
            }
            "@rmax" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @rmax".to_string()); }
                rdf_args.rmax = args[i].parse()
                    .map_err(|_| format!("Invalid value for @rmax: {}", args[i]))?;
            }
            "@bins" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @bins".to_string()); }
                rdf_args.bins = args[i].parse()
                    .map_err(|_| format!("Invalid value for @bins: {}", args[i]))?;
            }
            "@skip" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @skip".to_string()); }
                rdf_args.skip = args[i].parse()
                    .map_err(|_| format!("Invalid value for @skip: {}", args[i]))?;
            }
            "@every" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @every".to_string()); }
                rdf_args.every = args[i].parse()
                    .map_err(|_| format!("Invalid value for @every: {}", args[i]))?;
            }
            "@verbose" | "@v" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @verbose".to_string()); }
                rdf_args.verbose = args[i].parse()
                    .map_err(|_| format!("Invalid value for @verbose: {}", args[i]))?;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    // Validate required arguments
    if rdf_args.topo_file.is_empty() {
        return Err("Missing required argument: @topo".to_string());
    }
    if rdf_args.traj_file.is_empty() {
        return Err("Missing required argument: @traj".to_string());
    }
    if rdf_args.group1_spec.is_empty() {
        return Err("Missing required argument: @group1".to_string());
    }
    if rdf_args.group2_spec.is_empty() {
        return Err("Missing required argument: @group2".to_string());
    }

    Ok(rdf_args)
}

fn calculate_rdf(
    traj_reader: &mut TrajectoryReader,
    group1: &[usize],
    group2: &[usize],
    rmax: f32,
    bins: usize,
    skip: usize,
    every: usize,
) -> Result<(Vec<f32>, Vec<f64>), String> {
    let dr = rmax / bins as f32;
    let mut histogram = vec![0u64; bins];

    log_info!("Calculating RDF with {} bins, dr = {:.6} nm", bins, dr);
    log_debug!("Group 1: {} atoms", group1.len());
    log_debug!("Group 2: {} atoms", group2.len());

    // Read all frames
    log_info!("Reading all frames from trajectory");
    let all_frames = traj_reader.read_all_frames()
        .map_err(|e| format!("Failed to read frames: {}", e))?;

    let total_frames = all_frames.len();
    log_info!("Total frames in trajectory: {}", total_frames);
    log_info!("Skipping first {} frames, using every {} frame", skip, every);

    let mut progress = ProgressBar::new(total_frames);
    let mut frames_used = 0;
    let mut avg_volume = 0.0f32;

    // Process frames
    for (frame_idx, frame) in all_frames.iter().enumerate() {
        progress.update(frame_idx + 1);

        // Skip frames
        if frame_idx < skip {
            continue;
        }

        // Use every Nth frame
        if (frame_idx - skip) % every != 0 {
            continue;
        }

        log_debug!("Processing frame {} (time = {:.3} ps)", frame_idx, frame.time);

        // Get box volume
        let box_volume = frame.box_dims.x * frame.box_dims.y * frame.box_dims.z;
        avg_volume += box_volume;

        // Calculate distances between all pairs
        for &i in group1 {
            if i >= frame.positions.len() {
                return Err(format!("Atom index {} out of range (max {})",
                    i + 1, frame.positions.len()));
            }

            for &j in group2 {
                if j >= frame.positions.len() {
                    return Err(format!("Atom index {} out of range (max {})",
                        j + 1, frame.positions.len()));
                }

                // Skip self-pairs
                if i == j {
                    continue;
                }

                let r_ij = frame.positions[j] - frame.positions[i];
                let dist = r_ij.length();

                if dist < rmax {
                    let bin = (dist / dr) as usize;
                    if bin < bins {
                        histogram[bin] += 1;
                    }
                }
            }
        }

        frames_used += 1;
    }

    log_info!("Processed {} frames", frames_used);

    if frames_used == 0 {
        return Err("No frames were processed".to_string());
    }

    // Calculate average volume
    avg_volume /= frames_used as f32;
    log_debug!("Average box volume: {:.3} nm³", avg_volume);

    // Calculate density (number of particles per unit volume)
    let n1 = group1.len() as f64;
    let n2 = group2.len() as f64;
    let density = n2 / avg_volume as f64;

    log_debug!("Number density: {:.6} particles/nm³", density);

    // Normalize histogram to g(r)
    let mut r_values = Vec::with_capacity(bins);
    let mut g_r = Vec::with_capacity(bins);

    for i in 0..bins {
        let r = (i as f32 + 0.5) * dr;
        r_values.push(r);

        // Volume of spherical shell
        let r_inner = i as f32 * dr;
        let r_outer = (i + 1) as f32 * dr;
        let shell_volume = (4.0 / 3.0) * std::f32::consts::PI *
            (r_outer.powi(3) - r_inner.powi(3));

        // Expected number of particles in ideal gas
        let n_ideal = density * shell_volume as f64;

        // Normalize
        let n_observed = histogram[i] as f64 / (frames_used as f64 * n1);
        let g = if n_ideal > 0.0 {
            n_observed / n_ideal
        } else {
            0.0
        };

        g_r.push(g);
    }

    Ok((r_values, g_r))
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();

    if args.is_empty() || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        return;
    }

    let rdf_args = match parse_args(args) {
        Ok(args) => args,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    // Set up logging
    let log_level = if rdf_args.verbose > 0 {
        LogLevel::Debug
    } else {
        LogLevel::Info
    };
    set_log_level(log_level);

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║              RDF - Radial Distribution Function             ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    // Read topology
    log_info!("Reading topology: {}", rdf_args.topo_file);
    let blocks = match read_topology_file(&rdf_args.topo_file) {
        Ok(blocks) => blocks,
        Err(e) => {
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        }
    };

    let topo = build_topology(blocks);
    log_info!("  Total atoms: {}", topo.num_atoms());

    // Parse atom selections using AtomSelection (GROMOS++ compatible)
    let group1 = match AtomSelection::from_string(&rdf_args.group1_spec, &topo) {
        Ok(sel) => sel,
        Err(e) => {
            eprintln!("Error parsing group1 selection '{}': {}", rdf_args.group1_spec, e);
            process::exit(1);
        }
    };

    let group2 = match AtomSelection::from_string(&rdf_args.group2_spec, &topo) {
        Ok(sel) => sel,
        Err(e) => {
            eprintln!("Error parsing group2 selection '{}': {}", rdf_args.group2_spec, e);
            process::exit(1);
        }
    };

    log_info!("  Group 1: {} atoms", group1.len());
    log_info!("  Group 2: {} atoms", group2.len());

    log_info!("Reading trajectory: {}", rdf_args.traj_file);

    let mut traj_reader = match TrajectoryReader::new(&rdf_args.traj_file) {
        Ok(reader) => reader,
        Err(e) => {
            eprintln!("Error opening trajectory: {}", e);
            process::exit(1);
        }
    };

    println!("Group 1: {} atoms", group1.len());
    println!("Group 2: {} atoms", group2.len());
    println!("r_max:   {:.3} nm", rdf_args.rmax);
    println!("Bins:    {}", rdf_args.bins);
    println!("Output:  {}", rdf_args.output_file);
    println!();

    log_info!("Calculating RDF...");

    let (r_values, g_r) = match calculate_rdf(
        &mut traj_reader,
        group1.indices(),
        group2.indices(),
        rdf_args.rmax,
        rdf_args.bins,
        rdf_args.skip,
        rdf_args.every,
    ) {
        Ok(result) => result,
        Err(e) => {
            eprintln!("Error calculating RDF: {}", e);
            process::exit(1);
        }
    };

    // Write output
    log_info!("Writing output to: {}", rdf_args.output_file);

    let file = match File::create(&rdf_args.output_file) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error creating output file: {}", e);
            process::exit(1);
        }
    };

    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(writer, "# Radial Distribution Function g(r)").unwrap();
    writeln!(writer, "# Trajectory: {}", rdf_args.traj_file).unwrap();
    writeln!(writer, "# Group 1: {} atoms", group1.len()).unwrap();
    writeln!(writer, "# Group 2: {} atoms", group2.len()).unwrap();
    writeln!(writer, "# r_max: {:.3} nm", rdf_args.rmax).unwrap();
    writeln!(writer, "# Bins: {}", rdf_args.bins).unwrap();
    writeln!(writer, "#").unwrap();
    writeln!(writer, "# Column 1: r (nm)").unwrap();
    writeln!(writer, "# Column 2: g(r)").unwrap();
    writeln!(writer, "#").unwrap();

    // Write data
    for (r, g) in r_values.iter().zip(g_r.iter()) {
        writeln!(writer, "{:12.6} {:12.6}", r, g).unwrap();
    }

    writer.flush().unwrap();

    // Find first peak
    let mut max_g = 0.0;
    let mut max_r = 0.0;
    for (r, g) in r_values.iter().zip(g_r.iter()) {
        if *g > max_g {
            max_g = *g;
            max_r = *r;
        }
    }

    println!();
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                      RDF Calculation Complete                ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();
    println!("First peak: g({:.3} nm) = {:.3}", max_r, max_g);
    println!("Output written to: {}", rdf_args.output_file);
    println!();

    log_info!("RDF calculation complete");
}
