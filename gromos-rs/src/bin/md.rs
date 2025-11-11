//! md - Molecular Dynamics simulation engine
//!
//! Usage: md @topo <topology> @conf <coordinates> [@steps <n>] [@dt <timestep>] [@traj <output>] [@ene <energy>]
//!
//! The main GROMOS molecular dynamics simulation program.
//! Performs MD simulations with various integration algorithms.

use gromos_rs::{
    configuration::{Configuration, Box as SimBox},
    integrator::{Integrator, LeapFrog},
    math::Vec3,
    io::{
        topology::{read_topology_file, build_topology},
        trajectory::TrajectoryWriter,
        energy::{EnergyWriter, EnergyFrame},
    },
    interaction::bonded::calculate_bonded_forces,
    validation::{validate_topology, validate_coordinates, validate_configuration, validate_energy},
    logging::{LogLevel, set_log_level, start_timer, Timer},
    log_debug, log_info, log_warn, log_error,
};
use std::env;
use std::process;
use std::time::Instant;

fn print_usage() {
    eprintln!("md - Molecular Dynamics simulation");
    eprintln!();
    eprintln!("Usage: md @topo <topology> @conf <coordinates> [@steps <n>] [@dt <timestep>] [@traj <output>] [@ene <energy>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo     Topology file (.top)");
    eprintln!("  @conf     Initial coordinates (.g96)");
    eprintln!("  @steps    Number of MD steps (default: 1000)");
    eprintln!("  @dt       Time step in ps (default: 0.002)");
    eprintln!("  @traj     Output trajectory file (.trc, default: md.trc)");
    eprintln!("  @ene      Output energy file (.tre, default: md.tre)");
    eprintln!("  @temp     Temperature in K (default: 300.0)");
    eprintln!("  @verbose  Verbose output for debugging (0=normal, 1=verbose, 2=debug)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  md @topo system.top @conf initial.g96");
    eprintln!("  md @topo system.top @conf initial.g96 @steps 10000 @dt 0.002");
    eprintln!("  md @topo system.top @conf initial.g96 @traj output.trc @ene output.tre @verbose 1");
}

#[derive(Debug)]
struct MDArgs {
    topo_file: String,
    conf_file: String,
    traj_file: String,
    ene_file: String,
    n_steps: usize,
    dt: f64,
    temperature: f64,
    nstlog: usize,
    nstxout: usize,
    nstener: usize,
    verbose: usize,
}

impl Default for MDArgs {
    fn default() -> Self {
        Self {
            topo_file: String::new(),
            conf_file: String::new(),
            traj_file: "md.trc".to_string(),
            ene_file: "md.tre".to_string(),
            n_steps: 1000,
            dt: 0.002,
            temperature: 300.0,
            nstlog: 100,
            nstxout: 100,
            nstener: 10,
            verbose: 0,
        }
    }
}

fn parse_args(args: Vec<String>) -> Result<MDArgs, String> {
    let mut md_args = MDArgs::default();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @topo".to_string()); }
                md_args.topo_file = args[i].clone();
            }
            "@conf" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @conf".to_string()); }
                md_args.conf_file = args[i].clone();
            }
            "@traj" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @traj".to_string()); }
                md_args.traj_file = args[i].clone();
            }
            "@ene" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @ene".to_string()); }
                md_args.ene_file = args[i].clone();
            }
            "@steps" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @steps".to_string()); }
                md_args.n_steps = args[i].parse()
                    .map_err(|_| format!("Invalid value for @steps: {}", args[i]))?;
            }
            "@dt" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @dt".to_string()); }
                md_args.dt = args[i].parse()
                    .map_err(|_| format!("Invalid value for @dt: {}", args[i]))?;
            }
            "@temp" | "@temperature" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @temp".to_string()); }
                md_args.temperature = args[i].parse()
                    .map_err(|_| format!("Invalid value for @temp: {}", args[i]))?;
            }
            "@verbose" | "@v" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @verbose".to_string()); }
                md_args.verbose = args[i].parse()
                    .map_err(|_| format!("Invalid value for @verbose: {}", args[i]))?;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    if md_args.topo_file.is_empty() {
        return Err("Missing required argument @topo".to_string());
    }

    if md_args.conf_file.is_empty() {
        return Err("Missing required argument @conf".to_string());
    }

    Ok(md_args)
}

/// Simple coordinate file reader (reads first POSITION block from .g96)
fn read_coordinates(path: &str) -> Result<(Vec<Vec3>, Vec3), String> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path).map_err(|e| format!("Cannot open file: {}", e))?;
    let reader = BufReader::new(file);

    let mut positions = Vec::new();
    let mut box_dims = Vec3::new(3.0, 3.0, 3.0); // Default
    let mut in_position = false;
    let mut in_box = false;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Read error: {}", e))?;
        let trimmed = line.trim();

        if trimmed == "POSITION" || trimmed == "POSITIONRED" {
            in_position = true;
            continue;
        }

        if trimmed == "GENBOX" {
            in_box = true;
            in_position = false;
            continue;
        }

        if trimmed == "END" {
            if in_box {
                in_box = false;
            }
            if in_position {
                in_position = false;
            }
            continue;
        }

        if in_position {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 3 {
                // Try parsing last 3 elements as coordinates
                let len = parts.len();
                if let (Ok(x), Ok(y), Ok(z)) = (
                    parts[len-3].parse::<f32>(),
                    parts[len-2].parse::<f32>(),
                    parts[len-1].parse::<f32>(),
                ) {
                    positions.push(Vec3::new(x, y, z));
                }
            }
        }

        if in_box {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 3 {
                if let (Ok(x), Ok(y), Ok(z)) = (
                    parts[0].parse::<f32>(),
                    parts[1].parse::<f32>(),
                    parts[2].parse::<f32>(),
                ) {
                    box_dims = Vec3::new(x, y, z);
                }
            }
        }
    }

    if positions.is_empty() {
        return Err("No coordinates found in file".to_string());
    }

    Ok((positions, box_dims))
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    // Parse arguments
    let md_args = match parse_args(args) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    // Set up logging
    let log_level = match md_args.verbose {
        0 => LogLevel::Info,
        1 => LogLevel::Debug,
        _ => LogLevel::Debug,
    };
    set_log_level(log_level);
    start_timer();

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                   GROMOS-RS MD Engine                        ║");
    println!("║           Rust Implementation of GROMOS MD                   ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    log_info!("GROMOS-RS MD simulation starting");
    log_debug!("Verbose level: {}", md_args.verbose);

    // Load topology
    println!("Loading topology: {}", md_args.topo_file);
    log_debug!("Reading topology file: {}", md_args.topo_file);
    let _timer = Timer::new("Topology loading");

    let topo_data = match read_topology_file(&md_args.topo_file) {
        Ok(data) => data,
        Err(e) => {
            log_error!("Failed to read topology: {}", e);
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        }
    };

    log_debug!("Building topology data structures");
    let topo = build_topology(topo_data);

    println!("  Atoms: {}", topo.num_atoms());
    println!("  Bonds: {}", topo.solute.bonds.len());
    println!("  Angles: {}", topo.solute.angles.len());
    println!("  Dihedrals: {}", topo.solute.proper_dihedrals.len());
    println!();

    // Validate topology
    log_debug!("Validating topology");
    let topo_validation = validate_topology(&topo);
    if topo_validation.has_errors() {
        topo_validation.print();
        topo_validation.print_summary();
        if topo_validation.has_fatal() {
            log_error!("Fatal errors in topology - cannot continue");
            process::exit(1);
        }
        log_warn!("Topology has errors, but continuing");
    } else if !topo_validation.warnings.is_empty() {
        topo_validation.print();
        log_debug!("{} warnings in topology", topo_validation.warnings.len());
    } else {
        log_debug!("Topology validation passed");
    }

    // Load coordinates
    println!("Loading coordinates: {}", md_args.conf_file);
    log_debug!("Reading coordinate file: {}", md_args.conf_file);
    let _timer = Timer::new("Coordinate loading");

    let (positions, box_dims) = match read_coordinates(&md_args.conf_file) {
        Ok(data) => data,
        Err(e) => {
            log_error!("Failed to read coordinates: {}", e);
            eprintln!("Error reading coordinates: {}", e);
            process::exit(1);
        }
    };

    println!("  Positions loaded: {}", positions.len());
    println!("  Box: ({:.4}, {:.4}, {:.4}) nm", box_dims.x, box_dims.y, box_dims.z);
    println!();

    if positions.len() != topo.num_atoms() {
        log_error!("Atom count mismatch: topology={}, coordinates={}",
            topo.num_atoms(), positions.len());
        eprintln!("Error: Number of atoms in topology ({}) != coordinates ({})",
            topo.num_atoms(), positions.len());
        process::exit(1);
    }

    // Validate coordinates
    log_debug!("Validating coordinates");
    let coord_validation = validate_coordinates(&positions, Some(box_dims));
    if coord_validation.has_errors() {
        coord_validation.print();
        coord_validation.print_summary();
        if coord_validation.has_fatal() {
            log_error!("Fatal errors in coordinates - cannot continue");
            process::exit(1);
        }
        log_warn!("Coordinates have errors, but continuing");
    } else if !coord_validation.warnings.is_empty() {
        coord_validation.print();
        log_debug!("{} warnings in coordinates", coord_validation.warnings.len());
    } else {
        log_debug!("Coordinate validation passed");
    }

    // Create configuration
    log_debug!("Creating configuration");
    let mut conf = Configuration::new(topo.num_atoms(), 1, 1);
    conf.current_mut().pos = positions.clone();
    conf.current_mut().vel = vec![Vec3::ZERO; topo.num_atoms()]; // Zero initial velocities
    conf.current_mut().box_config = SimBox::rectangular(
        box_dims.x,
        box_dims.y,
        box_dims.z,
    );
    conf.copy_current_to_old();

    // Validate configuration
    log_debug!("Validating configuration (topology + coordinates)");
    let conf_validation = validate_configuration(&topo, &conf);
    if conf_validation.has_errors() {
        conf_validation.print();
        conf_validation.print_summary();
        if conf_validation.has_fatal() {
            log_error!("Fatal errors in configuration - cannot continue");
            process::exit(1);
        }
        log_warn!("Configuration has errors, but continuing");
    } else if !conf_validation.warnings.is_empty() {
        conf_validation.print();
        log_debug!("{} warnings in configuration", conf_validation.warnings.len());
    } else {
        log_debug!("Configuration validation passed");
    }

    // Setup integrator
    println!("Setting up integrator: Leap-Frog");
    let mut integrator = LeapFrog::new();
    println!();

    // Setup trajectory writer
    let mut traj_writer = match TrajectoryWriter::new(
        &md_args.traj_file,
        "GROMOS-RS MD trajectory",
        false,  // velocities
        false,  // forces
    ) {
        Ok(w) => w,
        Err(e) => {
            eprintln!("Error creating trajectory file: {}", e);
            process::exit(1);
        }
    };

    // Setup energy writer
    let mut ene_writer = match EnergyWriter::new(&md_args.ene_file, "GROMOS-RS MD energies") {
        Ok(w) => w,
        Err(e) => {
            eprintln!("Error creating energy file: {}", e);
            process::exit(1);
        }
    };

    // MD parameters summary
    println!("MD Parameters:");
    println!("  Steps:         {}", md_args.n_steps);
    println!("  Time step:     {} ps", md_args.dt);
    println!("  Total time:    {} ps", md_args.n_steps as f64 * md_args.dt);
    println!("  Temperature:   {} K", md_args.temperature);
    println!("  Traj output:   {}", md_args.traj_file);
    println!("  Energy output: {}", md_args.ene_file);
    println!();

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                   Starting MD Simulation                     ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    log_info!("Starting MD simulation: {} steps, dt={} ps", md_args.n_steps, md_args.dt);

    let start_time = Instant::now();
    let mut energy_history: Vec<(f64, f64, f64)> = Vec::new();

    // Main MD loop
    for step in 0..=md_args.n_steps {
        let time = step as f64 * md_args.dt;

        log_debug!("Step {}: time = {:.6} ps", step, time);

        // Calculate bonded forces (bonds + angles + dihedrals)
        let _force_timer = Timer::new("Force calculation");
        let bonded_result = calculate_bonded_forces(&topo, &conf, true);

        // Apply forces to configuration
        let state = conf.current_mut();
        for i in 0..topo.num_atoms() {
            state.force[i] = bonded_result.forces[i];
        }

        // Update energies
        state.energies.bond_total = bonded_result.energy;
        state.energies.update_potential_total();

        // Calculate kinetic energy
        state.calculate_kinetic_energy(&topo.mass);

        // Validate energy
        let temp = state.temperature(topo.num_atoms() * 3);
        let ene_validation = validate_energy(
            state.energies.kinetic_total,
            state.energies.potential_total,
            state.energies.total(),
            temp,
        );
        if ene_validation.has_errors() && md_args.verbose > 0 {
            ene_validation.print();
            log_warn!("Energy validation failed at step {}", step);
        }

        // Store energy for drift check
        energy_history.push((
            state.energies.kinetic_total,
            state.energies.potential_total,
            state.energies.total(),
        ));

        // Log progress
        if step % md_args.nstlog == 0 {
            println!("Step {:6}  Time: {:8.3} ps  E_pot: {:12.4}  E_kin: {:12.4}  E_tot: {:12.4}  T: {:6.1} K",
                step, time, state.energies.potential_total, state.energies.kinetic_total,
                state.energies.total(), temp);
            log_debug!("  Bond energy: {:.4}", state.energies.bond_total);
        }

        // Write trajectory
        if step % md_args.nstxout == 0 {
            if let Err(e) = traj_writer.write_frame(step, time, &conf) {
                eprintln!("Error writing trajectory: {}", e);
            }
        }

        // Write energies
        if step % md_args.nstener == 0 {
            // Get values without keeping the borrow
            let temp = {
                let state = conf.current();
                state.temperature(topo.num_atoms() * 3)
            };
            let volume = conf.current().box_config.volume();
            let pressure = conf.current().pressure();
            let energies = conf.current().energies.clone();

            let ene_frame = EnergyFrame {
                time,
                kinetic: energies.kinetic_total,
                potential: energies.potential_total,
                total: energies.total(),
                temperature: temp,
                volume,
                pressure,
                bond: energies.bond_total,
                angle: 0.0,
                improper: 0.0,
                dihedral: 0.0,
                lj: 0.0,
                coul_real: 0.0,
                coul_recip: 0.0,
                coul_self: 0.0,
                shake: 0.0,
                restraint: 0.0,
                extra: Vec::new(),
            };

            if let Err(e) = ene_writer.write_frame(&ene_frame) {
                eprintln!("Error writing energy: {}", e);
            }
        }

        // Integrate (skip last step)
        if step < md_args.n_steps {
            integrator.step(md_args.dt, &topo, &mut conf);
        }
    }

    log_info!("MD loop completed - {} steps", md_args.n_steps);

    // Check energy drift
    log_debug!("Checking energy drift over trajectory");
    use gromos_rs::validation::check_energy_drift;
    let drift_report = check_energy_drift(&energy_history);
    if drift_report.has_errors() || !drift_report.warnings.is_empty() {
        drift_report.print();
        drift_report.print_summary();
    } else {
        log_debug!("Energy drift check passed");
    }

    // Finalize output files
    log_debug!("Finalizing output files");
    if let Err(e) = traj_writer.flush() {
        log_error!("Failed to flush trajectory: {}", e);
        eprintln!("Error flushing trajectory: {}", e);
    } else {
        log_debug!("Trajectory file finalized: {}", md_args.traj_file);
    }

    if let Err(e) = ene_writer.finalize() {
        log_error!("Failed to finalize energy file: {}", e);
        eprintln!("Error finalizing energy file: {}", e);
    } else {
        log_debug!("Energy file finalized: {}", md_args.ene_file);
    }

    let elapsed = start_time.elapsed();
    log_info!("Simulation wall time: {:.2} s", elapsed.as_secs_f64());

    println!();
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                   Simulation Complete                        ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();
    println!("Statistics:");
    println!("  Total steps:     {}", md_args.n_steps);
    println!("  Simulation time: {:.3} ps", md_args.n_steps as f64 * md_args.dt);
    println!("  Wall time:       {:.2} s", elapsed.as_secs_f64());
    println!("  Performance:     {:.1} ns/day",
        (md_args.n_steps as f64 * md_args.dt * 1e-3) / elapsed.as_secs_f64() * 86400.0);
    println!();
    println!("Output files:");
    println!("  Trajectory: {}", md_args.traj_file);
    println!("  Energies:   {}", md_args.ene_file);
    println!();
    println!("Done!");
}
