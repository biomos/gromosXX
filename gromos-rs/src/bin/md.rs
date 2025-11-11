//! md - Molecular Dynamics simulation engine
//!
//! Usage: md @topo <topology> @conf <coordinates> [@steps <n>] [@dt <timestep>] [@traj <output>] [@ene <energy>]
//!
//! The main GROMOS molecular dynamics simulation program.
//! Performs MD simulations with various integration algorithms.

use gromos_rs::{
    configuration::{Configuration, Box as SimBox},
    integrator::{Integrator, LeapFrog},
    math::{Vec3, Rectangular},
    io::{
        topology::{read_topology_file, build_topology},
        trajectory::TrajectoryWriter,
        energy::{EnergyWriter, EnergyFrame},
    },
    interaction::{
        bonded::calculate_bonded_forces,
        nonbonded::{lj_crf_innerloop, CRFParameters, ForceStorage},
    },
    pairlist::{PairlistContainer, StandardPairlistAlgorithm},
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
    eprintln!("  @topo       Topology file (.top)");
    eprintln!("  @conf       Initial coordinates (.g96)");
    eprintln!("  @steps      Number of MD steps (default: 1000)");
    eprintln!("  @dt         Time step in ps (default: 0.002)");
    eprintln!("  @traj       Output trajectory file (.trc, default: md.trc)");
    eprintln!("  @ene        Output energy file (.tre, default: md.tre)");
    eprintln!("  @temp       Temperature in K (default: 300.0)");
    eprintln!("  @cutoff     Nonbonded cutoff distance in nm (default: 1.4)");
    eprintln!("  @epsilon    Dielectric constant inside cutoff (default: 1.0)");
    eprintln!("  @rf_epsilon Dielectric constant outside cutoff for CRF (default: 61.0)");
    eprintln!("  @rf_kappa   Screening parameter for CRF in nm^-1 (default: 0.0)");
    eprintln!("  @pairlist   Pairlist update frequency in steps (default: 5)");
    eprintln!("  @verbose    Verbose output for debugging (0=normal, 1=verbose, 2=debug)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  md @topo system.top @conf initial.g96");
    eprintln!("  md @topo system.top @conf initial.g96 @steps 10000 @dt 0.002");
    eprintln!("  md @topo system.top @conf initial.g96 @cutoff 1.4 @rf_epsilon 61.0");
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
    cutoff: f64,
    epsilon: f64,
    rf_epsilon: f64,
    rf_kappa: f64,
    pairlist_update: usize,
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
            cutoff: 1.4,
            epsilon: 1.0,
            rf_epsilon: 61.0,  // Typical for water
            rf_kappa: 0.0,     // No screening
            pairlist_update: 5,
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
            "@cutoff" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @cutoff".to_string()); }
                md_args.cutoff = args[i].parse()
                    .map_err(|_| format!("Invalid value for @cutoff: {}", args[i]))?;
            }
            "@epsilon" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @epsilon".to_string()); }
                md_args.epsilon = args[i].parse()
                    .map_err(|_| format!("Invalid value for @epsilon: {}", args[i]))?;
            }
            "@rf_epsilon" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @rf_epsilon".to_string()); }
                md_args.rf_epsilon = args[i].parse()
                    .map_err(|_| format!("Invalid value for @rf_epsilon: {}", args[i]))?;
            }
            "@rf_kappa" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @rf_kappa".to_string()); }
                md_args.rf_kappa = args[i].parse()
                    .map_err(|_| format!("Invalid value for @rf_kappa: {}", args[i]))?;
            }
            "@pairlist" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @pairlist".to_string()); }
                md_args.pairlist_update = args[i].parse()
                    .map_err(|_| format!("Invalid value for @pairlist: {}", args[i]))?;
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

/// Calculate CRF parameters from physical constants
///
/// Translated from md++/src/interaction/nonbonded/interaction/cuda_nonbonded_set.cc:197-204
fn calculate_crf_parameters(
    cutoff: f64,
    epsilon: f64,
    rf_epsilon: f64,
    rf_kappa: f64,
) -> CRFParameters {
    // CRF constant calculation (GROMOS formula)
    let kappa_cut = rf_kappa * cutoff;
    let kappa_cut2 = kappa_cut * kappa_cut;

    let crf = (2.0 * (epsilon - rf_epsilon) * (1.0 + kappa_cut) - rf_epsilon * kappa_cut2) /
              ((epsilon + 2.0 * rf_epsilon) * (1.0 + kappa_cut) + rf_epsilon * kappa_cut2);

    let crf_cut3i = 1.0 / (cutoff * cutoff * cutoff);
    let crf_2cut3i = crf * crf_cut3i / 2.0;

    CRFParameters {
        crf_cut: cutoff,
        crf_2cut3i,
        crf_cut3i,
    }
}

/// Convert topology LJ parameters to nonbonded module format
///
/// Topology already contains a matrix of LJ parameters.
/// This function converts them to the format expected by nonbonded.rs
fn convert_lj_parameters(
    topo: &gromos_rs::topology::Topology
) -> Vec<Vec<gromos_rs::interaction::nonbonded::LJParameters>> {
    use gromos_rs::interaction::nonbonded::LJParameters as NBLJParams;

    // Convert topology LJ parameters to nonbonded format
    topo.lj_parameters.iter()
        .map(|row| {
            row.iter()
                .map(|params| NBLJParams {
                    c6: params.c6,
                    c12: params.c12,
                })
                .collect()
        })
        .collect()
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

    // Setup nonbonded interactions
    println!("Setting up nonbonded interactions:");
    println!("  Cutoff:      {:.3} nm", md_args.cutoff);
    println!("  Epsilon:     {:.2}", md_args.epsilon);
    println!("  RF epsilon:  {:.2}", md_args.rf_epsilon);
    println!("  RF kappa:    {:.3} nm^-1", md_args.rf_kappa);
    println!();

    log_debug!("Calculating CRF parameters");
    let crf_params = calculate_crf_parameters(
        md_args.cutoff,
        md_args.epsilon,
        md_args.rf_epsilon,
        md_args.rf_kappa,
    );
    log_debug!("CRF parameters: crf_2cut3i={:.6}, crf_cut3i={:.6}",
        crf_params.crf_2cut3i, crf_params.crf_cut3i);

    log_debug!("Converting LJ parameter matrix");
    let lj_params = convert_lj_parameters(&topo);
    log_debug!("LJ parameter matrix: {}x{} atom types", lj_params.len(),
        if lj_params.is_empty() { 0 } else { lj_params[0].len() });

    log_debug!("Initializing pairlist");
    let mut pairlist = PairlistContainer::new(
        md_args.cutoff,  // short range cutoff
        md_args.cutoff,  // long range cutoff (same for now)
        0.0,             // skin (no extra distance)
    );
    pairlist.update_frequency = md_args.pairlist_update;
    log_debug!("Pairlist update frequency: {} steps", pairlist.update_frequency);

    let pairlist_algorithm = StandardPairlistAlgorithm::new(false);  // atom-based for now
    let periodicity = Rectangular::new(box_dims);

    // Initial pairlist generation
    log_debug!("Generating initial pairlist");
    pairlist_algorithm.update(&topo, &conf, &mut pairlist, &periodicity);
    println!("  Initial pairlist: {} pairs", pairlist.total_pairs());
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

        // Update pairlist if needed
        if pairlist.needs_update() {
            log_debug!("Updating pairlist at step {}", step);
            let _pl_timer = Timer::new("Pairlist update");
            pairlist_algorithm.update(&topo, &conf, &mut pairlist, &periodicity);
            log_debug!("Pairlist updated: {} pairs", pairlist.total_pairs());
        }
        pairlist.step();

        // Calculate bonded forces (bonds + angles + dihedrals)
        let _force_timer = Timer::new("Force calculation");
        let bonded_result = calculate_bonded_forces(&topo, &conf, true);

        // Calculate nonbonded forces (LJ + CRF)
        log_debug!("Calculating nonbonded forces");
        let mut nonbonded_storage = ForceStorage::new(topo.num_atoms());

        // Convert pairlist to (u32, u32) format
        let pairlist_short: Vec<(u32, u32)> = pairlist.solute_short.iter()
            .map(|&(i, j)| (i as u32, j as u32))
            .collect();

        // Convert charge from Vec<f64> to Vec<f32> for compatibility
        let charges_f32: Vec<f32> = topo.charge.iter().map(|&q| q as f32).collect();

        // Convert iac from Vec<usize> to Vec<u32> for compatibility
        let iac_u32: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();

        lj_crf_innerloop(
            &conf.current().pos,
            &charges_f32,
            &iac_u32,
            &pairlist_short,
            &lj_params,
            &crf_params,
            &periodicity,
            &mut nonbonded_storage,
        );

        log_debug!("Nonbonded energies: LJ={:.4}, CRF={:.4}",
            nonbonded_storage.e_lj, nonbonded_storage.e_crf);

        // Apply forces to configuration (bonded + nonbonded)
        let state = conf.current_mut();
        for i in 0..topo.num_atoms() {
            // Add bonded and nonbonded forces
            state.force[i] = bonded_result.forces[i] + nonbonded_storage.forces[i];
        }

        // Update energies
        state.energies.bond_total = bonded_result.energy;
        state.energies.lj_total = nonbonded_storage.e_lj;
        state.energies.crf_total = nonbonded_storage.e_crf;
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
            log_debug!("  Bond: {:.4}  LJ: {:.4}  CRF: {:.4}",
                state.energies.bond_total, state.energies.lj_total, state.energies.crf_total);
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
                lj: energies.lj_total,
                coul_real: energies.crf_total,  // CRF energy
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
