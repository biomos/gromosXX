//! md - Molecular Dynamics simulation engine
//!
//! Usage: md @topo <topology> @conf <coordinates> [@steps <n>] [@dt <timestep>] [@traj <output>] [@ene <energy>]
//!
//! The main GROMOS molecular dynamics simulation program.
//! Performs MD simulations with various integration algorithms.

use gromos_rs::{
    configuration::{Configuration, Box as SimBox},
    integrator::{Integrator, LeapFrog},
    math::{Vec3, Rectangular, Mat3},
    io::{
        topology::{read_topology_file, build_topology},
        trajectory::TrajectoryWriter,
        energy::{EnergyWriter, EnergyFrame},
        GamdBlock, EdsBlock,
    },
    gamd::GamdParameters,
    eds::EDSParameters,
    interaction::{
        bonded::calculate_bonded_forces,
        nonbonded::{lj_crf_innerloop, CRFParameters, ForceStorage},
    },
    pairlist::{PairlistContainer, StandardPairlistAlgorithm},
    algorithm::{
        shake, ShakeParameters,
        berendsen_thermostat, BerendsenThermostatParameters,
        berendsen_barostat, BerendsenBarostatParameters,
    },
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
    eprintln!("Usage: md @topo <topology> @conf <coordinates> [@input <parameters>] [@steps <n>] [@dt <timestep>] [@traj <output>] [@ene <energy>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo       Topology file (.top)");
    eprintln!("  @conf       Initial coordinates (.g96)");
    eprintln!("  @input      Input parameter file (.imd) for GAMD/EDS/advanced sampling");
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
    eprintln!("  @thermostat Thermostat: off, berendsen (default: off)");
    eprintln!("  @tau_t      Thermostat coupling time in ps (default: 0.1)");
    eprintln!("  @barostat   Barostat: off, berendsen (default: off)");
    eprintln!("  @tau_p      Barostat coupling time in ps (default: 0.5)");
    eprintln!("  @pres       Target pressure in bar (default: 1.0)");
    eprintln!("  @shake      SHAKE constraints: on, off (default: off)");
    eprintln!("  @shake_tol  SHAKE tolerance (default: 0.0001)");
    eprintln!("  @verbose    Verbose output for debugging (0=normal, 1=verbose, 2=debug)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  md @topo system.top @conf initial.g96");
    eprintln!("  md @topo system.top @conf initial.g96 @steps 10000 @dt 0.002");
    eprintln!("  md @topo system.top @conf initial.g96 @cutoff 1.4 @rf_epsilon 61.0");
    eprintln!("  md @topo system.top @conf initial.g96 @thermostat berendsen @temp 300");
    eprintln!("  md @topo system.top @conf initial.g96 @barostat berendsen @pres 1.0");
    eprintln!("  md @topo system.top @conf initial.g96 @shake on @shake_tol 0.0001");
    eprintln!("  md @topo system.top @conf initial.g96 @traj output.trc @ene output.tre @verbose 1");
}

#[derive(Debug)]
struct MDArgs {
    topo_file: String,
    conf_file: String,
    input_file: Option<String>,  // Optional .imd input file for GAMD/EDS/etc.
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
    // Thermostat
    thermostat: String,      // "off", "berendsen", "nose-hoover"
    tau_t: f64,              // Thermostat coupling time
    // Barostat
    barostat: String,        // "off", "berendsen", "parrinello-rahman"
    tau_p: f64,              // Barostat coupling time
    target_pressure: f64,    // Target pressure in bar
    // Constraints
    shake_enabled: bool,
    shake_tolerance: f64,
    // Output
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
            input_file: None,
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
            // Thermostat
            thermostat: "off".to_string(),
            tau_t: 0.1,        // 0.1 ps (typical)
            // Barostat
            barostat: "off".to_string(),
            tau_p: 0.5,        // 0.5 ps (typical)
            target_pressure: 1.0,  // 1 bar (atmospheric)
            // Constraints
            shake_enabled: false,
            shake_tolerance: 1e-4,  // GROMOS default
            // Output
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
            "@input" | "@imd" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @input".to_string()); }
                md_args.input_file = Some(args[i].clone());
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
            "@thermostat" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @thermostat".to_string()); }
                md_args.thermostat = args[i].to_lowercase();
                if !["off", "berendsen", "nose-hoover"].contains(&md_args.thermostat.as_str()) {
                    return Err(format!("Invalid thermostat: {}. Use: off, berendsen, nose-hoover", args[i]));
                }
            }
            "@tau_t" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @tau_t".to_string()); }
                md_args.tau_t = args[i].parse()
                    .map_err(|_| format!("Invalid value for @tau_t: {}", args[i]))?;
            }
            "@barostat" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @barostat".to_string()); }
                md_args.barostat = args[i].to_lowercase();
                if !["off", "berendsen", "parrinello-rahman"].contains(&md_args.barostat.as_str()) {
                    return Err(format!("Invalid barostat: {}. Use: off, berendsen, parrinello-rahman", args[i]));
                }
            }
            "@tau_p" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @tau_p".to_string()); }
                md_args.tau_p = args[i].parse()
                    .map_err(|_| format!("Invalid value for @tau_p: {}", args[i]))?;
            }
            "@pres" | "@pressure" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @pres".to_string()); }
                md_args.target_pressure = args[i].parse()
                    .map_err(|_| format!("Invalid value for @pres: {}", args[i]))?;
            }
            "@shake" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @shake".to_string()); }
                let shake_str = args[i].to_lowercase();
                md_args.shake_enabled = match shake_str.as_str() {
                    "on" | "yes" | "true" | "1" => true,
                    "off" | "no" | "false" | "0" => false,
                    _ => return Err(format!("Invalid value for @shake: {}. Use: on, off", args[i])),
                };
            }
            "@shake_tol" | "@shake_tolerance" => {
                i += 1;
                if i >= args.len() { return Err("Missing value for @shake_tol".to_string()); }
                md_args.shake_tolerance = args[i].parse()
                    .map_err(|_| format!("Invalid value for @shake_tol: {}", args[i]))?;
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

    // Parse input file for GAMD/EDS blocks if provided
    let gamd_block = if let Some(ref input_file) = md_args.input_file {
        log_debug!("Parsing input file for GAMD block: {}", input_file);
        match GamdBlock::parse_file(input_file) {
            Ok(block) => block,
            Err(e) => {
                log_warn!("Failed to parse GAMD block: {}", e);
                None
            }
        }
    } else {
        None
    };

    let eds_block = if let Some(ref input_file) = md_args.input_file {
        log_debug!("Parsing input file for EDS block: {}", input_file);
        match EdsBlock::parse_file(input_file) {
            Ok(block) => block,
            Err(e) => {
                log_warn!("Failed to parse EDS block: {}", e);
                None
            }
        }
    } else {
        None
    };

    if let Some(ref block) = gamd_block {
        println!();
        println!("GAMD Parameters detected:");
        println!("  Search mode:  {:?}", block.search_mode);
        println!("  Boost form:   {:?}", block.boost_form);
        println!("  Threshold:    {:?}", block.threshold_type);
        println!("  Sigma0 dih:   {:.2}", block.sigma0_dih);
        println!("  Sigma0 tot:   {:.2}", block.sigma0_tot);
        if let (Some(k), Some(e)) = (block.k_tot, block.e_tot) {
            println!("  K_tot:        {:.6}", k);
            println!("  E_tot:        {:.2}", e);
        }
        println!();
    }

    if let Some(ref block) = eds_block {
        println!();
        println!("EDS Parameters detected:");
        println!("  Num states:   {}", block.num_states);
        println!("  Form:         {:?}", block.form);
        println!("  S values:     {:?}", block.s_values);
        println!("  E offsets:    {:?}", block.e_offsets);
        println!("  Temperature:  {:.1} K", block.temperature);
        if block.search_enabled {
            println!("  AEDS enabled: E_max={:.2}, E_min={:.2}", block.e_max, block.e_min);
        }
        println!();
    }

    // Check for conflicting modes
    if gamd_block.is_some() && eds_block.is_some() {
        log_error!("Cannot enable both GAMD and EDS simultaneously");
        eprintln!("Error: Both GAMD and EDS blocks found in input file");
        eprintln!("       These methods cannot be used together in the same simulation");
        process::exit(1);
    }

    // Create GAMD parameters if enabled
    let mut gamd_params = gamd_block.as_ref().map(|block| {
        log_info!("Creating GAMD parameters from input block");
        block.to_parameters()
    });

    // Create EDS parameters if enabled (will be done after topology is fully loaded)
    let eds_params = eds_block.as_ref().map(|block| {
        log_info!("Creating EDS parameters from input block");
        // Note: We need num_atoms, which we'll get after creating configuration
        block
    });

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

    // Setup thermostat
    let thermostat_params = if md_args.thermostat == "berendsen" {
        Some(BerendsenThermostatParameters {
            target_temperature: md_args.temperature,
            coupling_time: md_args.tau_t,
        })
    } else {
        None
    };

    if let Some(ref params) = thermostat_params {
        println!("Setting up thermostat: Berendsen");
        println!("  Target temp:   {:.1} K", params.target_temperature);
        println!("  Coupling time: {:.3} ps", params.coupling_time);
        println!();
    }

    // Setup barostat
    let barostat_params = if md_args.barostat == "berendsen" {
        Some(BerendsenBarostatParameters {
            target_pressure: md_args.target_pressure,
            coupling_time: md_args.tau_p,
            compressibility: 4.5e-5,  // Water compressibility
            isotropic: true,
        })
    } else {
        None
    };

    if let Some(ref params) = barostat_params {
        println!("Setting up barostat: Berendsen");
        println!("  Target pres:   {:.1} bar", params.target_pressure);
        println!("  Coupling time: {:.3} ps", params.coupling_time);
        println!();
    }

    // Setup SHAKE constraints
    let shake_params = if md_args.shake_enabled {
        Some(ShakeParameters {
            tolerance: md_args.shake_tolerance,
            max_iterations: 1000,
        })
    } else {
        None
    };

    if let Some(ref params) = shake_params {
        println!("Setting up constraints: SHAKE");
        println!("  Tolerance:     {:.6}", params.tolerance);
        println!("  Max iter:      {}", params.max_iterations);
        println!();
    }

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

        // Virial tensor for barostat (will be updated by nonbonded forces)
        let mut virial = Mat3::ZERO;

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

        // Update virial for barostat
        virial = Mat3 {
            x_axis: Vec3::new(
                nonbonded_storage.virial[0][0] as f32,
                nonbonded_storage.virial[0][1] as f32,
                nonbonded_storage.virial[0][2] as f32,
            ),
            y_axis: Vec3::new(
                nonbonded_storage.virial[1][0] as f32,
                nonbonded_storage.virial[1][1] as f32,
                nonbonded_storage.virial[1][2] as f32,
            ),
            z_axis: Vec3::new(
                nonbonded_storage.virial[2][0] as f32,
                nonbonded_storage.virial[2][1] as f32,
                nonbonded_storage.virial[2][2] as f32,
            ),
        };

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

            // Apply SHAKE constraints to satisfy bond length constraints
            if let Some(ref params) = shake_params {
                log_debug!("Applying SHAKE constraints");
                let _shake_timer = Timer::new("SHAKE constraints");
                let shake_result = shake(&topo, &mut conf, md_args.dt, params);

                if !shake_result.converged {
                    log_warn!("SHAKE did not converge at step {}: {} iterations, error={:.6}",
                        step, shake_result.iterations, shake_result.max_error);
                }

                log_debug!("SHAKE converged in {} iterations, error={:.6}",
                    shake_result.iterations, shake_result.max_error);
            }

            // Apply thermostat to control temperature
            if let Some(ref params) = thermostat_params {
                log_debug!("Applying Berendsen thermostat");
                let _thermo_timer = Timer::new("Thermostat");
                berendsen_thermostat(&topo, &mut conf, md_args.dt, params);
            }

            // Apply barostat to control pressure
            if let Some(ref params) = barostat_params {
                log_debug!("Applying Berendsen barostat");
                let _baro_timer = Timer::new("Barostat");

                // Use virial from nonbonded calculations
                berendsen_barostat(&topo, &mut conf, md_args.dt, params, &virial);
            }
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
