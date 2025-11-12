/// EDS/AEDS Molecular Dynamics Binary
///
/// Performs Enveloping Distribution Sampling (EDS) or Accelerated EDS (AEDS)
/// molecular dynamics simulations for enhanced sampling.
///
/// ## Features
/// - Single-s, multi-s, and pair-s EDS forms
/// - AEDS with adaptive parameter search
/// - Round-trip detection and statistics
/// - Free energy estimation
///
/// ## Usage
/// ```bash
/// eds @topo system.top @conf initial.g96 @states 3 @s 0.5 \
///     @offsets 0.0,10.0,20.0 @temp 300 @steps 100000 @dt 0.002
/// ```

use gromos_rs::{
    topology::Topology,
    configuration::{Configuration, Box as SimBox},
    integrator::LeapFrog,
    algorithm::thermostats::BerendsenThermostatParameters,
    algorithm::constraints::ShakeParameters,
    eds::{EDSParameters, EDSForm, EDSState, AEDSParameters, EDSRunner},
    math::Vec3,
    io::topology::{read_topology_file, build_topology},
};

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Parse command-line arguments
#[derive(Debug)]
struct EDSArgs {
    topology_file: String,
    coord_file: String,
    num_states: usize,
    s_values: Vec<f64>,
    offsets: Vec<f64>,
    temperature: f64,
    timestep: f64,
    total_steps: usize,
    output_interval: usize,
    form: EDSForm,
    use_aeds: bool,
    e_max: f64,
    e_min: f64,
    search_enabled: bool,
    cutoff: f64,
    shake_tol: f64,
    tau: f64,
}

impl Default for EDSArgs {
    fn default() -> Self {
        EDSArgs {
            topology_file: String::new(),
            coord_file: String::new(),
            num_states: 2,
            s_values: vec![0.5],
            offsets: vec![0.0, 10.0],
            temperature: 300.0,
            timestep: 0.002,
            total_steps: 10000,
            output_interval: 100,
            form: EDSForm::SingleS,
            use_aeds: false,
            e_max: 10.0,
            e_min: -50.0,
            search_enabled: false,
            cutoff: 1.4,
            shake_tol: 1e-4,
            tau: 0.1,
        }
    }
}

fn parse_args() -> Result<EDSArgs, String> {
    let args: Vec<String> = env::args().collect();
    let mut eds_args = EDSArgs::default();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing topology file after @topo".to_string());
                }
                eds_args.topology_file = args[i].clone();
            }
            "@conf" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing coordinate file after @conf".to_string());
                }
                eds_args.coord_file = args[i].clone();
            }
            "@states" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing number of states after @states".to_string());
                }
                eds_args.num_states = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid number of states: {}", args[i]))?;
            }
            "@s" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing s value(s) after @s".to_string());
                }
                eds_args.s_values = args[i]
                    .split(',')
                    .map(|s| s.trim().parse::<f64>())
                    .collect::<Result<Vec<f64>, _>>()
                    .map_err(|_| format!("Invalid s values: {}", args[i]))?;
            }
            "@offsets" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing energy offsets after @offsets".to_string());
                }
                eds_args.offsets = args[i]
                    .split(',')
                    .map(|s| s.trim().parse::<f64>())
                    .collect::<Result<Vec<f64>, _>>()
                    .map_err(|_| format!("Invalid offsets: {}", args[i]))?;
            }
            "@temp" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing temperature after @temp".to_string());
                }
                eds_args.temperature = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid temperature: {}", args[i]))?;
            }
            "@dt" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing timestep after @dt".to_string());
                }
                eds_args.timestep = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid timestep: {}", args[i]))?;
            }
            "@steps" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing number of steps after @steps".to_string());
                }
                eds_args.total_steps = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid number of steps: {}", args[i]))?;
            }
            "@output" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing output interval after @output".to_string());
                }
                eds_args.output_interval = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid output interval: {}", args[i]))?;
            }
            "@form" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing form type after @form".to_string());
                }
                eds_args.form = match args[i].to_lowercase().as_str() {
                    "single" | "singles" => EDSForm::SingleS,
                    "multi" | "multis" => EDSForm::MultiS,
                    "pair" | "pairs" => EDSForm::PairS,
                    _ => return Err(format!("Invalid form type: {}", args[i])),
                };
            }
            "@aeds" => {
                eds_args.use_aeds = true;
            }
            "@emax" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing E_max after @emax".to_string());
                }
                eds_args.e_max = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid E_max: {}", args[i]))?;
            }
            "@emin" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing E_min after @emin".to_string());
                }
                eds_args.e_min = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid E_min: {}", args[i]))?;
            }
            "@search" => {
                eds_args.search_enabled = true;
            }
            "@cutoff" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing cutoff after @cutoff".to_string());
                }
                eds_args.cutoff = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid cutoff: {}", args[i]))?;
            }
            "@tau" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing tau after @tau".to_string());
                }
                eds_args.tau = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid tau: {}", args[i]))?;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    // Validate arguments
    if eds_args.topology_file.is_empty() {
        return Err("Topology file required (@topo)".to_string());
    }
    if eds_args.coord_file.is_empty() {
        return Err("Coordinate file required (@conf)".to_string());
    }
    if eds_args.num_states < 2 {
        return Err("At least 2 states required".to_string());
    }
    if eds_args.offsets.len() != eds_args.num_states {
        return Err(format!(
            "Number of offsets ({}) must match number of states ({})",
            eds_args.offsets.len(),
            eds_args.num_states
        ));
    }

    // Validate s_values based on form
    match eds_args.form {
        EDSForm::SingleS => {
            if eds_args.s_values.len() != 1 {
                return Err("SingleS form requires exactly 1 s value".to_string());
            }
        }
        EDSForm::MultiS => {
            let expected = (eds_args.num_states * (eds_args.num_states - 1)) / 2;
            if eds_args.s_values.len() != expected {
                return Err(format!(
                    "MultiS form requires {} s values for {} states",
                    expected, eds_args.num_states
                ));
            }
        }
        EDSForm::PairS => {
            // Variable number allowed
        }
    }

    Ok(eds_args)
}

/// Read coordinates from G96 file
fn read_coordinates(filename: &str) -> Result<Vec<Vec3>, String> {
    let file = File::open(filename)
        .map_err(|e| format!("Could not open coordinate file: {}", e))?;
    let reader = BufReader::new(file);

    let mut positions = Vec::new();
    let mut in_position_block = false;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading line: {}", e))?;
        let trimmed = line.trim();

        if trimmed.starts_with("POSITION") {
            in_position_block = true;
            continue;
        }

        if trimmed.starts_with("END") {
            if in_position_block {
                break;
            }
            continue;
        }

        if in_position_block && !trimmed.is_empty() {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 6 {
                // G96 format: resnum resname atomname atomnum x y z
                let x: f32 = parts[3]
                    .parse()
                    .map_err(|_| format!("Invalid x coordinate: {}", parts[3]))?;
                let y: f32 = parts[4]
                    .parse()
                    .map_err(|_| format!("Invalid y coordinate: {}", parts[4]))?;
                let z: f32 = parts[5]
                    .parse()
                    .map_err(|_| format!("Invalid z coordinate: {}", parts[5]))?;

                positions.push(Vec3::new(x, y, z));
            }
        }
    }

    if positions.is_empty() {
        return Err("No positions found in coordinate file".to_string());
    }

    Ok(positions)
}

fn print_header(args: &EDSArgs) {
    println!("================================================================================");
    println!("                        EDS/AEDS Molecular Dynamics");
    println!("================================================================================");
    println!("Topology file:     {}", args.topology_file);
    println!("Coordinate file:   {}", args.coord_file);
    println!("Number of states:  {}", args.num_states);
    println!("EDS form:          {:?}", args.form);
    println!("Smoothing param:   {:?}", args.s_values);
    println!("Energy offsets:    {:?}", args.offsets);
    println!("Temperature:       {} K", args.temperature);
    println!("Timestep:          {} ps", args.timestep);
    println!("Total steps:       {}", args.total_steps);
    println!("Output interval:   {}", args.output_interval);

    if args.use_aeds {
        println!("Using AEDS:        yes");
        println!("E_max:             {} kJ/mol", args.e_max);
        println!("E_min:             {} kJ/mol", args.e_min);
        println!("Parameter search:  {}", if args.search_enabled { "yes" } else { "no" });
    } else {
        println!("Using AEDS:        no");
    }

    println!("================================================================================\n");
}

fn print_progress(
    step: usize,
    runner: &EDSRunner,
    config: &Configuration,
) {
    let current_state = runner.current_state();
    let ref_energy = runner.reference_energy();
    let visit_counts = runner.visit_counts();
    let round_trips = runner.round_trips();
    let total_energy = config.current().energies.total();
    let kinetic_energy = config.current().energies.kinetic_total;
    let temp = config.current().temperature(config.current().pos.len() * 3);

    println!(
        "Step {:6}  State: {}  V_R: {:10.2}  E_tot: {:10.2}  E_kin: {:8.2}  T: {:6.1} K  RT: {}",
        step, current_state, ref_energy, total_energy, kinetic_energy, temp, round_trips
    );

    // Print visit counts
    print!("  Visits: [");
    for (i, &count) in visit_counts.iter().enumerate() {
        if i > 0 {
            print!(", ");
        }
        print!("{}: {}", i, count);
    }
    println!("]");
}

fn main() -> Result<(), String> {
    // Parse arguments
    let args = parse_args()?;
    print_header(&args);

    // Load topology
    println!("Loading topology...");
    let topo_data = read_topology_file(&args.topology_file)
        .map_err(|e| format!("Failed to read topology file: {:?}", e))?;
    let topology = build_topology(topo_data);
    println!("  {} atoms loaded", topology.num_atoms());

    // Load coordinates
    println!("Loading coordinates...");
    let positions = read_coordinates(&args.coord_file)?;
    if positions.len() != topology.num_atoms() {
        return Err(format!(
            "Number of atoms in coordinate file ({}) does not match topology ({})",
            positions.len(),
            topology.num_atoms()
        ));
    }
    println!("  {} coordinates loaded", positions.len());

    // Create configuration
    println!("Initializing configuration...");
    let mut configuration = Configuration::new(topology.num_atoms(), 3, 1);
    configuration.current_mut().pos = positions;

    // Initialize velocities (zero for now - could add thermal initialization later)
    configuration.current_mut().vel = vec![Vec3::ZERO; topology.num_atoms()];

    // Set box dimensions (default 3x3x3 nm, could read from coordinate file)
    configuration.current_mut().box_config = SimBox::rectangular(3.0, 3.0, 3.0);
    configuration.copy_current_to_old();

    // Create EDS parameters
    println!("Setting up EDS...");
    let eds = EDSParameters::new(
        args.form,
        args.s_values,
        args.offsets,
        args.temperature,
        topology.num_atoms(),
    )?;

    // Create AEDS if requested
    let aeds = if args.use_aeds {
        println!("  Using AEDS with acceleration");
        AEDSParameters::new(eds, args.e_max, args.e_min, args.search_enabled)
    } else {
        println!("  Using standard EDS");
        AEDSParameters::new(eds, args.e_max, args.e_min, false)
    };

    // Create EDS runner
    let mut runner = EDSRunner::new(aeds, args.cutoff, true);

    // Setup integrator
    let mut integrator = LeapFrog::new().with_parallel();

    // Setup thermostat
    let thermostat_params = BerendsenThermostatParameters {
        target_temperature: args.temperature,
        coupling_time: args.tau,
    };

    // Setup SHAKE
    let shake_params = ShakeParameters {
        tolerance: args.shake_tol,
        max_iterations: 1000,
    };

    println!("\n================================================================================");
    println!("                           Starting EDS Simulation");
    println!("================================================================================\n");

    // Main simulation loop
    for step in 0..=args.total_steps {
        if step % args.output_interval == 0 {
            print_progress(step, &runner, &configuration);
        }

        // Perform one MD step
        runner.md_step(
            &topology,
            &mut configuration,
            &mut integrator,
            args.timestep,
            Some(&thermostat_params),
            Some(&shake_params),
        );
    }

    println!("\n================================================================================");
    println!("                           EDS Simulation Complete");
    println!("================================================================================");

    // Final statistics
    let visit_counts = runner.visit_counts();
    let round_trips = runner.round_trips();

    println!("\nFinal Statistics:");
    println!("  Total round trips: {}", round_trips);
    println!("  State visit counts:");
    for (i, &count) in visit_counts.iter().enumerate() {
        let fraction = count as f64 / args.total_steps as f64;
        println!("    State {}: {} visits ({:.2}%)", i, count, fraction * 100.0);
    }

    if args.search_enabled {
        println!("\n  Final parameters (after search):");
        println!("    s values: {:?}", runner.aeds.eds.s_values);
        println!("    Energy offsets:");
        for (i, state) in runner.aeds.eds.states.iter().enumerate() {
            println!("      State {}: {:.2} kJ/mol", i, state.offset);
        }
        println!("    Free energies:");
        for (i, &fe) in runner.aeds.free_energy.iter().enumerate() {
            println!("      State {}: {:.2} kJ/mol (relative to state 0)", i, fe);
        }
    }

    println!("\n================================================================================\n");

    Ok(())
}
