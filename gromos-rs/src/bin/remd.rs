//! remd - Replica Exchange Molecular Dynamics
//!
//! Usage: remd @topo <topology> @conf <coordinates> @temps <T1,T2,...> [@steps <n>] [@dt <timestep>] ...
//!
//! Performs replica exchange MD simulations for enhanced sampling.
//! Supports temperature REMD (T-REMD) and 2D temperature-lambda REPEX.

use gromos_rs::{
    configuration::{Configuration, Box as SimBox},
    integrator::{Integrator, LeapFrog},
    math::Vec3,
    topology::Topology,
    io::{
        topology::{read_topology_file, build_topology},
    },
    replica::{Replica, ReplicaInfo},
    remd::{ReplicaController, ExchangeType, ExchangeScheme},
    algorithm::{
        thermostats::BerendsenThermostatParameters,
        constraints::ShakeParameters,
    },
    logging::{set_log_level, LogLevel},
};
use std::env;
use std::process;
use std::sync::Arc;

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

        if trimmed == "BOX" {
            in_box = true;
            continue;
        }

        if trimmed == "END" {
            in_position = false;
            in_box = false;
            continue;
        }

        if in_position {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 3 {
                if let (Ok(x), Ok(y), Ok(z)) = (
                    parts[0].parse::<f32>(),
                    parts[1].parse::<f32>(),
                    parts[2].parse::<f32>(),
                ) {
                    positions.push(Vec3::new(x, y, z));
                }
            }
        }

        if in_box && !trimmed.is_empty() {
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
        return Err("No positions found in file".to_string());
    }

    Ok((positions, box_dims))
}

fn print_usage() {
    eprintln!("remd - Replica Exchange Molecular Dynamics");
    eprintln!();
    eprintln!("Usage: remd @topo <topology> @conf <coordinates> @temps <T1,T2,...> [@steps <n>] [@dt <timestep>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo           Topology file (.top)");
    eprintln!("  @conf           Initial coordinates (.g96)");
    eprintln!("  @temps          Comma-separated temperatures in K (e.g., 280,300,320,340)");
    eprintln!("  @steps          Total MD steps per replica (default: 10000)");
    eprintln!("  @dt             Time step in ps (default: 0.002)");
    eprintln!("  @exchange_int   Exchange interval in steps (default: 100)");
    eprintln!("  @exchange_type  Exchange type: temperature, lambda, temp-lambda (default: temperature)");
    eprintln!("  @exchange_scheme Exchange scheme: sequential, odd-even, random (default: odd-even)");
    eprintln!("  @equil          Equilibration steps before exchanges (default: 1000)");
    eprintln!("  @rescale_vel    Rescale velocities after exchange: yes, no (default: yes)");
    eprintln!("  @cutoff         Nonbonded cutoff distance in nm (default: 1.4)");
    eprintln!("  @tau_t          Thermostat coupling time in ps (default: 0.1)");
    eprintln!("  @shake          SHAKE constraints: on, off (default: off)");
    eprintln!("  @shake_tol      SHAKE tolerance (default: 0.0001)");
    eprintln!("  @traj_prefix    Output trajectory prefix (default: remd)");
    eprintln!("  @ene_prefix     Output energy prefix (default: remd)");
    eprintln!("  @verbose        Verbose output: 0=normal, 1=verbose, 2=debug (default: 1)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  # Basic T-REMD with 4 replicas");
    eprintln!("  remd @topo system.top @conf initial.g96 @temps 280,300,320,340");
    eprintln!();
    eprintln!("  # T-REMD with custom parameters");
    eprintln!("  remd @topo system.top @conf initial.g96 @temps 280,300,320,340,360 \\");
    eprintln!("       @steps 50000 @exchange_int 200 @equil 5000");
    eprintln!();
    eprintln!("  # Sequential exchange scheme without velocity rescaling");
    eprintln!("  remd @topo system.top @conf initial.g96 @temps 300,310,320,330 \\");
    eprintln!("       @exchange_scheme sequential @rescale_vel no");
}

#[derive(Debug)]
struct REMDArgs {
    topo_file: String,
    conf_file: String,
    temperatures: Vec<f64>,
    n_steps: usize,
    dt: f64,
    exchange_interval: usize,
    exchange_type: ExchangeType,
    exchange_scheme: ExchangeScheme,
    equilibration_steps: usize,
    rescale_velocities: bool,
    cutoff: f64,
    tau_t: f64,
    shake_enabled: bool,
    shake_tolerance: f64,
    traj_prefix: String,
    ene_prefix: String,
    verbose: usize,
}

impl Default for REMDArgs {
    fn default() -> Self {
        REMDArgs {
            topo_file: String::new(),
            conf_file: String::new(),
            temperatures: Vec::new(),
            n_steps: 10000,
            dt: 0.002,
            exchange_interval: 100,
            exchange_type: ExchangeType::Temperature,
            exchange_scheme: ExchangeScheme::OddEven,
            equilibration_steps: 1000,
            rescale_velocities: true,
            cutoff: 1.4,
            tau_t: 0.1,
            shake_enabled: false,
            shake_tolerance: 0.0001,
            traj_prefix: "remd".to_string(),
            ene_prefix: "remd".to_string(),
            verbose: 1,
        }
    }
}

fn parse_args() -> Result<REMDArgs, String> {
    let args: Vec<String> = env::args().collect();
    let mut remd_args = REMDArgs::default();

    if args.len() < 2 {
        return Err("No arguments provided".to_string());
    }

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing topology file after @topo".to_string());
                }
                remd_args.topo_file = args[i].clone();
            }
            "@conf" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing coordinates file after @conf".to_string());
                }
                remd_args.conf_file = args[i].clone();
            }
            "@temps" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing temperatures after @temps".to_string());
                }
                // Parse comma-separated temperatures
                remd_args.temperatures = args[i]
                    .split(',')
                    .map(|s| s.trim().parse::<f64>())
                    .collect::<Result<Vec<f64>, _>>()
                    .map_err(|_| "Invalid temperature values".to_string())?;
            }
            "@steps" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing number of steps after @steps".to_string());
                }
                remd_args.n_steps = args[i]
                    .parse()
                    .map_err(|_| "Invalid number of steps".to_string())?;
            }
            "@dt" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing timestep after @dt".to_string());
                }
                remd_args.dt = args[i]
                    .parse()
                    .map_err(|_| "Invalid timestep".to_string())?;
            }
            "@exchange_int" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing exchange interval after @exchange_int".to_string());
                }
                remd_args.exchange_interval = args[i]
                    .parse()
                    .map_err(|_| "Invalid exchange interval".to_string())?;
            }
            "@exchange_type" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing exchange type after @exchange_type".to_string());
                }
                remd_args.exchange_type = match args[i].to_lowercase().as_str() {
                    "temperature" | "temp" | "t" => ExchangeType::Temperature,
                    "lambda" | "l" => ExchangeType::Lambda,
                    "temp-lambda" | "temperature-lambda" | "tl" => ExchangeType::TemperatureLambda,
                    _ => return Err(format!("Unknown exchange type: {}", args[i])),
                };
            }
            "@exchange_scheme" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing exchange scheme after @exchange_scheme".to_string());
                }
                remd_args.exchange_scheme = match args[i].to_lowercase().as_str() {
                    "sequential" | "seq" => ExchangeScheme::Sequential,
                    "odd-even" | "oddeven" | "oe" => ExchangeScheme::OddEven,
                    "random" | "rand" => ExchangeScheme::Random,
                    _ => return Err(format!("Unknown exchange scheme: {}", args[i])),
                };
            }
            "@equil" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing equilibration steps after @equil".to_string());
                }
                remd_args.equilibration_steps = args[i]
                    .parse()
                    .map_err(|_| "Invalid equilibration steps".to_string())?;
            }
            "@rescale_vel" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing rescale_vel value after @rescale_vel".to_string());
                }
                remd_args.rescale_velocities = match args[i].to_lowercase().as_str() {
                    "yes" | "on" | "true" | "1" => true,
                    "no" | "off" | "false" | "0" => false,
                    _ => return Err(format!("Invalid rescale_vel value: {}", args[i])),
                };
            }
            "@cutoff" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing cutoff after @cutoff".to_string());
                }
                remd_args.cutoff = args[i]
                    .parse()
                    .map_err(|_| "Invalid cutoff".to_string())?;
            }
            "@tau_t" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing tau_t after @tau_t".to_string());
                }
                remd_args.tau_t = args[i]
                    .parse()
                    .map_err(|_| "Invalid tau_t".to_string())?;
            }
            "@shake" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing shake value after @shake".to_string());
                }
                remd_args.shake_enabled = match args[i].to_lowercase().as_str() {
                    "on" | "yes" | "true" | "1" => true,
                    "off" | "no" | "false" | "0" => false,
                    _ => return Err(format!("Invalid shake value: {}", args[i])),
                };
            }
            "@shake_tol" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing shake tolerance after @shake_tol".to_string());
                }
                remd_args.shake_tolerance = args[i]
                    .parse()
                    .map_err(|_| "Invalid shake tolerance".to_string())?;
            }
            "@traj_prefix" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing trajectory prefix after @traj_prefix".to_string());
                }
                remd_args.traj_prefix = args[i].clone();
            }
            "@ene_prefix" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing energy prefix after @ene_prefix".to_string());
                }
                remd_args.ene_prefix = args[i].clone();
            }
            "@verbose" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing verbose level after @verbose".to_string());
                }
                remd_args.verbose = args[i]
                    .parse()
                    .map_err(|_| "Invalid verbose level".to_string())?;
            }
            arg => {
                return Err(format!("Unknown argument: {}", arg));
            }
        }
        i += 1;
    }

    // Validate required arguments
    if remd_args.topo_file.is_empty() {
        return Err("Missing required argument: @topo".to_string());
    }
    if remd_args.conf_file.is_empty() {
        return Err("Missing required argument: @conf".to_string());
    }
    if remd_args.temperatures.is_empty() {
        return Err("Missing required argument: @temps".to_string());
    }
    if remd_args.temperatures.len() < 2 {
        return Err("At least 2 temperatures required for replica exchange".to_string());
    }

    Ok(remd_args)
}

fn main() {
    // Parse command line arguments
    let args = match parse_args() {
        Ok(args) => args,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    // Set log level
    let log_level = match args.verbose {
        0 => LogLevel::Warn,
        1 => LogLevel::Info,
        _ => LogLevel::Debug,
    };
    set_log_level(log_level);

    println!("=== GROMOS-RS Replica Exchange MD ===");
    println!("Topology: {}", args.topo_file);
    println!("Coordinates: {}", args.conf_file);
    println!("Number of replicas: {}", args.temperatures.len());
    println!("Temperatures: {:?} K", args.temperatures);
    println!("Total steps: {}", args.n_steps);
    println!("Timestep: {} ps", args.dt);
    println!("Exchange interval: {} steps", args.exchange_interval);
    println!("Exchange type: {:?}", args.exchange_type);
    println!("Exchange scheme: {:?}", args.exchange_scheme);
    println!("Equilibration: {} steps", args.equilibration_steps);
    println!("Rescale velocities: {}", args.rescale_velocities);

    // Read topology
    println!("Reading topology...");
    let topo_data = match read_topology_file(&args.topo_file) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        }
    };

    let topology = Arc::new(build_topology(topo_data));
    println!("Topology loaded: {} atoms", topology.num_atoms());

    // Read initial configuration
    println!("Reading initial coordinates...");
    let (positions, box_dims) = match read_coordinates(&args.conf_file) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("Error reading coordinates: {}", e);
            process::exit(1);
        }
    };

    if positions.len() != topology.num_atoms() {
        eprintln!("Error: Number of atoms in topology ({}) != coordinates ({})",
            topology.num_atoms(), positions.len());
        process::exit(1);
    }

    // Create replicas
    println!("Creating {} replicas...", args.temperatures.len());
    let mut replicas = Vec::new();

    for (id, &temp) in args.temperatures.iter().enumerate() {
        // Clone configuration for each replica
        let mut conf = Configuration::new(
            topology.num_atoms(),
            1,
            1,
        );

        // Set positions, velocities, and box
        conf.current_mut().pos = positions.clone();
        conf.current_mut().vel = vec![Vec3::ZERO; topology.num_atoms()];
        conf.current_mut().box_config = SimBox::rectangular(box_dims.x, box_dims.y, box_dims.z);
        conf.copy_current_to_old();

        // Calculate initial kinetic energy
        conf.current_mut().calculate_kinetic_energy(&topology.mass);

        // Create replica info
        let info = ReplicaInfo::new(id, temp, 0.0, args.dt);

        // Create replica (with cutoff from args)
        let replica = Replica::new(info, conf, None, args.cutoff);
        replicas.push(replica);

        println!("  Replica {}: T = {} K", id, temp);
    }

    // Create replica controller
    println!("Initializing replica exchange controller...");
    let mut controller = ReplicaController::new(
        replicas,
        topology.clone(),
        args.exchange_type,
        args.exchange_scheme,
        args.exchange_interval,
        args.equilibration_steps,
        args.rescale_velocities,
    );

    // Setup integrator
    let integrator = LeapFrog::new().with_parallel();

    // Setup thermostat
    let thermostat_params = Some(BerendsenThermostatParameters {
        target_temperature: 0.0, // Will be overridden per replica
        coupling_time: args.tau_t,
    });

    // Setup SHAKE
    let shake_params = if args.shake_enabled {
        Some(ShakeParameters {
            tolerance: args.shake_tolerance,
            max_iterations: 1000,
        })
    } else {
        None
    };

    // Run replica exchange simulation
    println!("Starting replica exchange simulation...");
    println!("Total steps: {}", args.n_steps);
    println!("Exchange attempts: {}", args.n_steps / args.exchange_interval);

    let start_time = std::time::Instant::now();

    controller.run(
        &integrator,
        args.n_steps,
        thermostat_params.as_ref(),
        shake_params.as_ref(),
    );

    let elapsed = start_time.elapsed();

    println!("\nSimulation complete!");
    println!("Elapsed time: {:.2} s", elapsed.as_secs_f64());
    println!(
        "Performance: {:.2} ns/day",
        (args.n_steps as f64 * args.dt * 1e-3) / elapsed.as_secs_f64() * 86400.0
    );

    // Print statistics
    println!();
    controller.print_statistics();

    // Print final replica energies
    println!("\n=== Final Replica Energies ===");
    for (id, replica) in controller.replicas.iter().enumerate() {
        println!(
            "Replica {}: T = {:6.1} K, E_pot = {:12.3} kJ/mol, E_kin = {:12.3} kJ/mol, E_tot = {:12.3} kJ/mol",
            id,
            replica.info.temperature,
            replica.potential_energy(),
            replica.kinetic_energy(),
            replica.total_energy()
        );
    }

    println!("\nREMD simulation finished successfully.");
}
