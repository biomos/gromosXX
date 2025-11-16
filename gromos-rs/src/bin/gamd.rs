/// Gaussian Accelerated Molecular Dynamics (GaMD) Binary
///
/// Performs Gaussian accelerated MD simulations for enhanced sampling
/// by smoothing the potential energy landscape.
///
/// ## Features
/// - Three boost modes: dihedral, total, dual
/// - Three search phases: CMD, GaMD, production
/// - Adaptive parameter calculation
/// - Lower/upper bound threshold options
///
/// ## Usage
/// ```bash
/// # CMD search (classical MD statistics collection)
/// gamd @topo system.top @conf initial.g96 @temp 300 @steps 10000 \
///      @mode cmd_search @form total @thresh lower @sigma0 6.0
///
/// # GaMD search (accelerated with adaptive parameters)
/// gamd @topo system.top @conf initial.g96 @temp 300 @steps 50000 \
///      @mode gamd_search @form dual @thresh upper @sigma0_dih 6.0 @sigma0_tot 6.0
///
/// # Production (fixed parameters from previous search)
/// gamd @topo system.top @conf initial.g96 @temp 300 @steps 100000 \
///      @mode production @form total @k 0.001 @e 100.0
/// ```

use gromos_rs::{
    configuration::{Configuration, Box as SimBox},
    integrator::LeapFrog,
    algorithm::thermostats::BerendsenThermostatParameters,
    algorithm::constraints::ShakeParameters,
    gamd::{GamdParameters, GamdRunner, SearchMode, BoostForm, ThresholdType},
    math::Vec3,
    io::topology::{read_topology_file, build_topology},
};

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Parse command-line arguments
#[derive(Debug)]
struct GamdArgs {
    topology_file: String,
    coord_file: String,
    temperature: f64,
    timestep: f64,
    total_steps: usize,
    output_interval: usize,
    search_mode: SearchMode,
    boost_form: BoostForm,
    threshold_type: ThresholdType,
    sigma0_dih: f64,
    sigma0_tot: f64,
    k_dih: Option<f64>,
    k_tot: Option<f64>,
    e_dih: Option<f64>,
    e_tot: Option<f64>,
    equilibration_steps: usize,
    window_size: usize,
    cutoff: f64,
    shake_tol: f64,
    tau: f64,
}

impl Default for GamdArgs {
    fn default() -> Self {
        GamdArgs {
            topology_file: String::new(),
            coord_file: String::new(),
            temperature: 300.0,
            timestep: 0.002,
            total_steps: 10000,
            output_interval: 100,
            search_mode: SearchMode::CmdSearch,
            boost_form: BoostForm::TotalBoost,
            threshold_type: ThresholdType::LowerBound,
            sigma0_dih: 6.0,
            sigma0_tot: 6.0,
            k_dih: None,
            k_tot: None,
            e_dih: None,
            e_tot: None,
            equilibration_steps: 0,
            window_size: 0,
            cutoff: 1.4,
            shake_tol: 1e-4,
            tau: 0.1,
        }
    }
}

fn parse_args() -> Result<GamdArgs, String> {
    let args: Vec<String> = env::args().collect();
    let mut gamd_args = GamdArgs::default();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing topology file after @topo".to_string());
                }
                gamd_args.topology_file = args[i].clone();
            }
            "@conf" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing coordinate file after @conf".to_string());
                }
                gamd_args.coord_file = args[i].clone();
            }
            "@temp" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing temperature after @temp".to_string());
                }
                gamd_args.temperature = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid temperature: {}", args[i]))?;
            }
            "@dt" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing timestep after @dt".to_string());
                }
                gamd_args.timestep = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid timestep: {}", args[i]))?;
            }
            "@steps" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing number of steps after @steps".to_string());
                }
                gamd_args.total_steps = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid number of steps: {}", args[i]))?;
            }
            "@output" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing output interval after @output".to_string());
                }
                gamd_args.output_interval = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid output interval: {}", args[i]))?;
            }
            "@mode" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing search mode after @mode".to_string());
                }
                gamd_args.search_mode = match args[i].to_lowercase().as_str() {
                    "cmd_search" | "cmd" => SearchMode::CmdSearch,
                    "gamd_search" | "gamd" => SearchMode::GamdSearch,
                    "production" | "prod" | "no_search" => SearchMode::NoSearch,
                    _ => return Err(format!("Invalid search mode: {}", args[i])),
                };
            }
            "@form" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing boost form after @form".to_string());
                }
                gamd_args.boost_form = match args[i].to_lowercase().as_str() {
                    "dihedral" | "dih" => BoostForm::DihedralBoost,
                    "total" | "tot" => BoostForm::TotalBoost,
                    "dual" => BoostForm::DualBoost,
                    _ => return Err(format!("Invalid boost form: {}", args[i])),
                };
            }
            "@thresh" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing threshold type after @thresh".to_string());
                }
                gamd_args.threshold_type = match args[i].to_lowercase().as_str() {
                    "lower" | "lower_bound" => ThresholdType::LowerBound,
                    "upper" | "upper_bound" => ThresholdType::UpperBound,
                    _ => return Err(format!("Invalid threshold type: {}", args[i])),
                };
            }
            "@sigma0_dih" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing sigma0_dih after @sigma0_dih".to_string());
                }
                gamd_args.sigma0_dih = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid sigma0_dih: {}", args[i]))?;
            }
            "@sigma0_tot" | "@sigma0" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing sigma0_tot after @sigma0_tot".to_string());
                }
                gamd_args.sigma0_tot = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid sigma0_tot: {}", args[i]))?;
            }
            "@k_dih" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing k_dih after @k_dih".to_string());
                }
                gamd_args.k_dih = Some(
                    args[i]
                        .parse()
                        .map_err(|_| format!("Invalid k_dih: {}", args[i]))?,
                );
            }
            "@k_tot" | "@k" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing k_tot after @k_tot".to_string());
                }
                gamd_args.k_tot = Some(
                    args[i]
                        .parse()
                        .map_err(|_| format!("Invalid k_tot: {}", args[i]))?,
                );
            }
            "@e_dih" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing e_dih after @e_dih".to_string());
                }
                gamd_args.e_dih = Some(
                    args[i]
                        .parse()
                        .map_err(|_| format!("Invalid e_dih: {}", args[i]))?,
                );
            }
            "@e_tot" | "@e" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing e_tot after @e_tot".to_string());
                }
                gamd_args.e_tot = Some(
                    args[i]
                        .parse()
                        .map_err(|_| format!("Invalid e_tot: {}", args[i]))?,
                );
            }
            "@equil" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing equilibration steps after @equil".to_string());
                }
                gamd_args.equilibration_steps = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid equilibration steps: {}", args[i]))?;
            }
            "@window" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing window size after @window".to_string());
                }
                gamd_args.window_size = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid window size: {}", args[i]))?;
            }
            "@cutoff" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing cutoff after @cutoff".to_string());
                }
                gamd_args.cutoff = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid cutoff: {}", args[i]))?;
            }
            "@tau" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing tau after @tau".to_string());
                }
                gamd_args.tau = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid tau: {}", args[i]))?;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    // Validate
    if gamd_args.topology_file.is_empty() {
        return Err("Topology file required (@topo)".to_string());
    }
    if gamd_args.coord_file.is_empty() {
        return Err("Coordinate file required (@conf)".to_string());
    }

    Ok(gamd_args)
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

fn print_header(args: &GamdArgs) {
    println!("================================================================================");
    println!("                   Gaussian Accelerated Molecular Dynamics");
    println!("================================================================================");
    println!("Topology file:     {}", args.topology_file);
    println!("Coordinate file:   {}", args.coord_file);
    println!("Search mode:       {:?}", args.search_mode);
    println!("Boost form:        {:?}", args.boost_form);
    println!("Threshold type:    {:?}", args.threshold_type);
    println!("Temperature:       {} K", args.temperature);
    println!("Timestep:          {} ps", args.timestep);
    println!("Total steps:       {}", args.total_steps);
    println!("Output interval:   {}", args.output_interval);
    println!("Sigma0 (dih):      {} kJ/mol", args.sigma0_dih);
    println!("Sigma0 (tot):      {} kJ/mol", args.sigma0_tot);
    if let Some(k) = args.k_tot {
        println!("k (total):         {}", k);
    }
    if let Some(e) = args.e_tot {
        println!("E (total):         {} kJ/mol", e);
    }
    println!("================================================================================\n");
}

fn print_progress(
    step: usize,
    runner: &GamdRunner,
    config: &Configuration,
) {
    let boost = runner.boost_potential();
    let total_energy = config.current().energies.total();
    let kinetic_energy = config.current().energies.kinetic_total;
    let potential = config.current().energies.potential_total;
    let temp = config.current().temperature(config.current().pos.len() * 3);

    println!(
        "Step {:6}  V_tot: {:10.2}  V_boost: {:8.4}  E_kin: {:8.2}  E_tot: {:10.2}  T: {:6.1} K",
        step, potential - boost, boost, kinetic_energy, total_energy, temp
    );

    // Print statistics
    let tot_stats = runner.tot_stats();
    if tot_stats.steps > 0 {
        println!(
            "  Stats: V_mean={:.2}, V_min={:.2}, V_max={:.2}, σ_V={:.2} ({} samples)",
            tot_stats.v_mean, tot_stats.v_min, tot_stats.v_max, tot_stats.sigma_v, tot_stats.steps
        );

        if runner.gamd.k_tot > 0.0 {
            println!(
                "  Parameters: k0={:.4}, k={:.6}, E={:.2}",
                runner.gamd.k0_tot, runner.gamd.k_tot, runner.gamd.e_tot
            );
        }
    }
}

fn main() -> Result<(), String> {
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
    configuration.current_mut().vel = vec![Vec3::ZERO; topology.num_atoms()];
    configuration.current_mut().box_config = SimBox::rectangular(3.0, 3.0, 3.0);
    configuration.copy_current_to_old();

    // Create GaMD parameters
    println!("Setting up GaMD...");
    let mut gamd = GamdParameters::new(
        args.search_mode,
        args.boost_form,
        args.threshold_type,
        args.sigma0_dih,
        args.sigma0_tot,
    );

    // Set fixed parameters if in production mode
    if args.search_mode == SearchMode::NoSearch {
        if let Some(k) = args.k_tot {
            gamd.k_tot = k;
        }
        if let Some(e) = args.e_tot {
            gamd.e_tot = e;
        }
        if let Some(k) = args.k_dih {
            gamd.k_dih = k;
        }
        if let Some(e) = args.e_dih {
            gamd.e_dih = e;
        }
    }

    gamd.equilibration_steps = args.equilibration_steps;
    gamd.window_size = args.window_size;

    // Create GaMD runner
    let mut runner = GamdRunner::new(gamd, args.cutoff, true);

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
    println!("                          Starting GaMD Simulation");
    println!("================================================================================\n");

    // Main simulation loop
    for step in 0..=args.total_steps {
        if step % args.output_interval == 0 {
            print_progress(step, &runner, &configuration);
        }

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
    println!("                          GaMD Simulation Complete");
    println!("================================================================================");

    // Final statistics
    let tot_stats = runner.tot_stats();
    println!("\nFinal Statistics:");
    println!("  Total potential energy:");
    println!("    Mean:   {:.2} kJ/mol", tot_stats.v_mean);
    println!("    Min:    {:.2} kJ/mol", tot_stats.v_min);
    println!("    Max:    {:.2} kJ/mol", tot_stats.v_max);
    println!("    Std:    {:.2} kJ/mol", tot_stats.sigma_v);
    println!("    Samples: {}", tot_stats.steps);

    if runner.gamd.search_mode != SearchMode::NoSearch && tot_stats.is_ready() {
        println!("\n  Computed GaMD parameters:");
        println!("    k0 = {:.4}", runner.gamd.k0_tot);
        println!("    k  = {:.6}", runner.gamd.k_tot);
        println!("    E  = {:.2} kJ/mol", runner.gamd.e_tot);
        println!("\n  To use these parameters in production, run:");
        println!("    gamd ... @mode production @k {:.6} @e {:.2}", runner.gamd.k_tot, runner.gamd.e_tot);
    }

    println!("\n================================================================================\n");

    Ok(())
}
