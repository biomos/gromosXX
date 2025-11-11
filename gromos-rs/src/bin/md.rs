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
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  md @topo system.top @conf initial.g96");
    eprintln!("  md @topo system.top @conf initial.g96 @steps 10000 @dt 0.002");
    eprintln!("  md @topo system.top @conf initial.g96 @traj output.trc @ene output.tre");
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

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                   GROMOS-RS MD Engine                        ║");
    println!("║           Rust Implementation of GROMOS MD                   ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    // Load topology
    println!("Loading topology: {}", md_args.topo_file);
    let topo_data = match read_topology_file(&md_args.topo_file) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        }
    };

    let topo = build_topology(topo_data);

    println!("  Atoms: {}", topo.num_atoms());
    println!("  Bonds: {}", topo.solute.bonds.len());
    println!("  Angles: {}", topo.solute.angles.len());
    println!("  Dihedrals: {}", topo.solute.proper_dihedrals.len());
    println!();

    // Load coordinates
    println!("Loading coordinates: {}", md_args.conf_file);
    let (positions, box_dims) = match read_coordinates(&md_args.conf_file) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("Error reading coordinates: {}", e);
            process::exit(1);
        }
    };

    println!("  Positions loaded: {}", positions.len());
    println!("  Box: ({:.4}, {:.4}, {:.4}) nm", box_dims.x, box_dims.y, box_dims.z);
    println!();

    if positions.len() != topo.num_atoms() {
        eprintln!("Error: Number of atoms in topology ({}) != coordinates ({})",
            topo.num_atoms(), positions.len());
        process::exit(1);
    }

    // Create configuration
    let mut conf = Configuration::new(topo.num_atoms(), 1, 1);
    conf.current_mut().pos = positions.clone();
    conf.current_mut().vel = vec![Vec3::ZERO; topo.num_atoms()]; // Zero initial velocities
    conf.current_mut().box_config = SimBox::rectangular(
        box_dims.x,
        box_dims.y,
        box_dims.z,
    );
    conf.copy_current_to_old();

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

    let start_time = Instant::now();

    // Main MD loop
    for step in 0..=md_args.n_steps {
        let time = step as f64 * md_args.dt;

        // Calculate bonded forces (bonds + angles + dihedrals)
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

        // Log progress
        if step % md_args.nstlog == 0 {
            let temp = state.temperature(topo.num_atoms() * 3);
            println!("Step {:6}  Time: {:8.3} ps  E_pot: {:12.4}  E_kin: {:12.4}  E_tot: {:12.4}  T: {:6.1} K",
                step, time, state.energies.potential_total, state.energies.kinetic_total,
                state.energies.total(), temp);
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

    // Finalize output files
    if let Err(e) = traj_writer.flush() {
        eprintln!("Error flushing trajectory: {}", e);
    }

    if let Err(e) = ene_writer.finalize() {
        eprintln!("Error finalizing energy file: {}", e);
    }

    let elapsed = start_time.elapsed();

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
