//! make_pt_top - Create perturbation topology from two molecular topologies
//!
//! Usage: make_pt_top @topo <top_a> <top_b> [@softpar <alpha_lj> <alpha_crf>] [@out <output.ptp>]
//!
//! Creates a perturbation topology file (.ptp) for free energy calculations
//! by comparing two molecular topologies (states A and B).

use gromos_rs::io::{
    topology::{read_topology_file, build_topology},
    PtpWriter,
};
use std::env;
use std::process;

fn print_usage() {
    eprintln!("make_pt_top - Create perturbation topology for FEP calculations");
    eprintln!();
    eprintln!("Usage: make_pt_top @topo <topology_a> <topology_b> [@softpar <alpha_lj> <alpha_crf>] [@out <output.ptp>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo       Two topology files: state A (λ=0) and state B (λ=1)");
    eprintln!("  @softpar    Softness parameters (default: 1.51 0.5)");
    eprintln!("              alpha_lj:  Lennard-Jones softness (typ. 1.51)");
    eprintln!("              alpha_crf: Coulomb reaction field softness (typ. 0.5)");
    eprintln!("  @out        Output .ptp file (default: perturbation.ptp)");
    eprintln!();
    eprintln!("Description:");
    eprintln!("  Compares two molecular topologies and generates a perturbation topology");
    eprintln!("  file for free energy perturbation (FEP) calculations. Both topologies");
    eprintln!("  must have the same number of atoms.");
    eprintln!();
    eprintln!("  The tool identifies atoms with different:");
    eprintln!("  - Integer atom codes (IAC / atom types)");
    eprintln!("  - Masses");
    eprintln!("  - Charges");
    eprintln!();
    eprintln!("  Softness parameters control the soft-core potential:");
    eprintln!("  - alpha_lj:  ~1.5 for smooth LJ transitions");
    eprintln!("  - alpha_crf: ~0.5 for smooth electrostatic transitions");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  # Create perturbation topology with default softness");
    eprintln!("  make_pt_top @topo stateA.top stateB.top");
    eprintln!();
    eprintln!("  # Specify softness parameters");
    eprintln!("  make_pt_top @topo stateA.top stateB.top @softpar 1.51 0.5");
    eprintln!();
    eprintln!("  # Custom output file");
    eprintln!("  make_pt_top @topo stateA.top stateB.top @out custom.ptp");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  A .ptp file containing:");
    eprintln!("  - PERTURBEDATOM block: List of atoms that change between states");
    eprintln!("  - PERTURBATIONPARAMETERS block: Softness parameters");
    eprintln!();
    eprintln!("Reference:");
    eprintln!("  GROMOS manual section on free energy perturbation");
    eprintln!("  Liu & Berne, J. Chem. Phys. 118, 2977 (2003) - soft-core potentials");
}

#[derive(Debug)]
struct MakePtTopArgs {
    topo_a: String,
    topo_b: String,
    alpha_lj: f64,
    alpha_crf: f64,
    output: String,
}

impl Default for MakePtTopArgs {
    fn default() -> Self {
        Self {
            topo_a: String::new(),
            topo_b: String::new(),
            alpha_lj: 1.51,
            alpha_crf: 0.5,
            output: "perturbation.ptp".to_string(),
        }
    }
}

fn parse_args(args: Vec<String>) -> Result<MakePtTopArgs, String> {
    let mut pt_args = MakePtTopArgs::default();
    let mut topo_files = Vec::new();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                // Read two topology files
                i += 1;
                if i >= args.len() {
                    return Err("Missing topology file A for @topo".to_string());
                }
                topo_files.push(args[i].clone());

                i += 1;
                if i >= args.len() {
                    return Err("Missing topology file B for @topo".to_string());
                }
                topo_files.push(args[i].clone());
            }
            "@softpar" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing alpha_lj for @softpar".to_string());
                }
                pt_args.alpha_lj = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid alpha_lj: {}", args[i]))?;

                i += 1;
                if i >= args.len() {
                    return Err("Missing alpha_crf for @softpar".to_string());
                }
                pt_args.alpha_crf = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid alpha_crf: {}", args[i]))?;
            }
            "@out" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing output file for @out".to_string());
                }
                pt_args.output = args[i].clone();
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    if topo_files.len() != 2 {
        return Err("Must provide exactly two topology files with @topo".to_string());
    }

    pt_args.topo_a = topo_files[0].clone();
    pt_args.topo_b = topo_files[1].clone();

    Ok(pt_args)
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args[1] == "-h" || args[1] == "--help" {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let pt_args = match parse_args(args) {
        Ok(args) => args,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    println!("make_pt_top - Creating perturbation topology");
    println!("  State A:    {}", pt_args.topo_a);
    println!("  State B:    {}", pt_args.topo_b);
    println!("  alpha_LJ:   {}", pt_args.alpha_lj);
    println!("  alpha_CRF:  {}", pt_args.alpha_crf);
    println!("  Output:     {}", pt_args.output);
    println!();

    // Read topology A
    println!("Reading state A topology...");
    let blocks_a = match read_topology_file(&pt_args.topo_a) {
        Ok(blocks) => blocks,
        Err(e) => {
            eprintln!("Error reading topology A: {}", e);
            process::exit(1);
        }
    };

    let topo_a = build_topology(blocks_a);
    println!("  State A: {} atoms", topo_a.num_atoms());

    // Read topology B
    println!("Reading state B topology...");
    let blocks_b = match read_topology_file(&pt_args.topo_b) {
        Ok(blocks) => blocks,
        Err(e) => {
            eprintln!("Error reading topology B: {}", e);
            process::exit(1);
        }
    };

    let topo_b = build_topology(blocks_b);
    println!("  State B: {} atoms", topo_b.num_atoms());

    // Check compatibility
    if topo_a.num_atoms() != topo_b.num_atoms() {
        eprintln!(
            "Error: Topologies must have the same number of atoms ({} vs {})",
            topo_a.num_atoms(),
            topo_b.num_atoms()
        );
        process::exit(1);
    }

    // Count perturbed atoms
    let mut n_perturbed = 0;
    for i in 0..topo_a.num_atoms() {
        let is_perturbed = topo_a.iac[i] != topo_b.iac[i]
            || (topo_a.mass[i] - topo_b.mass[i]).abs() > 1e-6
            || (topo_a.charge[i] - topo_b.charge[i]).abs() > 1e-6;
        if is_perturbed {
            n_perturbed += 1;
        }
    }

    println!();
    println!("Perturbation analysis:");
    println!("  Total atoms:      {}", topo_a.num_atoms());
    println!("  Perturbed atoms:  {}", n_perturbed);
    println!();

    if n_perturbed == 0 {
        eprintln!("Warning: No perturbed atoms found! States A and B appear identical.");
    }

    // Write perturbation topology
    println!("Writing perturbation topology to {}...", pt_args.output);
    let mut writer = match PtpWriter::new(&pt_args.output) {
        Ok(w) => w,
        Err(e) => {
            eprintln!("Error creating output file: {}", e);
            process::exit(1);
        }
    };

    let title = format!(
        "make_pt_top generated perturbation topology\nState A: {}\nState B: {}",
        pt_args.topo_a, pt_args.topo_b
    );

    if let Err(e) = writer.write(
        &topo_a,
        &topo_b,
        pt_args.alpha_lj,
        pt_args.alpha_crf,
        &title,
    ) {
        eprintln!("Error writing perturbation topology: {}", e);
        process::exit(1);
    }

    println!("Done! Perturbation topology written to {}", pt_args.output);
    println!();
    println!("Next steps:");
    println!("  1. Use this .ptp file with the md binary: md @topo ... @ptp {}", pt_args.output);
    println!("  2. Run FEP simulations at different λ values (0.0 to 1.0)");
    println!("  3. Analyze with GROMOS++ bar or ext_ti_ana tools");
}
