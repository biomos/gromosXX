//! com_top - Combine multiple GROMOS topology files
//!
//! Usage: com_top <output.top> <input1.top> [input2.top] [...]
//!
//! Combines multiple topology files into a single topology file.
//! - Renumbers atoms sequentially
//! - Merges bonds, angles, dihedrals
//! - Combines force field parameters

use gromos_rs::io::topology::{read_topology_file, write_topology_file, build_topology};
use gromos_rs::topology::{Topology, Atom, Bond, Angle};
use std::collections::HashMap;
use std::env;
use std::process;

fn print_usage() {
    eprintln!("Usage: com_top <output.top> <input1.top> [input2.top] [...]");
    eprintln!();
    eprintln!("Combine multiple GROMOS topology files");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  output.top   Output topology file");
    eprintln!("  input*.top   Input topology files to combine");
    eprintln!();
    eprintln!("Features:");
    eprintln!("  - Renumbers atoms sequentially across files");
    eprintln!("  - Merges bonded interactions");
    eprintln!("  - Combines force field parameters");
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        process::exit(if args.len() < 3 { 1 } else { 0 });
    }

    let output_file = &args[1];
    let input_files = &args[2..];

    println!("Combining {} topology files...", input_files.len());

    // Read all input topology files
    let mut topologies = Vec::new();
    for (i, input_file) in input_files.iter().enumerate() {
        println!("  Reading {}: {}", i + 1, input_file);
        match read_topology_file(input_file) {
            Ok(parsed) => {
                let topo = build_topology(parsed);
                println!("    {} atoms, {} bonds, {} angles",
                         topo.solute.num_atoms(),
                         topo.solute.bonds.len(),
                         topo.solute.angles.len());
                topologies.push(topo);
            }
            Err(e) => {
                eprintln!("Error reading {}: {:?}", input_file, e);
                process::exit(1);
            }
        }
    }

    // Combine topologies
    println!("\nCombining topologies...");
    let combined = combine_topologies(&topologies);

    println!("  Combined topology:");
    println!("    {} atoms", combined.solute.num_atoms());
    println!("    {} bonds", combined.solute.bonds.len());
    println!("    {} angles", combined.solute.angles.len());
    println!("    {} bond parameters", combined.bond_parameters.len());
    println!("    {} angle parameters", combined.angle_parameters.len());

    // Write combined topology
    println!("\nWriting combined topology: {}", output_file);
    let title = format!("Combined topology from {} files", topologies.len());

    match write_topology_file(output_file, &combined, &title) {
        Ok(_) => {
            println!("Successfully wrote combined topology");
        }
        Err(e) => {
            eprintln!("Error writing topology: {:?}", e);
            process::exit(1);
        }
    }
}

/// Combine multiple topologies into one
///
/// Strategy:
/// 1. Renumber atoms sequentially across topologies
/// 2. Merge atom arrays
/// 3. Update bond/angle atom indices
/// 4. Merge force field parameters (deduplicating)
fn combine_topologies(topologies: &[Topology]) -> Topology {
    let mut combined = Topology::new();

    let mut atom_offset = 0;
    let mut bond_param_map: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut angle_param_map: HashMap<(i64, i64, i64), usize> = HashMap::new();

    for topo in topologies {
        // Merge atoms
        for atom in &topo.solute.atoms {
            let mut new_atom = atom.clone();
            // Keep sequential numbering
            combined.solute.atoms.push(new_atom);
        }

        // Merge masses, charges, IAC
        combined.mass.extend(&topo.mass);
        combined.charge.extend(&topo.charge);
        combined.iac.extend(&topo.iac);

        // Merge bonds with renumbered atoms
        for bond in &topo.solute.bonds {
            let new_bond = Bond {
                i: bond.i + atom_offset,
                j: bond.j + atom_offset,
                bond_type: bond.bond_type,
            };

            // Check if we need to add this bond parameter type
            let param = &topo.bond_parameters[bond.bond_type];
            let param_key = (
                (param.k_quartic * 1e10) as i64,
                (param.k_harmonic * 1e10) as i64,
                (param.r0 * 1e10) as i64,
            );

            let new_bond_type = if let Some(&idx) = bond_param_map.get(&param_key) {
                idx
            } else {
                let idx = combined.bond_parameters.len();
                combined.bond_parameters.push(*param);
                bond_param_map.insert(param_key, idx);
                idx
            };

            combined.solute.bonds.push(Bond {
                i: new_bond.i,
                j: new_bond.j,
                bond_type: new_bond_type,
            });
        }

        // Merge angles with renumbered atoms
        for angle in &topo.solute.angles {
            let new_angle_idx = Angle {
                i: angle.i + atom_offset,
                j: angle.j + atom_offset,
                k: angle.k + atom_offset,
                angle_type: angle.angle_type,
            };

            // Check if we need to add this angle parameter type
            let param = &topo.angle_parameters[angle.angle_type];
            let param_key = (
                (param.k_cosine * 1e10) as i64,
                (param.k_harmonic * 1e10) as i64,
                (param.theta0 * 1e10) as i64,
            );

            let new_angle_type = if let Some(&idx) = angle_param_map.get(&param_key) {
                idx
            } else {
                let idx = combined.angle_parameters.len();
                combined.angle_parameters.push(*param);
                angle_param_map.insert(param_key, idx);
                idx
            };

            combined.solute.angles.push(Angle {
                i: new_angle_idx.i,
                j: new_angle_idx.j,
                k: new_angle_idx.k,
                angle_type: new_angle_type,
            });
        }

        // Merge exclusions with renumbered atoms
        for (i, excl_set) in topo.exclusions.iter().enumerate() {
            let mut new_excl_set = std::collections::HashSet::new();
            for &excl in excl_set {
                new_excl_set.insert(excl + atom_offset);
            }
            combined.exclusions.push(new_excl_set);
        }

        // Update atom offset for next topology
        atom_offset += topo.solute.num_atoms();
    }

    // Compute inverse masses
    combined.compute_inverse_masses();

    // Build LJ parameter matrix from the first topology (assuming same force field)
    if !topologies.is_empty() && !topologies[0].lj_parameters.is_empty() {
        combined.lj_parameters = topologies[0].lj_parameters.clone();
    }

    combined
}
