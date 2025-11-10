//! Validation tests against C++ GROMOS
//!
//! These tests load GROMOS input files and perform energy calculations
//! to verify that the Rust implementation produces the same results as C++ GROMOS.

use gromos_rs::io::coordinate::read_coordinate_file;
use gromos_rs::io::topology::{read_topology_file, build_topology};
use gromos_rs::math::{Periodicity, Rectangular, Vacuum};
use gromos_rs::configuration::BoxType;

#[test]
fn test_cg16_single_point_energy() {
    // Load topology
    let parsed_topo = read_topology_file("../md++/src/check/data/cg16.topo");

    if parsed_topo.is_err() {
        println!("Skipping test: cg16.topo not found");
        return;
    }

    let parsed_topo = parsed_topo.unwrap();
    let topo = build_topology(parsed_topo);

    // Load coordinates
    let mut conf = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

    if conf.is_err() {
        println!("Skipping test: cg16.conf not found");
        return;
    }

    let mut conf = conf.unwrap();

    println!("\n========== CG16 System ==========");
    println!("Number of atoms: {}", topo.num_atoms());
    println!("Number of bonds: {}", topo.solute.bonds.len());
    println!("Number of angles: {}", topo.solute.angles.len());

    // Print atom info
    for i in 0..topo.num_atoms().min(4) {
        println!("Atom {}: IAC={}, mass={}, charge={}, pos=({:.3}, {:.3}, {:.3})",
                 i,
                 topo.iac[i],
                 topo.mass[i],
                 topo.charge[i],
                 conf.current().pos[i].x,
                 conf.current().pos[i].y,
                 conf.current().pos[i].z);
    }

    // Calculate nonbonded energy (LJ only, no electrostatics for CG)
    let box_config = conf.current().box_config.clone();
    let periodicity = match box_config.box_type {
        BoxType::Vacuum => Periodicity::Vacuum(Vacuum),
        BoxType::Rectangular => {
            let dims = box_config.dimensions();
            Periodicity::Rectangular(Rectangular::new(dims))
        }
        _ => Periodicity::Vacuum(Vacuum),  // Simplified for now
    };

    let cutoff = 1.4;  // Standard GROMOS cutoff
    let rf_epsilon = 1.0;  // Vacuum
    let rf_kappa = 0.0;

    // Calculate LJ energy for all pairs
    let mut total_lj = 0.0;

    for i in 0..topo.num_atoms() {
        for j in (i + 1)..topo.num_atoms() {
            // Check if excluded
            if topo.is_excluded(i, j) {
                continue;
            }

            let r = periodicity.nearest_image(
                conf.current().pos[i],
                conf.current().pos[j]
            );

            let dist2 = r.length_squared() as f64;
            let dist = dist2.sqrt();

            if dist > cutoff {
                continue;
            }

            // Get LJ parameters
            let iac_i = topo.iac[i];
            let iac_j = topo.iac[j];

            if iac_i == 0 || iac_j == 0 || iac_i > topo.lj_parameters.len() || iac_j > topo.lj_parameters[0].len() {
                continue;
            }

            let lj_params = &topo.lj_parameters[iac_i - 1][iac_j - 1];

            // V_LJ = C12/r^12 - C6/r^6
            let r6 = dist2 * dist2 * dist2;
            let r12 = r6 * r6;

            let lj_energy = lj_params.c12 / r12 - lj_params.c6 / r6;

            total_lj += lj_energy;

            println!("Pair {}-{}: dist={:.4}, C12={:.6e}, C6={:.6e}, E_LJ={:.6}",
                     i, j, dist, lj_params.c12, lj_params.c6, lj_energy);
        }
    }

    println!("\n========== Energies ==========");
    println!("Total LJ Energy: {:.6} kJ/mol", total_lj);
    println!("\nNOTE: Compare this with C++ GROMOS output for the same system");
    println!("      Run: md++ --topo cg16.topo --conf cg16.conf --verb 3");
}

#[test]
fn test_load_both_files() {
    // Simple test to verify we can load both files
    let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");
    let conf_result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

    if topo_result.is_err() || conf_result.is_err() {
        println!("Skipping: files not available");
        return;
    }

    let parsed_topo = topo_result.unwrap();
    let topo = build_topology(parsed_topo);
    let conf = conf_result.unwrap();

    // Verify consistency
    assert_eq!(topo.num_atoms(), conf.current().pos.len(),
               "Topology and coordinates have different atom counts");

    println!("✓ Loaded {} atoms from topology and coordinates", topo.num_atoms());
    println!("✓ Coordinates and topology are consistent");
}
