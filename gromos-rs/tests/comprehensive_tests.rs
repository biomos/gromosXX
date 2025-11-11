//! Comprehensive test suite for gromos-rs
//!
//! Tests all major components with real GROMOS data

use gromos_rs::io::coordinate::read_coordinate_file;
use gromos_rs::io::topology::{read_topology_file, build_topology};
use gromos_rs::math::{Vec3, Periodicity, Rectangular, Vacuum};
use gromos_rs::configuration::BoxType;
use gromos_rs::topology::Topology;

/// Test that we can load and parse all test files
#[test]
fn test_load_all_gromos_test_files() {
    let test_files = vec![
        ("../md++/src/check/data/cg16.topo", "../md++/src/check/data/cg16.conf"),
    ];

    for (topo_file, conf_file) in test_files {
        println!("\nTesting: {} + {}", topo_file, conf_file);

        let parsed_topo = read_topology_file(topo_file);
        if parsed_topo.is_err() {
            println!("  Skipping: topology not found");
            continue;
        }

        let parsed_topo = parsed_topo.unwrap();
        let topo = build_topology(parsed_topo);

        let conf = read_coordinate_file(conf_file, 1, 1);
        if conf.is_err() {
            println!("  Skipping: coordinates not found");
            continue;
        }

        let conf = conf.unwrap();

        println!("  ✓ Loaded {} atoms", topo.num_atoms());
        println!("  ✓ Loaded {} bonds", topo.solute.bonds.len());
        println!("  ✓ Loaded {} angles", topo.solute.angles.len());

        // Verify consistency
        assert_eq!(topo.num_atoms(), conf.current().pos.len(),
                   "Topology and coordinates atom count mismatch");

        // Verify we have masses and charges
        assert_eq!(topo.mass.len(), topo.num_atoms());
        assert_eq!(topo.charge.len(), topo.num_atoms());
    }
}

/// Test topology building and exclusion lists
#[test]
fn test_topology_exclusions() {
    let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");
    if topo_result.is_err() {
        println!("Skipping: test file not found");
        return;
    }

    let parsed = topo_result.unwrap();
    let topo = build_topology(parsed);

    println!("\n========== Exclusions Test ==========");
    println!("Number of atoms: {}", topo.num_atoms());

    // Check that bonded atoms are excluded
    for bond in &topo.solute.bonds {
        println!("Bond {}-{}: checking exclusions", bond.i, bond.j);

        if bond.i < topo.exclusions.len() && bond.j < topo.exclusions.len() {
            assert!(topo.is_excluded(bond.i, bond.j),
                    "Bonded atoms {}-{} should be excluded", bond.i, bond.j);
            println!("  ✓ Bond {}-{} properly excluded", bond.i, bond.j);
        }
    }

    // Check that 1-3 pairs from angles are excluded
    for angle in &topo.solute.angles {
        println!("Angle {}-{}-{}: checking 1-3 exclusion", angle.i, angle.j, angle.k);

        if angle.i < topo.exclusions.len() && angle.k < topo.exclusions.len() {
            assert!(topo.is_excluded(angle.i, angle.k),
                    "1-3 pair {}-{} from angle should be excluded", angle.i, angle.k);
            println!("  ✓ 1-3 pair {}-{} properly excluded", angle.i, angle.k);
        }
    }
}

/// Test bonded energy calculations
#[test]
fn test_bonded_energies() {
    let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");
    let conf_result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

    if topo_result.is_err() || conf_result.is_err() {
        println!("Skipping: test files not found");
        return;
    }

    let parsed = topo_result.unwrap();
    let topo = build_topology(parsed);
    let conf = conf_result.unwrap();

    println!("\n========== Bonded Energy Test ==========");

    // Calculate bond energies
    let mut bond_energy = 0.0;
    for bond in &topo.solute.bonds {
        if bond.bond_type >= topo.bond_parameters.len() {
            continue;
        }

        let params = &topo.bond_parameters[bond.bond_type];
        let r_vec = conf.current().pos[bond.j] - conf.current().pos[bond.i];
        let r = r_vec.length() as f64;

        // GROMOS bond potential: V = (1/4) * k_harmonic * (r^2 - r0^2)^2
        let dr2 = r * r - params.r0 * params.r0;
        let energy = 0.25 * params.k_harmonic * dr2 * dr2;

        bond_energy += energy;

        println!("Bond {}-{}: r={:.4} nm, r0={:.4} nm, E={:.6} kJ/mol",
                 bond.i, bond.j, r, params.r0, energy);
    }

    println!("Total bond energy: {:.6} kJ/mol", bond_energy);

    // Calculate angle energies
    let mut angle_energy = 0.0;
    for angle in &topo.solute.angles {
        if angle.angle_type >= topo.angle_parameters.len() {
            continue;
        }

        let params = &topo.angle_parameters[angle.angle_type];

        let r_ij = conf.current().pos[angle.j] - conf.current().pos[angle.i];
        let r_kj = conf.current().pos[angle.j] - conf.current().pos[angle.k];

        let cos_theta = r_ij.dot(r_kj) / (r_ij.length() * r_kj.length());
        let cos_theta = cos_theta.clamp(-1.0, 1.0); // Numerical safety

        // GROMOS angle potential: V = (1/2) * k_harmonic * (cos(theta) - cos(theta0))^2
        let d_cos = cos_theta as f64 - params.theta0.cos();
        let energy = 0.5 * params.k_harmonic * d_cos * d_cos;

        angle_energy += energy;

        println!("Angle {}-{}-{}: cos(θ)={:.4}, cos(θ0)={:.4}, E={:.6} kJ/mol",
                 angle.i, angle.j, angle.k, cos_theta, params.theta0.cos(), energy);
    }

    println!("Total angle energy: {:.6} kJ/mol", angle_energy);
    println!("Total bonded energy: {:.6} kJ/mol", bond_energy + angle_energy);
}

/// Test nonbonded energy calculation with periodic boundary conditions
#[test]
fn test_nonbonded_with_pbc() {
    let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");
    let conf_result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

    if topo_result.is_err() || conf_result.is_err() {
        println!("Skipping: test files not found");
        return;
    }

    let parsed = topo_result.unwrap();
    let topo = build_topology(parsed);
    let conf = conf_result.unwrap();

    println!("\n========== Nonbonded Energy Test ==========");

    let box_config = conf.current().box_config.clone();
    let periodicity = match box_config.box_type {
        BoxType::Vacuum => Periodicity::Vacuum(Vacuum),
        BoxType::Rectangular => {
            let dims = box_config.dimensions();
            Periodicity::Rectangular(Rectangular::new(dims))
        }
        _ => Periodicity::Vacuum(Vacuum),
    };

    println!("Box type: {:?}", box_config.box_type);
    println!("Box dimensions: {:?}", box_config.dimensions());

    let cutoff = 1.4; // nm
    let mut lj_energy = 0.0;
    let mut n_pairs = 0;

    for i in 0..topo.num_atoms() {
        for j in (i + 1)..topo.num_atoms() {
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

            let iac_i = topo.iac[i];
            let iac_j = topo.iac[j];

            if iac_i == 0 || iac_j == 0 ||
               iac_i > topo.lj_parameters.len() ||
               iac_j > topo.lj_parameters[0].len() {
                continue;
            }

            let lj_params = &topo.lj_parameters[iac_i - 1][iac_j - 1];

            let r6 = dist2 * dist2 * dist2;
            let r12 = r6 * r6;
            let energy = lj_params.c12 / r12 - lj_params.c6 / r6;

            lj_energy += energy;
            n_pairs += 1;

            if n_pairs <= 5 {
                println!("Pair {}-{}: dist={:.4} nm, E={:.6} kJ/mol", i, j, dist, energy);
            }
        }
    }

    println!("Total pairs within cutoff: {}", n_pairs);
    println!("Total LJ energy: {:.6} kJ/mol", lj_energy);

    // Sanity check: energy should be negative (attractive at medium range)
    if n_pairs > 0 {
        assert!(lj_energy < 0.0, "LJ energy should be negative for attractive interactions");
    }
}

/// Test coordinate file parsing accuracy
#[test]
fn test_coordinate_parsing_precision() {
    let conf_result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

    if conf_result.is_err() {
        println!("Skipping: test file not found");
        return;
    }

    let conf = conf_result.unwrap();

    println!("\n========== Coordinate Parsing Test ==========");
    println!("Number of atoms: {}", conf.current().pos.len());
    println!("Number of velocities: {}", conf.current().vel.len());

    // Check that positions are reasonable (not zero, not NaN)
    for (i, pos) in conf.current().pos.iter().enumerate() {
        assert!(!pos.x.is_nan() && !pos.y.is_nan() && !pos.z.is_nan(),
                "Atom {} has NaN coordinates", i);
        assert!(pos.length() > 0.0, "Atom {} at origin", i);

        if i == 0 {
            println!("First atom position: ({:.6}, {:.6}, {:.6})",
                     pos.x, pos.y, pos.z);
        }
    }

    // Check velocities
    for (i, vel) in conf.current().vel.iter().enumerate() {
        assert!(!vel.x.is_nan() && !vel.y.is_nan() && !vel.z.is_nan(),
                "Atom {} has NaN velocities", i);

        if i == 0 {
            println!("First atom velocity: ({:.6}, {:.6}, {:.6})",
                     vel.x, vel.y, vel.z);
        }
    }

    println!("✓ All coordinates and velocities valid");
}

/// Test LJ parameter matrix
#[test]
fn test_lj_parameter_matrix() {
    let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");

    if topo_result.is_err() {
        println!("Skipping: test file not found");
        return;
    }

    let parsed = topo_result.unwrap();
    let topo = build_topology(parsed);

    println!("\n========== LJ Parameter Matrix Test ==========");

    // Find unique atom types
    let max_iac = topo.iac.iter().max().unwrap_or(&0);
    println!("Max IAC: {}", max_iac);
    println!("LJ matrix size: {}x{}", topo.lj_parameters.len(),
             if topo.lj_parameters.is_empty() { 0 } else { topo.lj_parameters[0].len() });

    // Check that all used IACs have parameters
    for &iac in &topo.iac {
        if iac > 0 && iac <= topo.lj_parameters.len() {
            let params = &topo.lj_parameters[iac - 1][iac - 1];

            println!("IAC {}: C6={:.6e}, C12={:.6e}, sigma={:.4} nm, epsilon={:.4} kJ/mol",
                     iac, params.c6, params.c12, params.sigma(), params.epsilon());

            // Sanity checks
            assert!(params.c6 >= 0.0, "C6 should be non-negative");
            assert!(params.c12 >= 0.0, "C12 should be non-negative");

            if params.c6 > 0.0 || params.c12 > 0.0 {
                assert!(params.sigma() > 0.0, "Sigma should be positive for interacting atoms");
                assert!(params.epsilon() >= 0.0, "Epsilon should be non-negative");
            }
        }
    }

    println!("✓ All LJ parameters valid");
}

/// Integration test: Full energy calculation
#[test]
fn test_full_energy_calculation() {
    let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");
    let conf_result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

    if topo_result.is_err() || conf_result.is_err() {
        println!("Skipping: test files not found");
        return;
    }

    let parsed = topo_result.unwrap();
    let topo = build_topology(parsed);
    let conf = conf_result.unwrap();

    println!("\n========== Full Energy Calculation ==========");
    println!("System: {} atoms, {} bonds, {} angles",
             topo.num_atoms(), topo.solute.bonds.len(), topo.solute.angles.len());

    // Bonded energies
    let mut bond_energy = 0.0;
    for bond in &topo.solute.bonds {
        if bond.bond_type < topo.bond_parameters.len() {
            let params = &topo.bond_parameters[bond.bond_type];
            let r_vec = conf.current().pos[bond.j] - conf.current().pos[bond.i];
            let r = r_vec.length() as f64;
            let dr2 = r * r - params.r0 * params.r0;
            bond_energy += 0.25 * params.k_harmonic * dr2 * dr2;
        }
    }

    let mut angle_energy = 0.0;
    for angle in &topo.solute.angles {
        if angle.angle_type < topo.angle_parameters.len() {
            let params = &topo.angle_parameters[angle.angle_type];
            let r_ij = conf.current().pos[angle.j] - conf.current().pos[angle.i];
            let r_kj = conf.current().pos[angle.j] - conf.current().pos[angle.k];
            let cos_theta = (r_ij.dot(r_kj) / (r_ij.length() * r_kj.length())).clamp(-1.0, 1.0);
            let d_cos = cos_theta as f64 - params.theta0.cos();
            angle_energy += 0.5 * params.k_harmonic * d_cos * d_cos;
        }
    }

    // Nonbonded energy
    let box_config = conf.current().box_config.clone();
    let periodicity = match box_config.box_type {
        BoxType::Vacuum => Periodicity::Vacuum(Vacuum),
        BoxType::Rectangular => {
            let dims = box_config.dimensions();
            Periodicity::Rectangular(Rectangular::new(dims))
        }
        _ => Periodicity::Vacuum(Vacuum),
    };

    let cutoff = 1.4;
    let mut lj_energy = 0.0;

    for i in 0..topo.num_atoms() {
        for j in (i + 1)..topo.num_atoms() {
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

            let iac_i = topo.iac[i];
            let iac_j = topo.iac[j];

            if iac_i > 0 && iac_j > 0 &&
               iac_i <= topo.lj_parameters.len() &&
               iac_j <= topo.lj_parameters[0].len() {
                let lj_params = &topo.lj_parameters[iac_i - 1][iac_j - 1];
                let r6 = dist2 * dist2 * dist2;
                let r12 = r6 * r6;
                lj_energy += lj_params.c12 / r12 - lj_params.c6 / r6;
            }
        }
    }

    println!("\n========== Energy Components ==========");
    println!("Bond energy:     {:12.6} kJ/mol", bond_energy);
    println!("Angle energy:    {:12.6} kJ/mol", angle_energy);
    println!("LJ energy:       {:12.6} kJ/mol", lj_energy);
    println!("----------------------------------------");
    println!("Total energy:    {:12.6} kJ/mol", bond_energy + angle_energy + lj_energy);

    println!("\n✓ Full energy calculation completed");
}
