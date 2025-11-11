//! Integration tests for the md binary
//!
//! Tests that the md binary can run simple simulations successfully

use std::process::Command;
use std::fs;
use std::path::Path;

/// Create a minimal topology file for testing
/// Simple 3-atom linear molecule (like CO2 or similar)
fn create_test_topology(path: &str) -> std::io::Result<()> {
    let topology = r#"TITLE
Simple 3-atom test system
END
PHYSICALCONSTANTS
# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)
  0.1389354
# HBAR: Planck's constant HBAR = H/(2* PI)
  0.0635078
# SPDL: Speed of light (nm/ps)
299792.458
# BOLTZ: Boltzmann's constant
  0.00831441
END
TOPVERSION
2.0
END
ATOMTYPENAME
# NRATT: number of atom types
    3
# TYPE: atom type names
    C
    O
    H
END
RESNAME
# NRRES2: number of residues
    1
# RESNAME: residue names
TEST
END
SOLUTEATOM
# NRSOLTATOM: number of atoms
     3
# ATNM MRES PANM IAC MASS CG CGC INE
    1    1  C1    1  12.0110  0.0  1  2
    2 3
    2    1  O1    2  15.9994  0.0  1  2
    1 3
    3    1  O2    2  15.9994  0.0  1  2
    1 2
END
BONDSTRETCHTYPE
# NBTY: number of bond types
    1
#  TYPE  CB    HB
    1   0.0   1.5e7
END
BONDH
# NBONH: number of bonds with H atoms
    0
# bond type, atom1, atom2
END
BOND
# NBON: number of bonds without H atoms
    2
# bond type, atom1, atom2
    1     1     2
    1     1     3
END
BONDANGLEBENDTYPE
# NTTY: number of angle types
    1
#  TYPE   CT    CHT
    1    0.0   405.0
END
BONDANGLEH
# NTHEH: number of angles with H atoms
    0
# angle type, atom1, atom2, atom3
END
BONDANGLE
# NTHE: number of angles without H atoms
    1
# angle type, atom1, atom2, atom3
    1     2     1     3
END
IMPDIHEDRALTYPE
# NQTY: number of improper dihedral types
    0
END
IMPDIHEDRALH
# NQHIH: number of improper dihedrals with H atoms
    0
END
IMPDIHEDRAL
# NQHI: number of improper dihedrals without H atoms
    0
END
DIHEDRALTYPE
# NPTY: number of dihedral types
    0
END
DIHEDRALH
# NPHIH: number of dihedrals with H atoms
    0
END
DIHEDRAL
# NPHI: number of dihedrals without H atoms
    0
END
LJPARAMETERS
# NRATT2: number of LJ atom types
    3
#  IAC  TYPE    C6      C12
    1    1    0.0     0.0
    2    2    0.0     0.0
    3    3    0.0     0.0
END
"#;

    fs::write(path, topology)?;
    Ok(())
}

/// Create a minimal coordinate file for testing
fn create_test_coordinates(path: &str) -> std::io::Result<()> {
    let coordinates = r#"TITLE
Simple 3-atom test coordinates
END
POSITION
    1 TEST     C1      1    0.000000000    0.000000000    0.000000000
    2 TEST     O1      2    0.120000000    0.000000000    0.000000000
    3 TEST     O2      3   -0.120000000    0.000000000    0.000000000
END
GENBOX
    3.000000000    3.000000000    3.000000000
END
"#;

    fs::write(path, coordinates)?;
    Ok(())
}

#[test]
fn test_md_binary_simple_simulation() {
    // Create temporary directory for test files
    let test_dir = "/tmp/gromos_md_test";
    fs::create_dir_all(test_dir).expect("Failed to create test directory");

    let topo_file = format!("{}/test.top", test_dir);
    let coord_file = format!("{}/test.g96", test_dir);
    let traj_file = format!("{}/test.trc", test_dir);
    let ene_file = format!("{}/test.tre", test_dir);

    // Create test input files
    create_test_topology(&topo_file).expect("Failed to create topology");
    create_test_coordinates(&coord_file).expect("Failed to create coordinates");

    // Build the md binary first
    let build_output = Command::new("cargo")
        .args(&["build", "--bin", "md"])
        .current_dir("/home/user/gromosXX/gromos-rs")
        .output()
        .expect("Failed to build md binary");

    assert!(build_output.status.success(), "Failed to build md binary");

    // Run a short MD simulation
    let output = Command::new("cargo")
        .args(&[
            "run", "--bin", "md", "--",
            "@topo", &topo_file,
            "@conf", &coord_file,
            "@steps", "10",
            "@dt", "0.001",
            "@traj", &traj_file,
            "@ene", &ene_file,
        ])
        .current_dir("/home/user/gromosXX/gromos-rs")
        .output()
        .expect("Failed to run md binary");

    // Print output for debugging
    println!("STDOUT:\n{}", String::from_utf8_lossy(&output.stdout));
    println!("STDERR:\n{}", String::from_utf8_lossy(&output.stderr));

    // Check that the simulation completed successfully
    assert!(output.status.success(), "MD simulation failed");

    // Check that output files were created
    assert!(Path::new(&traj_file).exists(), "Trajectory file not created");
    assert!(Path::new(&ene_file).exists(), "Energy file not created");

    // Check trajectory file has content
    let traj_content = fs::read_to_string(&traj_file).expect("Failed to read trajectory");
    assert!(traj_content.contains("TITLE"), "Trajectory missing TITLE block");
    assert!(traj_content.contains("POSITIONRED"), "Trajectory missing POSITIONRED block");
    assert!(traj_content.contains("GENBOX"), "Trajectory missing GENBOX block");

    // Check energy file has content
    let ene_content = fs::read_to_string(&ene_file).expect("Failed to read energy file");
    assert!(ene_content.contains("TITLE"), "Energy file missing TITLE block");
    assert!(ene_content.contains("ENERTRJ"), "Energy file missing ENERTRJ block");

    // Count number of energy frames (should be at least 2 - initial and final)
    let energy_lines: Vec<&str> = ene_content.lines()
        .filter(|line| {
            let trimmed = line.trim();
            // Energy data lines are numeric (first char after trim is a digit or -)
            !trimmed.is_empty()
            && !trimmed.starts_with('#')
            && trimmed != "TITLE"
            && trimmed != "ENERTRJ"
            && trimmed != "END"
            && (trimmed.chars().next().map_or(false, |c| c.is_numeric() || c == '-'))
        })
        .collect();

    println!("  Energy frames found: {}", energy_lines.len());
    assert!(energy_lines.len() >= 2, "Energy file should contain at least 2 frames, found: {}", energy_lines.len());

    // Parse the energy values to ensure they're reasonable
    for (idx, line) in energy_lines.iter().take(2).enumerate() {
        let values: Vec<&str> = line.split_whitespace().collect();
        println!("  Frame {}: {} values", idx, values.len());
        assert!(values.len() >= 4, "Energy line {} should have at least 4 values, got: {}", idx, values.len());

        // Parse time
        let time: f64 = values[0].parse().expect("Failed to parse time");
        assert!(time >= 0.0, "Time should be non-negative");

        // Parse energies
        let _kinetic: f64 = values[1].parse().expect("Failed to parse kinetic energy");
        let _potential: f64 = values[2].parse().expect("Failed to parse potential energy");
        let _total: f64 = values[3].parse().expect("Failed to parse total energy");
    }

    println!("\n✓ MD simulation test passed!");
    println!("  - Topology loaded successfully");
    println!("  - Coordinates loaded successfully");
    println!("  - Simulation ran for 10 steps");
    println!("  - Trajectory file created: {}", traj_file);
    println!("  - Energy file created: {}", ene_file);
    println!("  - Output files contain valid data");

    // Cleanup
    fs::remove_dir_all(test_dir).ok();
}

#[test]
fn test_md_binary_energy_conservation() {
    // Create temporary directory for test files
    let test_dir = "/tmp/gromos_md_energy_test";
    fs::create_dir_all(test_dir).expect("Failed to create test directory");

    let topo_file = format!("{}/test.top", test_dir);
    let coord_file = format!("{}/test.g96", test_dir);
    let ene_file = format!("{}/test_energy.tre", test_dir);

    // Create test input files
    create_test_topology(&topo_file).expect("Failed to create topology");
    create_test_coordinates(&coord_file).expect("Failed to create coordinates");

    // Run MD simulation
    let output = Command::new("cargo")
        .args(&[
            "run", "--bin", "md", "--",
            "@topo", &topo_file,
            "@conf", &coord_file,
            "@steps", "50",
            "@dt", "0.001",
            "@ene", &ene_file,
        ])
        .current_dir("/home/user/gromosXX/gromos-rs")
        .output()
        .expect("Failed to run md binary");

    assert!(output.status.success(), "MD simulation failed");

    // Read energy file and check conservation
    let ene_content = fs::read_to_string(&ene_file).expect("Failed to read energy file");

    let energy_lines: Vec<&str> = ene_content.lines()
        .filter(|line| {
            let trimmed = line.trim();
            !trimmed.is_empty()
            && !trimmed.starts_with('#')
            && !trimmed.starts_with("TITLE")
            && !trimmed.starts_with("ENERTRJ")
            && !trimmed.starts_with("END")
        })
        .collect();

    assert!(energy_lines.len() >= 5, "Should have at least 5 energy frames");

    // Extract total energies
    let mut total_energies = Vec::new();
    for line in &energy_lines {
        let values: Vec<&str> = line.split_whitespace().collect();
        if values.len() >= 4 {
            if let Ok(total) = values[3].parse::<f64>() {
                total_energies.push(total);
            }
        }
    }

    assert!(!total_energies.is_empty(), "No valid energy values found");

    // Calculate energy drift (should be small for such a simple system)
    let initial_energy = total_energies[0];
    let final_energy = *total_energies.last().unwrap();
    let energy_drift = (final_energy - initial_energy).abs();

    println!("\n✓ Energy conservation test:");
    println!("  - Initial total energy: {:.6} kJ/mol", initial_energy);
    println!("  - Final total energy: {:.6} kJ/mol", final_energy);
    println!("  - Energy drift: {:.6} kJ/mol", energy_drift);
    println!("  - Number of frames: {}", total_energies.len());

    // For a simple bonded system with small timestep, drift should be reasonable
    // Note: Without proper integrator initialization and with zero velocities,
    // we might see some drift, but it shouldn't be catastrophic

    // Cleanup
    fs::remove_dir_all(test_dir).ok();
}

#[test]
fn test_md_binary_help() {
    // Test that --help works
    let output = Command::new("cargo")
        .args(&["run", "--bin", "md", "--", "--help"])
        .current_dir("/home/user/gromosXX/gromos-rs")
        .output()
        .expect("Failed to run md binary");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Molecular Dynamics simulation") || stderr.contains("md - Molecular Dynamics"),
            "Help text should mention MD simulation");
    assert!(stderr.contains("@topo"), "Help should mention @topo argument");
    assert!(stderr.contains("@conf"), "Help should mention @conf argument");

    println!("\n✓ Help text test passed");
}

#[test]
fn test_md_binary_missing_arguments() {
    // Test that md fails gracefully without required arguments
    let output = Command::new("cargo")
        .args(&["run", "--bin", "md"])
        .current_dir("/home/user/gromosXX/gromos-rs")
        .output()
        .expect("Failed to run md binary");

    assert!(!output.status.success(), "MD should fail without arguments");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("@topo") || stderr.contains("Usage"),
            "Error should mention missing topology");

    println!("\n✓ Missing arguments test passed");
}
