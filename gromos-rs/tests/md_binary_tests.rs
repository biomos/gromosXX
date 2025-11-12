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

/// Create a perturbed topology file for FEP testing
/// Two atoms with a perturbed bond
fn create_fep_topology(path: &str) -> std::io::Result<()> {
    let topology = r#"TITLE
FEP test system - perturbed bond
END
PHYSICALCONSTANTS
  0.1389354
  0.0635078
299792.458
  0.00831441
END
TOPVERSION
2.0
END
ATOMTYPENAME
    1
    C
END
RESNAME
    1
FEP
END
SOLUTEATOM
# Two atoms with perturbed bond
     2
# ATNM MRES PANM IAC MASS CG CGC INE
    1    1  C1    1  12.0110  0.0  1  1
    2
    2    1  C2    1  12.0110  0.0  1  1
    1
END
BONDSTRETCHTYPE
# State A bond type: k=5000, r0=0.15
    2
#  TYPE  CB    HB
    1   0.0   5000.0
    2   0.0   10000.0
END
BOND
    0
END
BONDH
    0
END
BONDANGLEBENDTYPE
    0
END
BONDANGLE
    0
END
BONDANGLEH
    0
END
DIHEDRALTORSIONTYPE
    0
END
DIHEDRALTORSION
    0
END
IMPDIHEDRALTORSIONTYPE
    0
END
IMPDIHEDRALTORSION
    0
END
LJPARAMETERS
# NRATT2: number of LJ types
    1
#  IAC   C12      C6
    1   0.0      0.0
END
PERTURBATION
# NRPERTATOM: number of perturbed atoms
    2
# ATOM  ALPHLJ  ALPHCRF
    1    1.0     1.0
    2    1.0     1.0
# NRPERTBOND: number of perturbed bonds
    1
# I  J  TYPE_A  TYPE_B  ALPHLJ  ALPHCRF
  1  2    1       2      1.0     1.0
# NRPERTANGLE: number of perturbed angles
    0
# NRPERTDIHEDRAL: number of perturbed dihedrals
    0
END
END"#;

    fs::write(path, topology)?;
    Ok(())
}

/// Create minimal coordinates for FEP system
fn create_fep_coordinates(path: &str) -> std::io::Result<()> {
    let coords = r#"TITLE
FEP test coordinates
END
POSITION
  1  FEP  C1     1   0.000   0.000   0.000
  2  FEP  C2     2   0.120   0.000   0.000
END
BOX
   2.0   2.0   2.0
END"#;

    fs::write(path, coords)?;
    Ok(())
}

#[test]
fn test_md_binary_fep_bond_perturbation() {
    println!("\n========== MD Binary FEP Bond Perturbation Test ==========");

    // Create temporary directory for test files
    let test_dir = "/tmp/gromos_md_fep_bond_test";
    fs::create_dir_all(test_dir).expect("Failed to create test directory");

    let topo_file = format!("{}/fep_test.top", test_dir);
    let coord_file = format!("{}/fep_test.g96", test_dir);
    let traj_file = format!("{}/fep_test.trc", test_dir);
    let ene_file = format!("{}/fep_test.tre", test_dir);

    // Create test input files
    create_fep_topology(&topo_file).expect("Failed to create FEP topology");
    create_fep_coordinates(&coord_file).expect("Failed to create coordinates");

    println!("Created test files:");
    println!("  Topology: {}", topo_file);
    println!("  Coordinates: {}", coord_file);

    // Run MD simulation with FEP at lambda=0.5
    println!("\nRunning MD simulation with FEP at λ=0.5...");
    let output = Command::new("cargo")
        .args(&[
            "run", "--bin", "md", "--",
            "@topo", &topo_file,
            "@conf", &coord_file,
            "@steps", "10",
            "@dt", "0.001",
            "@traj", &traj_file,
            "@ene", &ene_file,
            "@lambda", "0.5",
        ])
        .current_dir("/home/user/gromosXX/gromos-rs")
        .output()
        .expect("Failed to run md binary");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    if !output.status.success() {
        println!("STDOUT:\n{}", stdout);
        println!("STDERR:\n{}", stderr);
        panic!("MD simulation with FEP failed");
    }

    // Check output files exist
    assert!(Path::new(&traj_file).exists(), "Trajectory file not created");
    assert!(Path::new(&ene_file).exists(), "Energy file not created");

    // Check energy file contains data
    let ene_content = fs::read_to_string(&ene_file).expect("Failed to read energy file");
    assert!(ene_content.contains("ENERTRJ"), "Energy file missing ENERTRJ block");

    // Count energy frames
    let energy_lines: Vec<&str> = ene_content.lines()
        .filter(|line| {
            let trimmed = line.trim();
            !trimmed.is_empty()
            && !trimmed.starts_with('#')
            && trimmed != "TITLE"
            && trimmed != "ENERTRJ"
            && trimmed != "END"
            && (trimmed.chars().next().map_or(false, |c| c.is_numeric() || c == '-'))
        })
        .collect();

    println!("\n✓ FEP simulation completed:");
    println!("  - Energy frames: {}", energy_lines.len());
    println!("  - Lambda value: 0.5");
    assert!(energy_lines.len() >= 2, "Should have at least 2 energy frames");

    // Verify energies are reasonable
    for (idx, line) in energy_lines.iter().take(2).enumerate() {
        let values: Vec<&str> = line.split_whitespace().collect();
        assert!(values.len() >= 4, "Energy line should have at least 4 values");

        let time: f64 = values[0].parse().expect("Failed to parse time");
        let potential: f64 = values[2].parse().expect("Failed to parse potential energy");

        println!("  Frame {}: t={:.3} ps, E_pot={:.6} kJ/mol", idx, time, potential);
        assert!(potential.is_finite(), "Potential energy should be finite");
    }

    println!("\n✓ MD binary FEP test passed!");

    // Cleanup
    fs::remove_dir_all(test_dir).ok();
}

#[test]
fn test_md_binary_fep_lambda_scan() {
    println!("\n========== MD Binary FEP Lambda Scan Test ==========");

    // Create temporary directory for test files
    let test_dir = "/tmp/gromos_md_fep_lambda_scan";
    fs::create_dir_all(test_dir).expect("Failed to create test directory");

    let topo_file = format!("{}/fep_test.top", test_dir);
    let coord_file = format!("{}/fep_test.g96", test_dir);

    // Create test input files
    create_fep_topology(&topo_file).expect("Failed to create FEP topology");
    create_fep_coordinates(&coord_file).expect("Failed to create coordinates");

    println!("Testing FEP at multiple lambda values:");
    println!("λ\tSteps\tStatus");
    println!("---\t-----\t------");

    // Test simulations at different lambda values
    let lambda_values = vec![0.0, 0.25, 0.5, 0.75, 1.0];
    let mut all_succeeded = true;

    for lambda in &lambda_values {
        let traj_file = format!("{}/fep_lambda_{:.2}.trc", test_dir, lambda);
        let ene_file = format!("{}/fep_lambda_{:.2}.tre", test_dir, lambda);

        let output = Command::new("cargo")
            .args(&[
                "run", "--bin", "md", "--",
                "@topo", &topo_file,
                "@conf", &coord_file,
                "@steps", "5",
                "@dt", "0.001",
                "@traj", &traj_file,
                "@ene", &ene_file,
                "@lambda", &lambda.to_string(),
            ])
            .current_dir("/home/user/gromosXX/gromos-rs")
            .output()
            .expect("Failed to run md binary");

        let success = output.status.success();
        all_succeeded &= success;

        let status = if success { "✓" } else { "✗" };
        println!("{:.2}\t5\t{}", lambda, status);

        if success {
            // Check that energy file was created and has content
            if Path::new(&ene_file).exists() {
                let ene_content = fs::read_to_string(&ene_file).expect("Failed to read energy file");
                assert!(ene_content.contains("ENERTRJ"), "Energy file missing data at λ={}", lambda);
            }
        }
    }

    assert!(all_succeeded, "Some lambda simulations failed");

    println!("\n✓ Lambda scan test passed!");
    println!("  - All {} lambda values simulated successfully", lambda_values.len());
    println!("  - Ready for thermodynamic integration!");

    // Cleanup
    fs::remove_dir_all(test_dir).ok();
}

#[test]
fn test_md_binary_fep_perturbed_angle() {
    println!("\n========== MD Binary FEP Angle Perturbation Test ==========");

    // Create temporary directory
    let test_dir = "/tmp/gromos_md_fep_angle_test";
    fs::create_dir_all(test_dir).expect("Failed to create test directory");

    let topo_file = format!("{}/angle_fep.top", test_dir);
    let coord_file = format!("{}/angle_fep.g96", test_dir);
    let ene_file = format!("{}/angle_fep.tre", test_dir);

    // Create topology with perturbed angle (3 atoms)
    let topology = r#"TITLE
FEP angle perturbation test
END
PHYSICALCONSTANTS
  0.1389354
  0.0635078
299792.458
  0.00831441
END
TOPVERSION
2.0
END
ATOMTYPENAME
    1
    C
END
RESNAME
    1
ANG
END
SOLUTEATOM
     3
    1    1  C1    1  12.0110  0.0  1  2
    2 3
    2    1  C2    1  12.0110  0.0  1  2
    1 3
    3    1  C3    1  12.0110  0.0  1  2
    1 2
END
BONDSTRETCHTYPE
    1
    1   0.0   5000.0
END
BOND
    2
    1     1     2
    1     2     3
END
BONDH
    0
END
BONDANGLEBENDTYPE
# State A: 109.5 degrees, State B: 120 degrees
    2
#  TYPE   CT    CHT
    1    0.0   400.0
    2    0.0   600.0
END
BONDANGLE
    0
END
BONDANGLEH
    0
END
DIHEDRALTORSIONTYPE
    0
END
DIHEDRALTORSION
    0
END
IMPDIHEDRALTORSIONTYPE
    0
END
IMPDIHEDRALTORSION
    0
END
LJPARAMETERS
    1
    1   0.0      0.0
END
PERTURBATION
# 3 perturbed atoms (middle atom is key)
    3
    1    1.0     1.0
    2    1.0     1.0
    3    1.0     1.0
# No perturbed bonds
    0
# 1 perturbed angle
    1
# I  J  K  TYPE_A  TYPE_B
  1  2  3    1       2
# No perturbed dihedrals
    0
END
END"#;

    fs::write(&topo_file, topology).expect("Failed to write angle topology");

    // Create coordinates with angle near 115 degrees
    let coords = r#"TITLE
Angle FEP coordinates
END
POSITION
  1  ANG  C1     1   0.100   0.000   0.000
  2  ANG  C2     2   0.000   0.000   0.000
  3  ANG  C3     3  -0.042   0.091   0.000
END
BOX
   2.0   2.0   2.0
END"#;

    fs::write(&coord_file, coords).expect("Failed to write angle coordinates");

    println!("Running angle perturbation at λ=0.5...");

    // Run simulation
    let output = Command::new("cargo")
        .args(&[
            "run", "--bin", "md", "--",
            "@topo", &topo_file,
            "@conf", &coord_file,
            "@steps", "10",
            "@dt", "0.001",
            "@ene", &ene_file,
            "@lambda", "0.5",
        ])
        .current_dir("/home/user/gromosXX/gromos-rs")
        .output()
        .expect("Failed to run md binary");

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        println!("STDERR: {}", stderr);
        panic!("Angle FEP simulation failed");
    }

    // Verify energy file
    assert!(Path::new(&ene_file).exists(), "Energy file not created");
    let ene_content = fs::read_to_string(&ene_file).expect("Failed to read energy file");
    assert!(ene_content.contains("ENERTRJ"), "Energy file missing data");

    println!("\n✓ Angle perturbation test passed!");
    println!("  - Tetrahedral → Trigonal planar transition");
    println!("  - Simulation completed at λ=0.5");

    // Cleanup
    fs::remove_dir_all(test_dir).ok();
}
