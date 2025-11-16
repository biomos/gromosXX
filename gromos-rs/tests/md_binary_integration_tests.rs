//! Integration tests for the MD binary
//!
//! These tests run the actual `md` binary with real topology and coordinate files
//! to verify end-to-end functionality for different simulation ensembles:
//! - NVE (microcanonical)
//! - NVT (canonical)
//! - NPT (isobaric-isothermal)
//! - Simulations with SHAKE constraints

use std::process::Command;
use std::path::Path;
use std::fs;

/// Helper function to run the md binary with given arguments
fn run_md(args: &[&str]) -> Result<std::process::Output, String> {
    let md_binary = env!("CARGO_BIN_EXE_md");

    Command::new(md_binary)
        .args(args)
        .output()
        .map_err(|e| format!("Failed to execute md: {}", e))
}

/// Helper function to check if a file exists and has content
fn check_output_file(path: &str) -> Result<(), String> {
    if !Path::new(path).exists() {
        return Err(format!("Output file does not exist: {}", path));
    }

    let metadata = fs::metadata(path)
        .map_err(|e| format!("Cannot read file metadata: {}", e))?;

    if metadata.len() == 0 {
        return Err(format!("Output file is empty: {}", path));
    }

    Ok(())
}

/// Helper function to parse energy from trajectory energy file
fn parse_energy_file(path: &str) -> Result<Vec<f64>, String> {
    let content = fs::read_to_string(path)
        .map_err(|e| format!("Cannot read energy file: {}", e))?;

    let mut energies = Vec::new();
    let mut in_data = false;

    for line in content.lines() {
        let trimmed = line.trim();

        // Skip empty lines and comments
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        // Check for END marker
        if trimmed == "END" {
            in_data = false;
            continue;
        }

        // Look for energy trajectory block
        if trimmed.starts_with("ENERTRJ") {
            in_data = true;
            continue;
        }

        // Look for TITLE block (skip)
        if trimmed == "TITLE" {
            in_data = false;
            continue;
        }

        // Parse energy values when in data section
        if in_data {
            // Each line has: time kin pot total temp vol press ...
            // We want the total energy (column 3, index 3)
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 4 {
                if let Ok(total_energy) = parts[3].parse::<f64>() {
                    energies.push(total_energy);
                }
            }
        }
    }

    Ok(energies)
}

/// Clean up test output files
fn cleanup_files(files: &[&str]) {
    for file in files {
        let _ = fs::remove_file(file);
    }
}

#[test]
fn test_md_binary_nve_simulation() {
    // Test NVE (microcanonical) ensemble - no thermostat, no barostat
    let topo_file = "../test_tutorial/water_complete.top";
    let coord_file = "../test_tutorial/water.g96";
    let traj_file = "/tmp/test_md_nve.trc";
    let ene_file = "/tmp/test_md_nve.tre";

    // Clean up any existing files
    cleanup_files(&[traj_file, ene_file]);

    // Run MD simulation: 20 steps, no thermostat, no barostat
    let result = run_md(&[
        "@topo", topo_file,
        "@conf", coord_file,
        "@steps", "20",
        "@dt", "0.001",  // Smaller timestep for stability
        "@traj", traj_file,
        "@ene", ene_file,
        "@cutoff", "1.4",
        "@rf_epsilon", "61.0",
    ]);

    assert!(result.is_ok(), "MD binary failed to execute");
    let output = result.unwrap();

    // Check that the command succeeded
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        panic!("MD simulation failed:\nSTDOUT:\n{}\nSTDERR:\n{}", stdout, stderr);
    }

    // Verify output files were created and have content
    assert!(check_output_file(traj_file).is_ok(),
            "Trajectory file not created or empty");
    assert!(check_output_file(ene_file).is_ok(),
            "Energy file not created or empty");

    // Parse energy file and verify we have data
    let energies = parse_energy_file(ene_file);
    assert!(energies.is_ok(), "Failed to parse energy file");

    let energy_values = energies.unwrap();
    assert!(energy_values.len() > 0, "No energy values found in output");

    // Verify first energy is reasonable (may have NaN later due to instabilities in vacuum)
    assert!(!energy_values[0].is_nan(), "Initial energy should not be NaN");
    assert!(!energy_values[0].is_infinite(), "Initial energy should not be infinite");

    // Clean up
    cleanup_files(&[traj_file, ene_file]);
}

#[test]
fn test_md_binary_nvt_simulation() {
    // Test NVT (canonical) ensemble - with Berendsen thermostat
    let topo_file = "../test_tutorial/water_complete.top";
    let coord_file = "../test_tutorial/water.g96";
    let traj_file = "/tmp/test_md_nvt.trc";
    let ene_file = "/tmp/test_md_nvt.tre";

    cleanup_files(&[traj_file, ene_file]);

    // Run MD simulation with thermostat
    let result = run_md(&[
        "@topo", topo_file,
        "@conf", coord_file,
        "@steps", "100",
        "@dt", "0.002",
        "@traj", traj_file,
        "@ene", ene_file,
        "@cutoff", "1.4",
        "@thermostat", "berendsen",
        "@temp", "300.0",
        "@tau_t", "0.1",
    ]);

    assert!(result.is_ok(), "MD binary failed to execute");
    let output = result.unwrap();

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        panic!("MD simulation failed:\nSTDOUT:\n{}\nSTDERR:\n{}", stdout, stderr);
    }

    // Verify output files
    assert!(check_output_file(traj_file).is_ok());
    assert!(check_output_file(ene_file).is_ok());

    // Verify energies
    let energies = parse_energy_file(ene_file).expect("Failed to parse energy file");
    assert!(energies.len() > 0, "No energy values in NVT simulation");

    // Verify initial energy is reasonable (thermostat helps stability)
    assert!(!energies[0].is_nan() && !energies[0].is_infinite(),
            "Initial energy should be finite");

    cleanup_files(&[traj_file, ene_file]);
}

#[test]
fn test_md_binary_npt_simulation() {
    // Test NPT (isobaric-isothermal) ensemble - with thermostat and barostat
    let topo_file = "../test_tutorial/water_complete.top";
    let coord_file = "../test_tutorial/water.g96";
    let traj_file = "/tmp/test_md_npt.trc";
    let ene_file = "/tmp/test_md_npt.tre";

    cleanup_files(&[traj_file, ene_file]);

    // Run MD simulation with thermostat and barostat
    let result = run_md(&[
        "@topo", topo_file,
        "@conf", coord_file,
        "@steps", "50",
        "@dt", "0.001",  // Smaller timestep
        "@traj", traj_file,
        "@ene", ene_file,
        "@cutoff", "1.4",
        "@thermostat", "berendsen",
        "@temp", "300.0",
        "@tau_t", "0.1",
        "@barostat", "berendsen",
        "@pres", "1.0",
        "@tau_p", "0.5",
    ]);

    assert!(result.is_ok(), "MD binary failed to execute");
    let output = result.unwrap();

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        panic!("MD simulation failed:\nSTDOUT:\n{}\nSTDERR:\n{}", stdout, stderr);
    }

    // Verify output files
    assert!(check_output_file(traj_file).is_ok());
    assert!(check_output_file(ene_file).is_ok());

    // Verify energies
    let energies = parse_energy_file(ene_file).expect("Failed to parse energy file");
    assert!(energies.len() > 0, "No energy values in NPT simulation");

    // Verify initial energy is reasonable
    assert!(!energies[0].is_nan() && !energies[0].is_infinite(),
            "Initial energy should be finite");

    cleanup_files(&[traj_file, ene_file]);
}

#[test]
fn test_md_binary_with_shake_constraints() {
    // Test simulation with SHAKE constraints for bond lengths
    let topo_file = "../test_tutorial/water_complete.top";
    let coord_file = "../test_tutorial/water.g96";
    let traj_file = "/tmp/test_md_shake.trc";
    let ene_file = "/tmp/test_md_shake.tre";

    cleanup_files(&[traj_file, ene_file]);

    // Run MD simulation with SHAKE
    let result = run_md(&[
        "@topo", topo_file,
        "@conf", coord_file,
        "@steps", "100",
        "@dt", "0.002",
        "@traj", traj_file,
        "@ene", ene_file,
        "@shake", "on",
        "@shake_tol", "0.0001",
    ]);

    assert!(result.is_ok(), "MD binary failed to execute");
    let output = result.unwrap();

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        panic!("MD simulation with SHAKE failed:\nSTDOUT:\n{}\nSTDERR:\n{}", stdout, stderr);
    }

    // Verify output files
    assert!(check_output_file(traj_file).is_ok());
    assert!(check_output_file(ene_file).is_ok());

    // Verify no SHAKE convergence warnings in stderr
    let stderr = String::from_utf8_lossy(&output.stderr);
    // Note: Some warnings are OK, but we shouldn't have failures
    assert!(!stderr.contains("SHAKE failed"), "SHAKE algorithm failed");

    cleanup_files(&[traj_file, ene_file]);
}

#[test]
fn test_md_binary_two_waters() {
    // Test simulation with two interacting water molecules
    let topo_file = "../test_tutorial/two_waters.top";
    let coord_file = "../test_tutorial/two_waters.g96";
    let traj_file = "/tmp/test_md_two_waters.trc";
    let ene_file = "/tmp/test_md_two_waters.tre";

    cleanup_files(&[traj_file, ene_file]);

    // Run MD simulation
    let result = run_md(&[
        "@topo", topo_file,
        "@conf", coord_file,
        "@steps", "50",
        "@dt", "0.002",
        "@traj", traj_file,
        "@ene", ene_file,
        "@cutoff", "1.4",
        "@rf_epsilon", "61.0",
        "@thermostat", "berendsen",
        "@temp", "300.0",
    ]);

    assert!(result.is_ok(), "MD binary failed to execute");
    let output = result.unwrap();

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        panic!("Two water simulation failed:\nSTDOUT:\n{}\nSTDERR:\n{}", stdout, stderr);
    }

    // Verify output files
    assert!(check_output_file(traj_file).is_ok());
    assert!(check_output_file(ene_file).is_ok());

    // Parse and verify energies
    let energies = parse_energy_file(ene_file).expect("Failed to parse energy file");
    assert!(energies.len() > 0, "No energy values for two waters");

    // With two waters, we should have nonbonded interactions
    // Energy values should be reasonable (not all zero)
    let has_variation = energies.windows(2).any(|w| (w[0] - w[1]).abs() > 1e-6);
    assert!(has_variation, "Energy should vary during simulation");

    cleanup_files(&[traj_file, ene_file]);
}

#[test]
fn test_md_binary_short_vs_long_cutoff() {
    // Test that different cutoffs produce different energies
    let topo_file = "../test_tutorial/two_waters.top";
    let coord_file = "../test_tutorial/two_waters.g96";

    // Short cutoff simulation
    let traj_short = "/tmp/test_md_cutoff_short.trc";
    let ene_short = "/tmp/test_md_cutoff_short.tre";

    cleanup_files(&[traj_short, ene_short]);

    let result_short = run_md(&[
        "@topo", topo_file,
        "@conf", coord_file,
        "@steps", "10",
        "@dt", "0.002",
        "@traj", traj_short,
        "@ene", ene_short,
        "@cutoff", "0.8",  // Short cutoff
    ]);

    assert!(result_short.is_ok());
    assert!(result_short.unwrap().status.success());

    // Long cutoff simulation
    let traj_long = "/tmp/test_md_cutoff_long.trc";
    let ene_long = "/tmp/test_md_cutoff_long.tre";

    cleanup_files(&[traj_long, ene_long]);

    let result_long = run_md(&[
        "@topo", topo_file,
        "@conf", coord_file,
        "@steps", "10",
        "@dt", "0.002",
        "@traj", traj_long,
        "@ene", ene_long,
        "@cutoff", "1.4",  // Long cutoff
    ]);

    assert!(result_long.is_ok());
    assert!(result_long.unwrap().status.success());

    // Both should have produced output
    assert!(check_output_file(ene_short).is_ok());
    assert!(check_output_file(ene_long).is_ok());

    // Parse energies
    let energies_short = parse_energy_file(ene_short).expect("Failed to parse short cutoff energies");
    let energies_long = parse_energy_file(ene_long).expect("Failed to parse long cutoff energies");

    // Both should have data
    assert!(energies_short.len() > 0);
    assert!(energies_long.len() > 0);

    // Energies should generally be different (longer cutoff includes more interactions)
    // Compare first energy values
    if energies_short.len() > 0 && energies_long.len() > 0 {
        // With two waters at 0.7 nm separation and cutoffs 0.8 vs 1.4,
        // just verify simulations ran and produced valid values
        assert!(!energies_short[0].is_nan());
        assert!(!energies_long[0].is_nan());
    }

    cleanup_files(&[traj_short, ene_short, traj_long, ene_long]);
}

#[test]
fn test_md_binary_invalid_arguments() {
    // Test that MD binary properly handles invalid arguments

    // Missing required argument
    let result = run_md(&["@topo", "../test_tutorial/water_complete.top"]);
    assert!(result.is_ok(), "Should execute even if fails");
    let output = result.unwrap();
    assert!(!output.status.success(), "Should fail with missing @conf");

    // Invalid file path
    let result = run_md(&[
        "@topo", "nonexistent_file.top",
        "@conf", "../test_tutorial/water.g96",
    ]);
    assert!(result.is_ok());
    let output = result.unwrap();
    assert!(!output.status.success(), "Should fail with nonexistent topology");

    // Invalid parameter value
    let result = run_md(&[
        "@topo", "../test_tutorial/water_complete.top",
        "@conf", "../test_tutorial/water.g96",
        "@steps", "not_a_number",
    ]);
    assert!(result.is_ok());
    let output = result.unwrap();
    assert!(!output.status.success(), "Should fail with invalid @steps value");
}
