//! Integration tests for mk_script binary
//!
//! Tests the mk_script command-line tool functionality

use std::process::Command;
use std::fs;
use std::path::Path;

fn get_mk_script_path() -> String {
    // Check if we're in debug or release mode
    if Path::new("target/release/mk_script").exists() {
        "target/release/mk_script".to_string()
    } else if Path::new("target/debug/mk_script").exists() {
        "target/debug/mk_script".to_string()
    } else {
        // Try to build it
        let output = Command::new("cargo")
            .args(&["build", "--bin", "mk_script"])
            .output()
            .expect("Failed to build mk_script");

        if !output.status.success() {
            panic!("Failed to build mk_script: {}", String::from_utf8_lossy(&output.stderr));
        }

        "target/debug/mk_script".to_string()
    }
}

#[test]
fn test_mk_script_help() {
    let mk_script = get_mk_script_path();

    let output = Command::new(&mk_script)
        .arg("--help")
        .output()
        .expect("Failed to run mk_script --help");

    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Generate GROMOS simulation scripts"));
    assert!(stdout.contains("--topo"));
    assert!(stdout.contains("--coord"));
    assert!(stdout.contains("--output"));
}

#[test]
fn test_mk_script_basic_generation() {
    let mk_script = get_mk_script_path();

    // Create dummy input files
    fs::write("/tmp/test_mk.top", "").expect("Failed to create topology");
    fs::write("/tmp/test_mk.cnf", "").expect("Failed to create coordinates");

    let output = Command::new(&mk_script)
        .args(&[
            "--topo", "/tmp/test_mk.top",
            "--coord", "/tmp/test_mk.cnf",
            "--output", "/tmp/test_mk.imd",
        ])
        .output()
        .expect("Failed to run mk_script");

    assert!(output.status.success(),
            "mk_script failed: {}", String::from_utf8_lossy(&output.stderr));

    // Verify IMD file was created
    assert!(Path::new("/tmp/test_mk.imd").exists());

    // Read and verify IMD content
    let imd_content = fs::read_to_string("/tmp/test_mk.imd")
        .expect("Failed to read generated IMD");

    assert!(imd_content.contains("TITLE"));
    assert!(imd_content.contains("SYSTEM"));
    assert!(imd_content.contains("STEP"));
    assert!(imd_content.contains("BOUNDCOND"));
    assert!(imd_content.contains("MULTIBATH"));
    assert!(imd_content.contains("CONSTRAINT"));
    assert!(imd_content.contains("PAIRLIST"));
    assert!(imd_content.contains("NONBONDED"));
    assert!(imd_content.contains("INITIALISE"));
    assert!(imd_content.contains("WRITETRAJ"));
    assert!(imd_content.contains("PRINTOUT"));
    assert!(imd_content.contains("END"));

    // Cleanup
    fs::remove_file("/tmp/test_mk.top").ok();
    fs::remove_file("/tmp/test_mk.cnf").ok();
    fs::remove_file("/tmp/test_mk.imd").ok();
}

#[test]
fn test_mk_script_custom_parameters() {
    let mk_script = get_mk_script_path();

    fs::write("/tmp/test_custom.top", "").ok();
    fs::write("/tmp/test_custom.cnf", "").ok();

    let output = Command::new(&mk_script)
        .args(&[
            "--topo", "/tmp/test_custom.top",
            "--coord", "/tmp/test_custom.cnf",
            "--output", "/tmp/test_custom.imd",
            "--title", "Custom test",
            "--nstlim", "5000",
            "--dt", "0.001",
            "--temp", "310.0",
            "--nlrele", "2",  // PME
            "--ntc", "3",     // All constraints
        ])
        .output()
        .expect("Failed to run mk_script");

    assert!(output.status.success());

    let imd_content = fs::read_to_string("/tmp/test_custom.imd")
        .expect("Failed to read IMD");

    // Verify custom parameters
    assert!(imd_content.contains("Custom test"));
    assert!(imd_content.contains("5000"));       // NSTLIM
    assert!(imd_content.contains("0.00100"));    // DT
    assert!(imd_content.contains("310.00"));     // TEMP
    assert!(imd_content.contains("NLRELE"));

    // Check for NTC = 3
    let ntc_line = imd_content.lines()
        .find(|l| l.contains("NTC") && !l.starts_with('#'))
        .expect("NTC line not found");
    assert!(ntc_line.contains("3"));

    // Cleanup
    fs::remove_file("/tmp/test_custom.top").ok();
    fs::remove_file("/tmp/test_custom.cnf").ok();
    fs::remove_file("/tmp/test_custom.imd").ok();
}

#[test]
fn test_mk_script_pressure_coupling() {
    let mk_script = get_mk_script_path();

    fs::write("/tmp/test_press.top", "").ok();
    fs::write("/tmp/test_press.cnf", "").ok();

    let output = Command::new(&mk_script)
        .args(&[
            "--topo", "/tmp/test_press.top",
            "--coord", "/tmp/test_press.cnf",
            "--output", "/tmp/test_press.imd",
            "--pcouple", "1",      // Isotropic
            "--pres", "1.0",
            "--tau-p", "0.5",
        ])
        .output()
        .expect("Failed to run mk_script");

    assert!(output.status.success());

    let imd_content = fs::read_to_string("/tmp/test_press.imd")
        .expect("Failed to read IMD");

    // Should have PRESSURESCALE block
    assert!(imd_content.contains("PRESSURESCALE"));
    assert!(imd_content.contains("COUPLE"));

    // Cleanup
    fs::remove_file("/tmp/test_press.top").ok();
    fs::remove_file("/tmp/test_press.cnf").ok();
    fs::remove_file("/tmp/test_press.imd").ok();
}

#[test]
fn test_mk_script_missing_input_files() {
    let mk_script = get_mk_script_path();

    let output = Command::new(&mk_script)
        .args(&[
            "--topo", "/tmp/nonexistent.top",
            "--coord", "/tmp/nonexistent.cnf",
            "--output", "/tmp/test_err.imd",
        ])
        .output()
        .expect("Failed to run mk_script");

    // Should fail with error
    assert!(!output.status.success());

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("not found") || stderr.contains("Error"));
}

#[test]
fn test_mk_script_reaction_field_settings() {
    let mk_script = get_mk_script_path();

    fs::write("/tmp/test_rf.top", "").ok();
    fs::write("/tmp/test_rf.cnf", "").ok();

    let output = Command::new(&mk_script)
        .args(&[
            "--topo", "/tmp/test_rf.top",
            "--coord", "/tmp/test_rf.cnf",
            "--output", "/tmp/test_rf.imd",
            "--nlrele", "1",       // Reaction Field
            "--epsrf", "78.4",     // Water
            "--rcutl", "1.2",
        ])
        .output()
        .expect("Failed to run mk_script");

    assert!(output.status.success());

    let imd_content = fs::read_to_string("/tmp/test_rf.imd")
        .expect("Failed to read IMD");

    // Check RF parameters
    assert!(imd_content.contains("NLRELE"));
    assert!(imd_content.contains("78.4"));
    assert!(imd_content.contains("1.20"));

    // Cleanup
    fs::remove_file("/tmp/test_rf.top").ok();
    fs::remove_file("/tmp/test_rf.cnf").ok();
    fs::remove_file("/tmp/test_rf.imd").ok();
}

#[test]
fn test_mk_script_output_frequencies() {
    let mk_script = get_mk_script_path();

    fs::write("/tmp/test_freq.top", "").ok();
    fs::write("/tmp/test_freq.cnf", "").ok();

    let output = Command::new(&mk_script)
        .args(&[
            "--topo", "/tmp/test_freq.top",
            "--coord", "/tmp/test_freq.cnf",
            "--output", "/tmp/test_freq.imd",
            "--ntwx", "250",
            "--ntwe", "50",
        ])
        .output()
        .expect("Failed to run mk_script");

    assert!(output.status.success());

    let imd_content = fs::read_to_string("/tmp/test_freq.imd")
        .expect("Failed to read IMD");

    // Check WRITETRAJ block
    assert!(imd_content.contains("WRITETRAJ"));
    assert!(imd_content.contains("250"));
    assert!(imd_content.contains("50"));

    // Cleanup
    fs::remove_file("/tmp/test_freq.top").ok();
    fs::remove_file("/tmp/test_freq.cnf").ok();
    fs::remove_file("/tmp/test_freq.imd").ok();
}
