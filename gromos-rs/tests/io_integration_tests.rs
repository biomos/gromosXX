//! Integration tests for GROMOS I/O modules
//!
//! Tests the complete I/O workflow:
//! - IMD file reading and parsing
//! - Trajectory writing (TRC)
//! - Energy writing (TRE)
//! - Force writing (TRF)
//! - mk_script binary integration

use gromos_rs::io::{
    imd, trajectory, energy, force,
    ImdParameters, TrajectoryWriter, EnergyWriter, EnergyFrame, ForceWriter,
};
use gromos_rs::configuration::{Configuration, Box as SimBox};
use gromos_rs::math::Vec3;
use std::fs;
use std::path::Path;

#[test]
fn test_imd_write_and_read() {
    let test_file = "/tmp/test_imd_integration.imd";

    // Write a test IMD file manually
    let imd_content = r#"TITLE
  Integration test simulation
END

SYSTEM
#  NPM    NSM
     1      0
END

STEP
#  NSTLIM      T        DT
    1000    0.0     0.00200
END

BOUNDCOND
#   NTB  NDFMIN
      1       0
END

MULTIBATH
#  ALGORITHM
          1
#  NUM
      1
#  TEMP0         TAU
  300.00      0.1000
#  DOFSET: num, last_atom_index
      1     0
END

CONSTRAINT
#  NTC  NTCP  NTCP0(1)  NTCS  NTCS0(1)
    2     0         0      1         0
END

PAIRLIST
#  ALGORITHM  NSNB  RCUTP   RCUTL     SIZE  TYPE
          0     5   0.80    1.40     0.40     0
END

NONBONDED
#  NLRELE  APPAK    RCRF   EPSRF  NSLFEXCL
       1    0.0    1.40     0.0         1
END

INITIALISE
#  NTIVEL  NTISHK  NTINHT  NTINHB  NTISHI     NTIRTC  NTICOM
        1       0       0       0       0          0       0
#  NTIR    NTIG      IG     TEMPI
      0       0   12345  300.00
END

WRITETRAJ
  NTWX     100
  NTWE      10
  NTWV       0
  NTWF       0
END

PRINTOUT
#  NTPR
     10
END
"#;

    fs::write(test_file, imd_content).expect("Failed to write test IMD file");

    // Read and parse the IMD file
    let params = imd::read_imd_file(test_file).expect("Failed to read IMD file");

    // Validate parsed parameters
    assert_eq!(params.title, "Integration test simulation");
    assert_eq!(params.npm, 1);
    assert_eq!(params.nsm, 0);
    assert_eq!(params.nstlim, 1000);
    assert_eq!(params.dt, 0.002);
    assert_eq!(params.ntc, 2);
    assert_eq!(params.nlrele, 1);
    assert_eq!(params.rcutl, 1.4);
    assert_eq!(params.ntwx, 100);
    assert_eq!(params.ntwe, 10);
    assert_eq!(params.tempi, 300.0);

    // Cleanup
    fs::remove_file(test_file).ok();
}

#[test]
fn test_trajectory_writer_integration() {
    let test_file = "/tmp/test_trajectory_integration.trc";

    // Create a simple configuration
    let mut config = Configuration::new(3, 1, 1);
    config.current_mut().pos = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.1, 0.0, 0.0),
        Vec3::new(0.0, 0.1, 0.0),
    ];
    config.current_mut().vel = vec![
        Vec3::new(0.01, 0.0, 0.0),
        Vec3::new(0.0, 0.01, 0.0),
        Vec3::new(0.0, 0.0, 0.01),
    ];
    config.current_mut().box_config = SimBox::rectangular(3.0, 3.0, 3.0);

    // Write multiple frames
    let mut writer = TrajectoryWriter::new(
        test_file,
        "Integration test trajectory",
        true,  // Write velocities
        false,
    ).expect("Failed to create trajectory writer");

    for step in 0..5 {
        let time = step as f64 * 0.002;
        writer.write_frame(step, time, &config).expect("Failed to write frame");
    }

    writer.flush().expect("Failed to flush trajectory");

    // Verify file was created and has content
    assert!(Path::new(test_file).exists());
    let content = fs::read_to_string(test_file).expect("Failed to read trajectory file");

    assert!(content.contains("TITLE"));
    assert!(content.contains("Integration test trajectory"));
    assert!(content.contains("TIMESTEP"));
    assert!(content.contains("POSITIONRED"));
    assert!(content.contains("VELOCITYRED"));
    assert!(content.contains("GENBOX"));

    // Should have 5 frames
    assert_eq!(content.matches("TIMESTEP").count(), 5);
    assert_eq!(writer.frame_count(), 5);

    // Cleanup
    fs::remove_file(test_file).ok();
}

#[test]
fn test_energy_writer_integration() {
    let test_file = "/tmp/test_energy_integration.tre";

    let mut writer = EnergyWriter::new(test_file, "Integration test energy")
        .expect("Failed to create energy writer");

    // Write multiple energy frames
    for step in 0..10 {
        let time = step as f64 * 0.002;
        let mut frame = EnergyFrame::new(
            time,
            100.0 + step as f64,      // Kinetic
            -200.0 - step as f64,     // Potential
            300.0,                     // Temperature
        );

        frame.bond = -50.0;
        frame.angle = -30.0;
        frame.lj = -80.0;
        frame.coul_real = -40.0;
        frame.update_potential();

        writer.write_frame(&frame).expect("Failed to write energy frame");
    }

    writer.finalize().expect("Failed to finalize energy file");

    // Verify file contents
    assert!(Path::new(test_file).exists());
    let content = fs::read_to_string(test_file).expect("Failed to read energy file");

    assert!(content.contains("TITLE"));
    assert!(content.contains("Integration test energy"));
    assert!(content.contains("ENERTRJ"));
    assert!(content.contains("Block definitions"));
    assert!(content.contains("Kinetic energy"));

    // Should have 10 data lines (excluding header/comments/title text)
    let data_lines = content.lines()
        .filter(|l| {
            let trimmed = l.trim();
            !trimmed.is_empty() &&
            !l.starts_with('#') &&
            !l.contains("TITLE") &&
            !l.contains("END") &&
            !l.contains("ENERTRJ") &&
            !l.contains("Integration test energy")
        })
        .count();
    assert_eq!(data_lines, 10);

    // Cleanup
    fs::remove_file(test_file).ok();
}

#[test]
fn test_force_writer_integration() {
    let test_file = "/tmp/test_force_integration.trf";

    let forces = vec![
        Vec3::new(1.0, 2.0, 3.0),
        Vec3::new(4.0, 5.0, 6.0),
        Vec3::new(7.0, 8.0, 9.0),
    ];

    let mut writer = ForceWriter::new(test_file, "Integration test forces")
        .expect("Failed to create force writer");

    // Write multiple force frames
    for step in 0..3 {
        let time = step as f64 * 0.002;
        writer.write_frame(step, time, &forces).expect("Failed to write force frame");
    }

    writer.flush().expect("Failed to flush forces");

    // Verify file contents
    assert!(Path::new(test_file).exists());
    let content = fs::read_to_string(test_file).expect("Failed to read force file");

    assert!(content.contains("TITLE"));
    assert!(content.contains("Integration test forces"));
    assert!(content.contains("TIMESTEP"));
    assert!(content.contains("FORCE"));
    assert_eq!(content.matches("TIMESTEP").count(), 3);

    // Cleanup
    fs::remove_file(test_file).ok();
}

#[test]
fn test_force_writer_detailed() {
    let test_file = "/tmp/test_force_detailed.trf";

    let n_atoms = 3;
    let bonded = vec![Vec3::new(1.0, 0.0, 0.0); n_atoms];
    let nonbonded = vec![Vec3::new(0.0, 1.0, 0.0); n_atoms];
    let constraint = vec![Vec3::new(0.0, 0.0, 1.0); n_atoms];

    let mut writer = ForceWriter::new(test_file, "Detailed force test")
        .expect("Failed to create force writer");

    writer.write_frame_detailed(0, 0.0, &bonded, &nonbonded, &constraint)
        .expect("Failed to write detailed frame");

    writer.flush().expect("Failed to flush");

    // Verify detailed breakdown
    let content = fs::read_to_string(test_file).expect("Failed to read force file");
    assert!(content.contains("FORCE"));
    assert!(content.contains("FORCE_BONDED"));
    assert!(content.contains("FORCE_NONBONDED"));
    assert!(content.contains("FORCE_CONSTRAINT"));

    // Cleanup
    fs::remove_file(test_file).ok();
}

#[test]
fn test_complete_io_workflow() {
    // Simulate a complete MD workflow
    let imd_file = "/tmp/workflow_test.imd";
    let trc_file = "/tmp/workflow_test.trc";
    let tre_file = "/tmp/workflow_test.tre";

    // 1. Write IMD file
    let imd_content = r#"TITLE
  Complete workflow test
END
SYSTEM
  NPM    NSM
     1      0
END
STEP
  NSTLIM    100
  T         0.0
  DT        0.002
END
BOUNDCOND
   NTB  NDFMIN
      1       0
END
WRITETRAJ
  NTWX       5
  NTWE       5
  NTWV       0
  NTWF       0
END
PRINTOUT
  NTPR
     5
END
"#;
    fs::write(imd_file, imd_content).expect("Failed to write IMD");

    // 2. Read IMD parameters
    let params = imd::read_imd_file(imd_file).expect("Failed to read IMD");
    assert_eq!(params.nstlim, 100);
    assert_eq!(params.ntwx, 5);
    assert_eq!(params.ntwe, 5);

    // 3. Create configuration
    let mut config = Configuration::new(10, 1, 1);
    for i in 0..10 {
        config.current_mut().pos[i] = Vec3::new(i as f32 * 0.1, 0.0, 0.0);
    }
    config.current_mut().box_config = SimBox::rectangular(5.0, 5.0, 5.0);

    // 4. Write trajectory
    let mut traj = TrajectoryWriter::new(trc_file, &params.title, false, false)
        .expect("Failed to create trajectory");

    // 5. Write energy
    let mut energy = EnergyWriter::new(tre_file, &params.title)
        .expect("Failed to create energy writer");

    // 6. Simulate MD steps
    for step in 0..params.nstlim {
        let time = step as f64 * params.dt;

        // Write trajectory every ntwx steps
        if step % params.ntwx == 0 {
            traj.write_frame(step, time, &config).expect("Failed to write trajectory");
        }

        // Write energy every ntwe steps
        if step % params.ntwe == 0 {
            let frame = EnergyFrame::new(time, 100.0, -200.0, 300.0);
            energy.write_frame(&frame).expect("Failed to write energy");
        }
    }

    traj.flush().expect("Failed to flush trajectory");
    energy.finalize().expect("Failed to finalize energy");

    // 7. Verify outputs
    assert!(Path::new(imd_file).exists());
    assert!(Path::new(trc_file).exists());
    assert!(Path::new(tre_file).exists());

    // Trajectory should have 100/5 = 20 frames
    assert_eq!(traj.frame_count(), 20);

    // Energy should have 100/5 = 20 frames
    assert_eq!(energy.frame_count(), 20);

    // Cleanup
    fs::remove_file(imd_file).ok();
    fs::remove_file(trc_file).ok();
    fs::remove_file(tre_file).ok();
}

#[test]
fn test_imd_parameter_defaults() {
    let params = ImdParameters::default();

    // Check GROMOS defaults
    assert_eq!(params.nstlim, 1000);
    assert_eq!(params.dt, 0.002);
    assert_eq!(params.ntc, 2);  // SHAKE on H-bonds
    assert_eq!(params.nlrele, 1);  // Reaction Field
    assert_eq!(params.epsrf, 0.0);  // Conducting boundary
    assert_eq!(params.tempi, 300.0);
}

#[test]
fn test_energy_frame_updates() {
    let mut frame = EnergyFrame::default();

    // Set components
    frame.kinetic = 100.0;
    frame.bond = -50.0;
    frame.angle = -30.0;
    frame.lj = -80.0;
    frame.coul_real = -40.0;

    // Update totals
    frame.update_potential();
    frame.update_total();

    assert_eq!(frame.potential, -200.0);
    assert_eq!(frame.total, -100.0);
}
