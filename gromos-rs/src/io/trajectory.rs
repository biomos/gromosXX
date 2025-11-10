//! GROMOS trajectory file writer (.trc/.trj)
//!
//! Writes coordinates (and optionally velocities and forces) at regular intervals.
//!
//! # Format
//! ```text
//! TITLE
//!   System trajectory
//! END
//! TIMESTEP
//!        1     0.000000
//! END
//! POSITIONRED
//!     1 RES    ATOM      1    x        y        z
//!     2 RES    ATOM      2    x        y        z
//! END
//! VELOCITYRED (optional)
//!     1 RES    ATOM      1    vx       vy       vz
//! END
//! FORCERED (optional, added in GROMOS-RS)
//!     1 RES    ATOM      1    fx       fy       fz
//! #         10  (comment every 10 atoms)
//! END
//! LATTICESHIFTS (optional)
//!     1    0    0    0
//! END
//! GENBOX
//!     lx   ly   lz
//! END
//! ```
//!
//! Note: For dedicated force trajectory output with FREEFORCERED/CONSFORCERED blocks,
//! use the `ForceWriter` from `io::force` module instead.

use crate::configuration::{Configuration, Box as SimBox};
use crate::math::Vec3;
use crate::io::IoError;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;

/// GROMOS trajectory writer
pub struct TrajectoryWriter {
    writer: BufWriter<File>,
    write_velocities: bool,
    write_forces: bool,
    step_count: usize,
}

impl TrajectoryWriter {
    /// Create a new trajectory writer
    pub fn new<P: AsRef<Path>>(
        path: P,
        title: &str,
        write_velocities: bool,
        write_forces: bool,
    ) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(writer, "TITLE")?;
        writeln!(writer, "  {}", title)?;
        writeln!(writer, "END")?;

        Ok(Self {
            writer,
            write_velocities,
            write_forces,
            step_count: 0,
        })
    }

    /// Write a trajectory frame
    pub fn write_frame(
        &mut self,
        step: usize,
        time: f64,
        config: &Configuration,
    ) -> Result<(), IoError> {
        let state = config.current();

        // TIMESTEP block
        writeln!(self.writer, "TIMESTEP")?;
        writeln!(self.writer, "{:>10} {:15.6}", step, time)?;
        writeln!(self.writer, "END")?;

        // POSITIONRED block (reduced precision for smaller files)
        writeln!(self.writer, "POSITIONRED")?;
        for (i, pos) in state.pos.iter().enumerate() {
            writeln!(
                self.writer,
                "{:>6} {:>6} {:>6} {:>6} {:12.6} {:12.6} {:12.6}",
                i + 1,         // Atom number (1-indexed)
                "RES",         // Residue name (placeholder)
                "ATOM",        // Atom name (placeholder)
                i + 1,         // Atom serial
                pos.x,
                pos.y,
                pos.z
            )?;
        }
        writeln!(self.writer, "END")?;

        // VELOCITYRED block (optional)
        if self.write_velocities {
            writeln!(self.writer, "VELOCITYRED")?;
            for (i, vel) in state.vel.iter().enumerate() {
                writeln!(
                    self.writer,
                    "{:>6} {:>6} {:>6} {:>6} {:12.6} {:12.6} {:12.6}",
                    i + 1,
                    "RES",
                    "ATOM",
                    i + 1,
                    vel.x,
                    vel.y,
                    vel.z
                )?;
            }
            writeln!(self.writer, "END")?;
        }

        // FORCERED block (optional) - uses higher precision for forces
        if self.write_forces {
            writeln!(self.writer, "FORCERED")?;
            for (i, force) in state.force.iter().enumerate() {
                writeln!(
                    self.writer,
                    "{:>6} {:>6} {:>6} {:>6} {:18.9} {:18.9} {:18.9}",
                    i + 1,
                    "RES",
                    "ATOM",
                    i + 1,
                    force.x,
                    force.y,
                    force.z
                )?;
                // Add comment every 10 atoms (GROMOS convention)
                if (i + 1) % 10 == 0 {
                    writeln!(self.writer, "#{:>10}", i + 1)?;
                }
            }
            writeln!(self.writer, "END")?;
        }

        // GENBOX block (box dimensions)
        let dims = state.box_config.dimensions();
        writeln!(self.writer, "GENBOX")?;
        writeln!(
            self.writer,
            " {:15.9} {:15.9} {:15.9}",
            dims.x,
            dims.y,
            dims.z
        )?;
        writeln!(self.writer, "END")?;

        self.step_count += 1;
        Ok(())
    }

    /// Flush buffered data to disk
    pub fn flush(&mut self) -> Result<(), IoError> {
        self.writer.flush()?;
        Ok(())
    }

    /// Get number of frames written
    pub fn frame_count(&self) -> usize {
        self.step_count
    }
}

/// Write a single trajectory frame to a file (convenience function)
pub fn write_trajectory_frame<P: AsRef<Path>>(
    path: P,
    step: usize,
    time: f64,
    config: &Configuration,
    title: &str,
) -> Result<(), IoError> {
    let mut writer = TrajectoryWriter::new(path, title, false, false)?;
    writer.write_frame(step, time, config)?;
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::configuration::Box as SimBox;

    #[test]
    fn test_trajectory_writer_creation() {
        let temp_file = "/tmp/test_traj.trc";
        let writer = TrajectoryWriter::new(temp_file, "Test trajectory", false, false);
        assert!(writer.is_ok());
    }

    #[test]
    fn test_write_frame() {
        let temp_file = "/tmp/test_traj2.trc";
        let mut config = Configuration::new(3, 1, 1);
        config.current_mut().pos = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::new(0.0, 0.1, 0.0),
        ];
        config.current_mut().box_config = SimBox::rectangular(3.0, 3.0, 3.0);

        let mut writer = TrajectoryWriter::new(temp_file, "Test", false, false).unwrap();
        let result = writer.write_frame(0, 0.0, &config);
        assert!(result.is_ok());
        assert_eq!(writer.frame_count(), 1);
    }

    #[test]
    fn test_write_frame_with_forces() {
        let temp_file = "/tmp/test_traj_forces.trc";
        let mut config = Configuration::new(3, 1, 1);
        config.current_mut().pos = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::new(0.0, 0.1, 0.0),
        ];
        config.current_mut().force = vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        config.current_mut().box_config = SimBox::rectangular(3.0, 3.0, 3.0);

        let mut writer = TrajectoryWriter::new(temp_file, "Test with forces", false, true).unwrap();
        let result = writer.write_frame(0, 0.0, &config);
        assert!(result.is_ok());
        assert_eq!(writer.frame_count(), 1);
    }
}
