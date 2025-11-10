//! GROMOS force file writer (.trf)
//!
//! Writes atomic forces at regular intervals (typically for analysis or debugging).
//!
//! # Format
//! ```text
//! TITLE
//!   Force trajectory
//! END
//! TIMESTEP
//!        1     0.000000
//! END
//! FORCE
//!     1 RES    ATOM      1    fx       fy       fz
//!     2 RES    ATOM      2    fx       fy       fz
//! END
//! ```

use crate::math::Vec3;
use crate::io::IoError;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;

/// GROMOS force file writer
pub struct ForceWriter {
    writer: BufWriter<File>,
    frame_count: usize,
}

impl ForceWriter {
    /// Create a new force writer
    pub fn new<P: AsRef<Path>>(path: P, title: &str) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write TITLE block
        writeln!(writer, "TITLE")?;
        writeln!(writer, "  {}", title)?;
        writeln!(writer, "END")?;

        Ok(Self {
            writer,
            frame_count: 0,
        })
    }

    /// Write a force frame
    pub fn write_frame(
        &mut self,
        step: usize,
        time: f64,
        forces: &[Vec3],
    ) -> Result<(), IoError> {
        // TIMESTEP block
        writeln!(self.writer, "TIMESTEP")?;
        writeln!(self.writer, "{:>10} {:15.6}", step, time)?;
        writeln!(self.writer, "END")?;

        // FORCE block
        writeln!(self.writer, "FORCE")?;
        for (i, force) in forces.iter().enumerate() {
            writeln!(
                self.writer,
                "{:>6} {:>6} {:>6} {:>6} {:15.8} {:15.8} {:15.8}",
                i + 1,         // Atom number (1-indexed)
                "RES",         // Residue name (placeholder)
                "ATOM",        // Atom name (placeholder)
                i + 1,         // Atom serial
                force.x,
                force.y,
                force.z
            )?;
        }
        writeln!(self.writer, "END")?;

        self.frame_count += 1;
        Ok(())
    }

    /// Write force components (bonded, nonbonded, constraint)
    ///
    /// This provides a detailed breakdown of forces by type
    pub fn write_frame_detailed(
        &mut self,
        step: usize,
        time: f64,
        bonded_forces: &[Vec3],
        nonbonded_forces: &[Vec3],
        constraint_forces: &[Vec3],
    ) -> Result<(), IoError> {
        // TIMESTEP block
        writeln!(self.writer, "TIMESTEP")?;
        writeln!(self.writer, "{:>10} {:15.6}", step, time)?;
        writeln!(self.writer, "END")?;

        // Total forces
        let n_atoms = bonded_forces.len();
        let mut total_forces = vec![Vec3::ZERO; n_atoms];
        for i in 0..n_atoms {
            total_forces[i] = bonded_forces[i] + nonbonded_forces[i] + constraint_forces[i];
        }

        writeln!(self.writer, "FORCE")?;
        for (i, force) in total_forces.iter().enumerate() {
            writeln!(
                self.writer,
                "{:>6} {:>6} {:>6} {:>6} {:15.8} {:15.8} {:15.8}",
                i + 1, "RES", "ATOM", i + 1,
                force.x, force.y, force.z
            )?;
        }
        writeln!(self.writer, "END")?;

        // Bonded forces breakdown
        writeln!(self.writer, "FORCE_BONDED")?;
        for (i, force) in bonded_forces.iter().enumerate() {
            writeln!(
                self.writer,
                "{:>6} {:>6} {:>6} {:>6} {:15.8} {:15.8} {:15.8}",
                i + 1, "RES", "ATOM", i + 1,
                force.x, force.y, force.z
            )?;
        }
        writeln!(self.writer, "END")?;

        // Nonbonded forces breakdown
        writeln!(self.writer, "FORCE_NONBONDED")?;
        for (i, force) in nonbonded_forces.iter().enumerate() {
            writeln!(
                self.writer,
                "{:>6} {:>6} {:>6} {:>6} {:15.8} {:15.8} {:15.8}",
                i + 1, "RES", "ATOM", i + 1,
                force.x, force.y, force.z
            )?;
        }
        writeln!(self.writer, "END")?;

        // Constraint forces breakdown
        writeln!(self.writer, "FORCE_CONSTRAINT")?;
        for (i, force) in constraint_forces.iter().enumerate() {
            writeln!(
                self.writer,
                "{:>6} {:>6} {:>6} {:>6} {:15.8} {:15.8} {:15.8}",
                i + 1, "RES", "ATOM", i + 1,
                force.x, force.y, force.z
            )?;
        }
        writeln!(self.writer, "END")?;

        self.frame_count += 1;
        Ok(())
    }

    /// Flush buffered data
    pub fn flush(&mut self) -> Result<(), IoError> {
        self.writer.flush()?;
        Ok(())
    }

    /// Get number of frames written
    pub fn frame_count(&self) -> usize {
        self.frame_count
    }
}

/// Write a single force frame to a file (convenience function)
pub fn write_force_frame<P: AsRef<Path>>(
    path: P,
    step: usize,
    time: f64,
    forces: &[Vec3],
    title: &str,
) -> Result<(), IoError> {
    let mut writer = ForceWriter::new(path, title)?;
    writer.write_frame(step, time, forces)?;
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_force_writer_creation() {
        let temp_file = "/tmp/test_force.trf";
        let writer = ForceWriter::new(temp_file, "Test forces");
        assert!(writer.is_ok());
    }

    #[test]
    fn test_write_force_frame() {
        let temp_file = "/tmp/test_force2.trf";
        let mut writer = ForceWriter::new(temp_file, "Test").unwrap();

        let forces = vec![
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(4.0, 5.0, 6.0),
            Vec3::new(7.0, 8.0, 9.0),
        ];

        let result = writer.write_frame(0, 0.0, &forces);
        assert!(result.is_ok());
        assert_eq!(writer.frame_count(), 1);
    }
}
