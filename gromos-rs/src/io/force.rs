//! GROMOS force trajectory file writer (.trf)
//!
//! Writes atomic forces at regular intervals matching the GROMOS md++ TRF format.
//!
//! # Real GROMOS TRF Format
//! ```text
//! TITLE
//!   Force trajectory
//! END
//! TIMESTEP
//!        1     0.000000
//! END
//! FREEFORCERED
//!    fx1       fy1       fz1
//!    fx2       fy2       fz2
//! #         10  (comment every 10 atoms)
//!    ...
//! END
//! CONSFORCERED (optional, constraint forces)
//!    cfx1      cfy1      cfz1
//!    ...
//! END
//! ```
//!
//! This matches the format from md++/src/io/configuration/out_configuration.cc

use crate::math::Vec3;
use crate::io::IoError;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;

/// GROMOS force trajectory writer (.trf files)
///
/// Matches the official GROMOS md++ format with FREEFORCERED and optional CONSFORCERED blocks
pub struct ForceWriter {
    writer: BufWriter<File>,
    write_constraint_forces: bool,
    frame_count: usize,
    force_precision: usize,  // Default: 9 decimals (matches md++)
    force_width: usize,      // Default: 18 characters (matches md++)
}

impl ForceWriter {
    /// Create a new force trajectory writer
    ///
    /// # Arguments
    /// * `path` - Output file path (.trf extension recommended)
    /// * `title` - Trajectory title
    /// * `write_constraint_forces` - Whether to write CONSFORCERED block
    pub fn new<P: AsRef<Path>>(
        path: P,
        title: &str,
        write_constraint_forces: bool,
    ) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write TITLE block
        writeln!(writer, "TITLE")?;
        writeln!(writer, "  {}", title)?;
        writeln!(writer, "END")?;

        Ok(Self {
            writer,
            write_constraint_forces,
            frame_count: 0,
            force_precision: 9,   // Match md++ default
            force_width: 18,      // Match md++ default
        })
    }

    /// Create a force writer with default settings (no constraint forces)
    pub fn new_simple<P: AsRef<Path>>(path: P, title: &str) -> Result<Self, IoError> {
        Self::new(path, title, false)
    }

    /// Write a force trajectory frame
    ///
    /// # Arguments
    /// * `step` - MD step number
    /// * `time` - Simulation time (ps)
    /// * `forces` - Free forces on all atoms
    /// * `constraint_forces` - Optional constraint forces (required if write_constraint_forces is true)
    pub fn write_frame(
        &mut self,
        step: usize,
        time: f64,
        forces: &[Vec3],
        constraint_forces: Option<&[Vec3]>,
    ) -> Result<(), IoError> {
        // TIMESTEP block
        writeln!(self.writer, "TIMESTEP")?;
        writeln!(self.writer, "{:>10} {:15.6}", step, time)?;
        writeln!(self.writer, "END")?;

        // FREEFORCERED block (normal forces)
        writeln!(self.writer, "FREEFORCERED")?;

        for (i, force) in forces.iter().enumerate() {
            writeln!(
                self.writer,
                "{:width$.prec$} {:width$.prec$} {:width$.prec$}",
                force.x,
                force.y,
                force.z,
                width = self.force_width,
                prec = self.force_precision
            )?;

            // Add comment every 10 atoms (GROMOS md++ convention)
            if (i + 1) % 10 == 0 {
                writeln!(self.writer, "#{:>10}", i + 1)?;
            }
        }
        writeln!(self.writer, "END")?;

        // CONSFORCERED block (constraint forces, optional)
        if self.write_constraint_forces {
            if let Some(cons_forces) = constraint_forces {
                writeln!(self.writer, "CONSFORCERED")?;

                for (i, force) in cons_forces.iter().enumerate() {
                    writeln!(
                        self.writer,
                        "{:width$.prec$} {:width$.prec$} {:width$.prec$}",
                        force.x,
                        force.y,
                        force.z,
                        width = self.force_width,
                        prec = self.force_precision
                    )?;

                    if (i + 1) % 10 == 0 {
                        writeln!(self.writer, "#{:>10}", i + 1)?;
                    }
                }
                writeln!(self.writer, "END")?;
            }
        }

        self.frame_count += 1;
        Ok(())
    }

    /// Write frame with only free forces (no constraints)
    pub fn write_frame_simple(
        &mut self,
        step: usize,
        time: f64,
        forces: &[Vec3],
    ) -> Result<(), IoError> {
        self.write_frame(step, time, forces, None)
    }

    /// Set custom force output precision and width
    ///
    /// Default is precision=9, width=18 (matches GROMOS md++)
    pub fn set_precision(&mut self, precision: usize, width: usize) {
        self.force_precision = precision;
        self.force_width = width;
    }

    /// Flush buffered data to disk
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
    let mut writer = ForceWriter::new_simple(path, title)?;
    writer.write_frame_simple(step, time, forces)?;
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_force_writer_creation() {
        let temp_file = "/tmp/test_force.trf";
        let writer = ForceWriter::new_simple(temp_file, "Test forces");
        assert!(writer.is_ok());
    }

    #[test]
    fn test_write_force_frame() {
        let temp_file = "/tmp/test_force2.trf";
        let mut writer = ForceWriter::new_simple(temp_file, "Test").unwrap();

        let forces = vec![
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(-1.0, -2.0, -3.0),
            Vec3::new(0.5, 0.5, 0.5),
        ];

        let result = writer.write_frame_simple(0, 0.0, &forces);
        assert!(result.is_ok());
        assert_eq!(writer.frame_count(), 1);
    }

    #[test]
    fn test_write_force_frame_with_constraints() {
        let temp_file = "/tmp/test_force_cons.trf";
        let forces = vec![
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(-1.0, -2.0, -3.0),
        ];
        let constraint_forces = vec![
            Vec3::new(0.1, 0.2, 0.3),
            Vec3::new(-0.1, -0.2, -0.3),
        ];

        let mut writer = ForceWriter::new(temp_file, "Test with constraints", true).unwrap();
        let result = writer.write_frame(0, 0.0, &forces, Some(&constraint_forces));
        assert!(result.is_ok());
        assert_eq!(writer.frame_count(), 1);
    }

    #[test]
    fn test_write_many_atoms() {
        // Test the comment insertion every 10 atoms
        let temp_file = "/tmp/test_force_many.trf";
        let forces = vec![Vec3::new(1.0, 2.0, 3.0); 25];

        let mut writer = ForceWriter::new_simple(temp_file, "Test many atoms").unwrap();
        let result = writer.write_frame_simple(1, 0.002, &forces);
        assert!(result.is_ok());
    }
}
