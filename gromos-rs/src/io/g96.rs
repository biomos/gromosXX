//! GROMOS96 (.g96) file format reader/writer
//!
//! GROMOS96 format is a block-based format with keywords like TITLE, POSITION, VELOCITY, BOX

use crate::math::Vec3;
use crate::configuration::Configuration;
use crate::topology::Topology;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;

/// Write a GROMOS96 format file
pub struct G96Writer {
    pub writer: BufWriter<File>,
}

impl G96Writer {
    /// Create a new G96 writer
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let file = File::create(path)
            .map_err(|e| format!("Cannot create .g96 file: {}", e))?;
        Ok(G96Writer {
            writer: BufWriter::new(file),
        })
    }

    /// Write title block
    pub fn write_title(&mut self, title: &str) -> Result<(), String> {
        writeln!(self.writer, "TITLE")
            .map_err(|e| format!("Write error: {}", e))?;
        writeln!(self.writer, "{}", title)
            .map_err(|e| format!("Write error: {}", e))?;
        writeln!(self.writer, "END")
            .map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write position block
    ///
    /// GROMOS96 position format:
    /// ```text
    /// POSITION
    /// %5d %-5s %-5s%7d%15.9f%15.9f%15.9f
    /// (resnum, resname, atomname, atomnum, x, y, z)
    /// END
    /// ```
    pub fn write_positions(
        &mut self,
        positions: &[Vec3],
        topology: Option<&Topology>,
    ) -> Result<(), String> {
        writeln!(self.writer, "POSITION")
            .map_err(|e| format!("Write error: {}", e))?;

        for (i, pos) in positions.iter().enumerate() {
            // Default values if no topology provided
            let resnum = 1;
            let resname = "UNK";
            let atomname = format!("AT{}", i + 1);
            let atomnum = i + 1;

            // Convert f32 to f64 for formatting
            let x = pos.x as f64;
            let y = pos.y as f64;
            let z = pos.z as f64;

            writeln!(
                self.writer,
                "{:5} {:5} {:5}{:7}{:15.9}{:15.9}{:15.9}",
                resnum, resname, atomname, atomnum, x, y, z
            )
            .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END")
            .map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write velocity block (optional)
    pub fn write_velocities(&mut self, velocities: &[Vec3]) -> Result<(), String> {
        writeln!(self.writer, "VELOCITY")
            .map_err(|e| format!("Write error: {}", e))?;

        for vel in velocities {
            let vx = vel.x as f64;
            let vy = vel.y as f64;
            let vz = vel.z as f64;

            writeln!(
                self.writer,
                "{:15.9}{:15.9}{:15.9}",
                vx, vy, vz
            )
            .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END")
            .map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write box dimensions (for periodic boundary conditions)
    pub fn write_box(&mut self, box_dims: Vec3) -> Result<(), String> {
        writeln!(self.writer, "BOX")
            .map_err(|e| format!("Write error: {}", e))?;

        let x = box_dims.x as f64;
        let y = box_dims.y as f64;
        let z = box_dims.z as f64;

        writeln!(
            self.writer,
            "{:15.9}{:15.9}{:15.9}",
            x, y, z
        )
        .map_err(|e| format!("Write error: {}", e))?;

        writeln!(self.writer, "END")
            .map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Flush and close the writer
    pub fn close(mut self) -> Result<(), String> {
        self.writer.flush()
            .map_err(|e| format!("Flush error: {}", e))
    }
}

/// Write a configuration to a GROMOS96 file
pub fn write_g96<P: AsRef<Path>>(
    path: P,
    title: &str,
    positions: &[Vec3],
    velocities: Option<&[Vec3]>,
    box_dims: Option<Vec3>,
    topology: Option<&Topology>,
) -> Result<(), String> {
    let mut writer = G96Writer::new(path)?;

    writer.write_title(title)?;
    writer.write_positions(positions, topology)?;

    if let Some(vels) = velocities {
        writer.write_velocities(vels)?;
    }

    if let Some(box_size) = box_dims {
        writer.write_box(box_size)?;
    }

    writer.close()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_g96_write() {
        let positions = vec![
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(4.0, 5.0, 6.0),
        ];

        let result = write_g96(
            "/tmp/test.g96",
            "Test structure",
            &positions,
            None,
            Some(Vec3::new(10.0, 10.0, 10.0)),
            None,
        );

        assert!(result.is_ok());
    }
}
