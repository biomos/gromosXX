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

use crate::configuration::Configuration;
use crate::math::Vec3;
use crate::io::IoError;
use std::fs::File;
use std::io::{Write, BufWriter, BufReader, BufRead, Seek};
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

/// A single frame from a trajectory file
#[derive(Debug, Clone)]
pub struct TrajectoryFrame {
    /// Simulation step number
    pub step: usize,
    /// Simulation time (ps)
    pub time: f64,
    /// Atomic positions (nm)
    pub positions: Vec<Vec3>,
    /// Atomic velocities (nm/ps), if present
    pub velocities: Option<Vec<Vec3>>,
    /// Atomic forces (kJ/(mol·nm)), if present
    pub forces: Option<Vec<Vec3>>,
    /// Lattice shifts (integer vectors), if present
    pub lattice_shifts: Option<Vec<(i32, i32, i32)>>,
    /// Box dimensions (nm)
    pub box_dims: Vec3,
}

/// GROMOS trajectory reader (.trc/.trj files)
///
/// Reads trajectory files frame by frame, parsing POSITIONRED, VELOCITYRED, FORCERED,
/// LATTICESHIFTS, and GENBOX blocks.
pub struct TrajectoryReader {
    reader: BufReader<File>,
    title: String,
    frames_read: usize,
    buffer: String,
}

impl TrajectoryReader {
    /// Open a trajectory file for reading
    ///
    /// # Arguments
    /// * `path` - Path to the .trc or .trj file
    ///
    /// # Returns
    /// * `Ok(TrajectoryReader)` - Successfully opened the file
    /// * `Err(IoError)` - Failed to open or parse the header
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, IoError> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        let mut buffer = String::new();

        // Read TITLE block
        let title = Self::read_title_block(&mut reader, &mut buffer)?;

        Ok(Self {
            reader,
            title,
            frames_read: 0,
            buffer,
        })
    }

    /// Get the title from the trajectory file
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Get the number of frames read so far
    pub fn frames_read(&self) -> usize {
        self.frames_read
    }

    /// Read the next frame from the trajectory
    ///
    /// # Returns
    /// * `Ok(Some(TrajectoryFrame))` - Successfully read a frame
    /// * `Ok(None)` - End of file reached
    /// * `Err(IoError)` - Error reading the frame
    pub fn read_frame(&mut self) -> Result<Option<TrajectoryFrame>, IoError> {
        // Try to read TIMESTEP block
        let (step, time) = match Self::read_timestep_block(&mut self.reader, &mut self.buffer)? {
            Some(st) => st,
            None => return Ok(None), // End of file
        };

        // Read POSITIONRED block (required)
        let positions = Self::read_position_block(&mut self.reader, &mut self.buffer)?;

        // Try to read optional blocks
        let velocities = Self::try_read_velocity_block(&mut self.reader, &mut self.buffer)?;
        let forces = Self::try_read_force_block(&mut self.reader, &mut self.buffer)?;
        let lattice_shifts = Self::try_read_lattice_shifts(&mut self.reader, &mut self.buffer)?;

        // Read GENBOX block (required)
        let box_dims = Self::read_genbox_block(&mut self.reader, &mut self.buffer)?;

        self.frames_read += 1;

        Ok(Some(TrajectoryFrame {
            step,
            time,
            positions,
            velocities,
            forces,
            lattice_shifts,
            box_dims,
        }))
    }

    /// Read all frames from the trajectory
    ///
    /// # Returns
    /// * `Ok(Vec<TrajectoryFrame>)` - All frames successfully read
    /// * `Err(IoError)` - Error reading frames
    pub fn read_all_frames(&mut self) -> Result<Vec<TrajectoryFrame>, IoError> {
        let mut frames = Vec::new();
        while let Some(frame) = self.read_frame()? {
            frames.push(frame);
        }
        Ok(frames)
    }

    // Helper functions for reading blocks

    fn read_title_block(reader: &mut BufReader<File>, buffer: &mut String) -> Result<String, IoError> {
        buffer.clear();
        reader.read_line(buffer)?;
        if !buffer.trim().starts_with("TITLE") {
            return Err(IoError::FormatError(
                "Expected TITLE block at start of trajectory".to_string(),
            ));
        }

        buffer.clear();
        reader.read_line(buffer)?;
        let title = buffer.trim().to_string();

        buffer.clear();
        reader.read_line(buffer)?;
        if !buffer.trim().starts_with("END") {
            return Err(IoError::FormatError("Expected END after TITLE".to_string()));
        }

        Ok(title)
    }

    fn read_timestep_block(
        reader: &mut BufReader<File>,
        buffer: &mut String,
    ) -> Result<Option<(usize, f64)>, IoError> {
        buffer.clear();
        let bytes_read = reader.read_line(buffer)?;
        if bytes_read == 0 {
            return Ok(None); // EOF
        }

        let line = buffer.trim();
        if line.is_empty() {
            return Self::read_timestep_block(reader, buffer); // Skip empty lines
        }

        if !line.starts_with("TIMESTEP") {
            return Err(IoError::FormatError(format!(
                "Expected TIMESTEP block, got: {}",
                line
            )));
        }

        // Read step and time
        buffer.clear();
        reader.read_line(buffer)?;
        let parts: Vec<&str> = buffer.trim().split_whitespace().collect();
        if parts.len() < 2 {
            return Err(IoError::FormatError(
                "TIMESTEP data should have step and time".to_string(),
            ));
        }

        let step = parts[0]
            .parse::<usize>()
            .map_err(|e| IoError::FormatError(format!("Invalid step number: {}", e)))?;
        let time = parts[1]
            .parse::<f64>()
            .map_err(|e| IoError::FormatError(format!("Invalid time: {}", e)))?;

        // Read END
        buffer.clear();
        reader.read_line(buffer)?;
        if !buffer.trim().starts_with("END") {
            return Err(IoError::FormatError("Expected END after TIMESTEP".to_string()));
        }

        Ok(Some((step, time)))
    }

    fn read_position_block(
        reader: &mut BufReader<File>,
        buffer: &mut String,
    ) -> Result<Vec<Vec3>, IoError> {
        buffer.clear();
        reader.read_line(buffer)?;
        if !buffer.trim().starts_with("POSITIONRED") {
            return Err(IoError::FormatError(
                "Expected POSITIONRED block".to_string(),
            ));
        }

        let mut positions = Vec::new();
        loop {
            buffer.clear();
            reader.read_line(buffer)?;
            let line = buffer.trim();

            if line.starts_with("END") {
                break;
            }

            if line.starts_with('#') || line.is_empty() {
                continue; // Skip comments and empty lines
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 7 {
                // Format: atom_num RES ATOM serial x y z
                let x = parts[4]
                    .parse::<f32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid x coordinate: {}", e)))?;
                let y = parts[5]
                    .parse::<f32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid y coordinate: {}", e)))?;
                let z = parts[6]
                    .parse::<f32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid z coordinate: {}", e)))?;
                positions.push(Vec3::new(x, y, z));
            }
        }

        Ok(positions)
    }

    fn try_read_velocity_block(
        reader: &mut BufReader<File>,
        buffer: &mut String,
    ) -> Result<Option<Vec<Vec3>>, IoError> {
        buffer.clear();
        let position = reader.stream_position()?;
        reader.read_line(buffer)?;

        if !buffer.trim().starts_with("VELOCITYRED") {
            // Not a velocity block, rewind
            reader.seek(std::io::SeekFrom::Start(position))?;
            return Ok(None);
        }

        let mut velocities = Vec::new();
        loop {
            buffer.clear();
            reader.read_line(buffer)?;
            let line = buffer.trim();

            if line.starts_with("END") {
                break;
            }

            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 7 {
                let vx = parts[4].parse::<f32>().map_err(|e| {
                    IoError::FormatError(format!("Invalid vx velocity: {}", e))
                })?;
                let vy = parts[5].parse::<f32>().map_err(|e| {
                    IoError::FormatError(format!("Invalid vy velocity: {}", e))
                })?;
                let vz = parts[6].parse::<f32>().map_err(|e| {
                    IoError::FormatError(format!("Invalid vz velocity: {}", e))
                })?;
                velocities.push(Vec3::new(vx, vy, vz));
            }
        }

        Ok(Some(velocities))
    }

    fn try_read_force_block(
        reader: &mut BufReader<File>,
        buffer: &mut String,
    ) -> Result<Option<Vec<Vec3>>, IoError> {
        buffer.clear();
        let position = reader.stream_position()?;
        reader.read_line(buffer)?;

        if !buffer.trim().starts_with("FORCERED") {
            // Not a force block, rewind
            reader.seek(std::io::SeekFrom::Start(position))?;
            return Ok(None);
        }

        let mut forces = Vec::new();
        loop {
            buffer.clear();
            reader.read_line(buffer)?;
            let line = buffer.trim();

            if line.starts_with("END") {
                break;
            }

            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 7 {
                let fx = parts[4]
                    .parse::<f32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid fx force: {}", e)))?;
                let fy = parts[5]
                    .parse::<f32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid fy force: {}", e)))?;
                let fz = parts[6]
                    .parse::<f32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid fz force: {}", e)))?;
                forces.push(Vec3::new(fx, fy, fz));
            }
        }

        Ok(Some(forces))
    }

    fn try_read_lattice_shifts(
        reader: &mut BufReader<File>,
        buffer: &mut String,
    ) -> Result<Option<Vec<(i32, i32, i32)>>, IoError> {
        buffer.clear();
        let position = reader.stream_position()?;
        reader.read_line(buffer)?;

        if !buffer.trim().starts_with("LATTICESHIFTS") {
            // Not a lattice shifts block, rewind
            reader.seek(std::io::SeekFrom::Start(position))?;
            return Ok(None);
        }

        let mut shifts = Vec::new();
        loop {
            buffer.clear();
            reader.read_line(buffer)?;
            let line = buffer.trim();

            if line.starts_with("END") {
                break;
            }

            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                let sx = parts[1]
                    .parse::<i32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid shift x: {}", e)))?;
                let sy = parts[2]
                    .parse::<i32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid shift y: {}", e)))?;
                let sz = parts[3]
                    .parse::<i32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid shift z: {}", e)))?;
                shifts.push((sx, sy, sz));
            }
        }

        Ok(Some(shifts))
    }

    fn read_genbox_block(
        reader: &mut BufReader<File>,
        buffer: &mut String,
    ) -> Result<Vec3, IoError> {
        buffer.clear();
        reader.read_line(buffer)?;
        if !buffer.trim().starts_with("GENBOX") {
            return Err(IoError::FormatError("Expected GENBOX block".to_string()));
        }

        buffer.clear();
        reader.read_line(buffer)?;
        let parts: Vec<&str> = buffer.trim().split_whitespace().collect();
        if parts.len() < 3 {
            return Err(IoError::FormatError(
                "GENBOX should have 3 dimensions".to_string(),
            ));
        }

        let lx = parts[0]
            .parse::<f32>()
            .map_err(|e| IoError::FormatError(format!("Invalid box x: {}", e)))?;
        let ly = parts[1]
            .parse::<f32>()
            .map_err(|e| IoError::FormatError(format!("Invalid box y: {}", e)))?;
        let lz = parts[2]
            .parse::<f32>()
            .map_err(|e| IoError::FormatError(format!("Invalid box z: {}", e)))?;

        buffer.clear();
        reader.read_line(buffer)?;
        if !buffer.trim().starts_with("END") {
            return Err(IoError::FormatError("Expected END after GENBOX".to_string()));
        }

        Ok(Vec3::new(lx, ly, lz))
    }
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
