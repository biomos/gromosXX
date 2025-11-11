//! GROMOS energy file writer (.tre)
//!
//! Writes energy time series with ENE (energy) and ENA (energy analysis) blocks.
//!
//! # Format
//! ```text
//! TITLE
//!   Energy trajectory
//! END
//! ENERTRJ
//! # block 0 = time
//! # block 1 = kinetic energy
//! # block 2 = potential energy
//! # block 3 = total energy
//! # block 4 = temperature
//! # ...
//!      0.000    123.456    -567.890    -444.434    300.00 ...
//!      0.002    124.567    -568.901    -444.334    301.00 ...
//! END
//! ```
//!
//! # Energy Blocks (GROMOS convention)
//! The GROMOS .tre format uses specific block numbers for different energy components:
//!
//! ## ENE blocks (standard energies):
//! - Block 0: Time (ps)
//! - Block 1: Kinetic energy (kJ/mol)
//! - Block 2: Potential energy (kJ/mol)
//! - Block 3: Total energy (kJ/mol)
//! - Block 4: Temperature (K)
//! - Block 5: Volume (nm³)
//! - Block 6: Pressure (bar/GPa)
//!
//! ## Bonded energies:
//! - Block 7: Bond energy
//! - Block 8: Angle energy
//! - Block 9: Improper dihedral energy
//! - Block 10: Proper dihedral energy
//!
//! ## Nonbonded energies:
//! - Block 11: Lennard-Jones energy
//! - Block 12: Electrostatic energy (real space)
//! - Block 13: Electrostatic energy (reciprocal space, PME)
//! - Block 14: Self-energy correction
//!
//! ## Constraint energies:
//! - Block 15: SHAKE energy
//! - Block 16: Distance restraint energy
//!
//! ## Special energies:
//! - Block 17: Kinetic energy (translational)
//! - Block 18: Kinetic energy (rotational)
//! - Block 19: Kinetic energy (internal)

use crate::io::IoError;
use std::fs::File;
use std::io::{Write, BufWriter, BufReader, BufRead};
use std::path::Path;

/// Energy component identifiers (GROMOS convention)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyBlock {
    Time = 0,
    Kinetic = 1,
    Potential = 2,
    Total = 3,
    Temperature = 4,
    Volume = 5,
    Pressure = 6,

    // Bonded
    Bond = 7,
    Angle = 8,
    ImproperDihedral = 9,
    ProperDihedral = 10,

    // Nonbonded
    LennardJones = 11,
    ElectrostaticReal = 12,
    ElectrostaticReciprocal = 13,
    ElectrostaticSelf = 14,

    // Constraints
    Shake = 15,
    DistanceRestraint = 16,

    // Kinetic components
    KineticTrans = 17,
    KineticRot = 18,
    KineticInt = 19,
}

/// Energy data for a single timestep
#[derive(Debug, Clone)]
pub struct EnergyFrame {
    pub time: f64,
    pub kinetic: f64,
    pub potential: f64,
    pub total: f64,
    pub temperature: f64,
    pub volume: f64,
    pub pressure: f64,

    // Bonded energies
    pub bond: f64,
    pub angle: f64,
    pub improper: f64,
    pub dihedral: f64,

    // Nonbonded energies
    pub lj: f64,
    pub coul_real: f64,
    pub coul_recip: f64,
    pub coul_self: f64,

    // Constraints
    pub shake: f64,
    pub restraint: f64,

    // Additional components (optional)
    pub extra: Vec<f64>,
}

impl Default for EnergyFrame {
    fn default() -> Self {
        Self {
            time: 0.0,
            kinetic: 0.0,
            potential: 0.0,
            total: 0.0,
            temperature: 0.0,
            volume: 0.0,
            pressure: 0.0,
            bond: 0.0,
            angle: 0.0,
            improper: 0.0,
            dihedral: 0.0,
            lj: 0.0,
            coul_real: 0.0,
            coul_recip: 0.0,
            coul_self: 0.0,
            shake: 0.0,
            restraint: 0.0,
            extra: Vec::new(),
        }
    }
}

impl EnergyFrame {
    /// Create energy frame from components
    pub fn new(time: f64, kinetic: f64, potential: f64, temperature: f64) -> Self {
        Self {
            time,
            kinetic,
            potential,
            total: kinetic + potential,
            temperature,
            ..Default::default()
        }
    }

    /// Update total energy (sum of kinetic and potential)
    pub fn update_total(&mut self) {
        self.total = self.kinetic + self.potential;
    }

    /// Update potential energy (sum of all components)
    pub fn update_potential(&mut self) {
        self.potential = self.bond + self.angle + self.improper + self.dihedral
            + self.lj + self.coul_real + self.coul_recip + self.coul_self
            + self.shake + self.restraint;
    }
}

/// GROMOS energy file writer
pub struct EnergyWriter {
    writer: BufWriter<File>,
    frame_count: usize,
    write_header: bool,
}

impl EnergyWriter {
    /// Create a new energy writer
    pub fn new<P: AsRef<Path>>(path: P, title: &str) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write TITLE block
        writeln!(writer, "TITLE")?;
        writeln!(writer, "  {}", title)?;
        writeln!(writer, "END")?;

        // Start ENERTRJ block
        writeln!(writer, "ENERTRJ")?;

        Ok(Self {
            writer,
            frame_count: 0,
            write_header: true,
        })
    }

    /// Write energy frame in ENE format (standard GROMOS)
    pub fn write_frame(&mut self, frame: &EnergyFrame) -> Result<(), IoError> {
        // Write header comment on first frame
        if self.write_header {
            writeln!(self.writer, "# Block definitions (ENE format):")?;
            writeln!(self.writer, "# {:>3} = Time (ps)", 0)?;
            writeln!(self.writer, "# {:>3} = Kinetic energy (kJ/mol)", 1)?;
            writeln!(self.writer, "# {:>3} = Potential energy (kJ/mol)", 2)?;
            writeln!(self.writer, "# {:>3} = Total energy (kJ/mol)", 3)?;
            writeln!(self.writer, "# {:>3} = Temperature (K)", 4)?;
            writeln!(self.writer, "# {:>3} = Volume (nm³)", 5)?;
            writeln!(self.writer, "# {:>3} = Pressure (bar)", 6)?;
            writeln!(self.writer, "# {:>3} = Bond energy", 7)?;
            writeln!(self.writer, "# {:>3} = Angle energy", 8)?;
            writeln!(self.writer, "# {:>3} = Improper dihedral", 9)?;
            writeln!(self.writer, "# {:>3} = Proper dihedral", 10)?;
            writeln!(self.writer, "# {:>3} = Lennard-Jones", 11)?;
            writeln!(self.writer, "# {:>3} = Coulomb (real)", 12)?;
            writeln!(self.writer, "# {:>3} = Coulomb (reciprocal)", 13)?;
            writeln!(self.writer, "# {:>3} = Coulomb (self)", 14)?;
            writeln!(self.writer, "# {:>3} = SHAKE constraint", 15)?;
            writeln!(self.writer, "# {:>3} = Distance restraint", 16)?;
            self.write_header = false;
        }

        // Write energy values (20 columns, GROMOS standard)
        write!(
            self.writer,
            "{:15.6} {:15.6} {:15.6} {:15.6} {:15.6}",
            frame.time,
            frame.kinetic,
            frame.potential,
            frame.total,
            frame.temperature
        )?;

        write!(
            self.writer,
            " {:15.6} {:15.6}",
            frame.volume,
            frame.pressure
        )?;

        // Bonded energies
        write!(
            self.writer,
            " {:15.6} {:15.6} {:15.6} {:15.6}",
            frame.bond,
            frame.angle,
            frame.improper,
            frame.dihedral
        )?;

        // Nonbonded energies
        write!(
            self.writer,
            " {:15.6} {:15.6} {:15.6} {:15.6}",
            frame.lj,
            frame.coul_real,
            frame.coul_recip,
            frame.coul_self
        )?;

        // Constraints
        write!(
            self.writer,
            " {:15.6} {:15.6}",
            frame.shake,
            frame.restraint
        )?;

        writeln!(self.writer)?;

        self.frame_count += 1;
        Ok(())
    }

    /// Write energy frame in ENA format (energy analysis, more detail)
    ///
    /// ENA format includes additional breakdown of energy terms
    /// by energy group (useful for free energy calculations)
    pub fn write_frame_ena(
        &mut self,
        frame: &EnergyFrame,
        energy_groups: &[Vec<f64>],
    ) -> Result<(), IoError> {
        // Write standard frame first
        self.write_frame(frame)?;

        // Write additional energy group data
        writeln!(self.writer, "# ENA - Energy group analysis:")?;
        for (i, group) in energy_groups.iter().enumerate() {
            write!(self.writer, "# Group {:3}:", i)?;
            for value in group {
                write!(self.writer, " {:15.6}", value)?;
            }
            writeln!(self.writer)?;
        }

        Ok(())
    }

    /// Finalize and close the energy file
    pub fn finalize(&mut self) -> Result<(), IoError> {
        writeln!(self.writer, "END")?;
        self.writer.flush()?;
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

/// Write a single energy frame to a file (convenience function)
pub fn write_energy_frame<P: AsRef<Path>>(
    path: P,
    frame: &EnergyFrame,
    title: &str,
) -> Result<(), IoError> {
    let mut writer = EnergyWriter::new(path, title)?;
    writer.write_frame(frame)?;
    writer.finalize()?;
    Ok(())
}

/// GROMOS energy file reader
pub struct EnergyReader {
    reader: BufReader<File>,
    title: String,
    frames_read: usize,
}

impl EnergyReader {
    /// Open an energy trajectory file for reading
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, IoError> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Read TITLE block
        let title = Self::read_title_block(&mut reader)?;

        Ok(Self {
            reader,
            title,
            frames_read: 0,
        })
    }

    /// Get the title from the energy file
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Get the number of frames read so far
    pub fn frames_read(&self) -> usize {
        self.frames_read
    }

    /// Read a single energy frame
    ///
    /// Returns Ok(Some(frame)) if a frame was read, Ok(None) if end of file reached
    pub fn read_frame(&mut self) -> Result<Option<EnergyFrame>, IoError> {
        let mut line = String::new();

        loop {
            line.clear();
            let bytes_read = self.reader.read_line(&mut line)?;

            if bytes_read == 0 {
                // End of file
                return Ok(None);
            }

            let trimmed = line.trim();

            // Skip empty lines and comments
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            // Check for END block
            if trimmed == "END" {
                return Ok(None);
            }

            // Check for ENERTRJ block
            if trimmed == "ENERTRJ" {
                continue;
            }

            // Parse energy values
            let values: Result<Vec<f64>, _> = trimmed
                .split_whitespace()
                .map(|s| s.parse::<f64>())
                .collect();

            match values {
                Ok(vals) if vals.len() >= 17 => {
                    // Successfully parsed energy frame
                    let frame = EnergyFrame {
                        time: vals[0],
                        kinetic: vals[1],
                        potential: vals[2],
                        total: vals[3],
                        temperature: vals[4],
                        volume: vals[5],
                        pressure: vals[6],
                        bond: vals[7],
                        angle: vals[8],
                        improper: vals[9],
                        dihedral: vals[10],
                        lj: vals[11],
                        coul_real: vals[12],
                        coul_recip: vals[13],
                        coul_self: vals[14],
                        shake: vals[15],
                        restraint: vals[16],
                        extra: if vals.len() > 17 {
                            vals[17..].to_vec()
                        } else {
                            Vec::new()
                        },
                    };

                    self.frames_read += 1;
                    return Ok(Some(frame));
                }
                Ok(vals) if vals.is_empty() => {
                    // Empty line after splitting, continue
                    continue;
                }
                Ok(vals) => {
                    return Err(IoError::ParseError(format!(
                        "Incomplete energy frame: expected at least 17 values, got {}",
                        vals.len()
                    )));
                }
                Err(e) => {
                    return Err(IoError::ParseError(format!(
                        "Failed to parse energy values: {}",
                        e
                    )));
                }
            }
        }
    }

    /// Read all energy frames from the file
    pub fn read_all_frames(&mut self) -> Result<Vec<EnergyFrame>, IoError> {
        let mut frames = Vec::new();
        while let Some(frame) = self.read_frame()? {
            frames.push(frame);
        }
        Ok(frames)
    }

    /// Read TITLE block
    fn read_title_block(reader: &mut BufReader<File>) -> Result<String, IoError> {
        let mut line = String::new();
        let mut title_lines = Vec::new();
        let mut in_title = false;

        loop {
            line.clear();
            let bytes_read = reader.read_line(&mut line)?;

            if bytes_read == 0 {
                return Err(IoError::ParseError("Unexpected EOF while reading TITLE".to_string()));
            }

            let trimmed = line.trim();

            if trimmed == "TITLE" {
                in_title = true;
                continue;
            }

            if trimmed == "END" {
                if in_title {
                    return Ok(title_lines.join(" ").trim().to_string());
                }
                continue;
            }

            if in_title {
                title_lines.push(trimmed.to_string());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_energy_frame_creation() {
        let frame = EnergyFrame::new(0.0, 100.0, -200.0, 300.0);
        assert_eq!(frame.time, 0.0);
        assert_eq!(frame.kinetic, 100.0);
        assert_eq!(frame.potential, -200.0);
        assert_eq!(frame.total, -100.0);
        assert_eq!(frame.temperature, 300.0);
    }

    #[test]
    fn test_energy_writer_creation() {
        let temp_file = "/tmp/test_energy.tre";
        let writer = EnergyWriter::new(temp_file, "Test energy");
        assert!(writer.is_ok());
    }

    #[test]
    fn test_write_energy_frame() {
        let temp_file = "/tmp/test_energy2.tre";
        let mut writer = EnergyWriter::new(temp_file, "Test").unwrap();

        let mut frame = EnergyFrame::new(0.0, 150.0, -300.0, 298.15);
        frame.bond = -50.0;
        frame.angle = -30.0;
        frame.lj = -150.0;
        frame.coul_real = -70.0;
        frame.update_potential();

        let result = writer.write_frame(&frame);
        assert!(result.is_ok());
        assert_eq!(writer.frame_count(), 1);

        let result = writer.finalize();
        assert!(result.is_ok());
    }
}
