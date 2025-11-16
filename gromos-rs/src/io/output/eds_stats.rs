/// EDS statistics output writer
///
/// Writes EDS statistics to file for analysis.
/// Format follows GROMOS convention for compatibility with analysis tools.

use std::fs::File;
use std::io::{Write, BufWriter};
use crate::eds::EDSParameters;

/// Writer for EDS state visit statistics
pub struct EdsStatsWriter {
    writer: BufWriter<File>,
    write_interval: usize,
}

impl EdsStatsWriter {
    /// Create a new EDS statistics writer
    pub fn new(path: &str, title: &str) -> Result<Self, String> {
        let file = File::create(path)
            .map_err(|e| format!("Failed to create EDS stats file {}: {}", path, e))?;

        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(writer, "# {}", title)
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# EDS Statistics File")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "#")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Columns:")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Step     : MD step number")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# State_i  : Visit count for state i")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_i      : Energy of state i (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_R      : Reference energy (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "#")
            .map_err(|e| format!("Failed to write header: {}", e))?;

        Ok(EdsStatsWriter {
            writer,
            write_interval: 100,
        })
    }

    /// Set the write interval
    pub fn set_write_interval(&mut self, interval: usize) {
        self.write_interval = interval;
    }

    /// Write a statistics frame
    pub fn write_frame(
        &mut self,
        step: usize,
        eds: &EDSParameters,
    ) -> Result<(), String> {
        // Only write at specified intervals
        if step % self.write_interval != 0 {
            return Ok(());
        }

        // Write step and reference energy
        write!(self.writer, "  {:8} {:12.4}", step, eds.reference_energy)
            .map_err(|e| format!("Failed to write EDS stats: {}", e))?;

        // Write state energies
        for state in &eds.states {
            write!(self.writer, " {:12.4}", state.energy)
                .map_err(|e| format!("Failed to write state energy: {}", e))?;
        }

        // Write state visit counts
        for state in &eds.states {
            write!(self.writer, " {:8}", state.visit_count)
                .map_err(|e| format!("Failed to write visit count: {}", e))?;
        }

        writeln!(self.writer)
            .map_err(|e| format!("Failed to write newline: {}", e))?;

        Ok(())
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<(), String> {
        self.writer.flush()
            .map_err(|e| format!("Failed to flush EDS stats writer: {}", e))
    }

    /// Finalize and close the file
    pub fn finalize(&mut self) -> Result<(), String> {
        self.flush()?;
        writeln!(self.writer, "# END")
            .map_err(|e| format!("Failed to write END marker: {}", e))?;
        self.flush()
    }
}

/// Writer for EDS reference energy trajectory
///
/// Writes the reference energy V_R at each step for analysis.
pub struct EdsVrWriter {
    writer: BufWriter<File>,
}

impl EdsVrWriter {
    /// Create a new EDS V_R writer
    pub fn new(path: &str, title: &str) -> Result<Self, String> {
        let file = File::create(path)
            .map_err(|e| format!("Failed to create EDS V_R file {}: {}", path, e))?;

        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(writer, "# {}", title)
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# EDS Reference Energy Trajectory")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "#")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Contains the reference energy V_R at each MD step")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# for computing ensemble averages and free energy differences.")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "#")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Columns:")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Step : MD step number")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Time : Simulation time (ps)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_R  : Reference energy (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "#")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# {:>8} {:>12} {:>12}", "Step", "Time", "V_R")
            .map_err(|e| format!("Failed to write column headers: {}", e))?;

        Ok(EdsVrWriter { writer })
    }

    /// Write a V_R frame
    pub fn write_frame(
        &mut self,
        step: usize,
        time: f64,
        v_r: f64,
    ) -> Result<(), String> {
        writeln!(
            self.writer,
            "  {:8} {:12.6} {:12.6}",
            step, time, v_r
        )
        .map_err(|e| format!("Failed to write V_R frame: {}", e))
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<(), String> {
        self.writer.flush()
            .map_err(|e| format!("Failed to flush V_R writer: {}", e))
    }

    /// Finalize and close the file
    pub fn finalize(&mut self) -> Result<(), String> {
        self.flush()?;
        writeln!(self.writer, "# END")
            .map_err(|e| format!("Failed to write END marker: {}", e))?;
        self.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eds::{EDSForm, EDSState};
    use std::fs;

    #[test]
    fn test_eds_stats_writer() {
        let path = "/tmp/test_eds_stats.dat";

        let mut writer = EdsStatsWriter::new(path, "Test EDS Stats").unwrap();
        writer.set_write_interval(1);

        // Create test EDS parameters
        let mut eds = EDSParameters {
            form: EDSForm::SingleS,
            num_states: 2,
            s_values: vec![0.5],
            temperature: 300.0,
            states: vec![
                EDSState::new(0, 0.0, 10),
                EDSState::new(1, 10.0, 10),
            ],
            reference_energy: -100.0,
        };

        eds.states[0].energy = -95.0;
        eds.states[1].energy = -105.0;

        writer.write_frame(0, &eds).unwrap();
        writer.write_frame(1, &eds).unwrap();

        writer.finalize().unwrap();

        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("EDS Statistics"));
        assert!(content.contains("-100.0"));
        assert!(content.contains("END"));

        fs::remove_file(path).ok();
    }

    #[test]
    fn test_eds_vr_writer() {
        let path = "/tmp/test_eds_vr.dat";

        let mut writer = EdsVrWriter::new(path, "Test V_R").unwrap();

        writer.write_frame(0, 0.0, -100.5).unwrap();
        writer.write_frame(1, 0.002, -101.2).unwrap();

        writer.finalize().unwrap();

        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("Reference Energy"));
        assert!(content.contains("-100.5"));
        assert!(content.contains("-101.2"));

        fs::remove_file(path).ok();
    }
}
