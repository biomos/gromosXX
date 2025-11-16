/// GaMD statistics output writer
///
/// Writes GaMD statistics to file for parameter estimation and analysis.
/// Format follows GROMOS convention for compatibility with analysis tools.
///
/// Output format:
/// ```text
/// # GaMD Statistics
/// # Step  V_dih       V_tot       V_max_dih  V_min_dih  V_avg_dih  σ_dih   V_max_tot  V_min_tot  V_avg_tot  σ_tot    k_dih    k_tot      E_dih    E_tot
/// 1000    -125.3      -15234.5    -100.2     -145.8     -123.4     12.3    -15100.0   -15350.0   -15220.0   58.4     0.0      0.001234   0.0      -15100.0
/// ```

use std::fs::File;
use std::io::{Write, BufWriter};
use crate::gamd::{GamdParameters, GamdStatistics};

/// Writer for GaMD statistics file
pub struct GamdStatsWriter {
    writer: BufWriter<File>,
    write_interval: usize,
}

impl GamdStatsWriter {
    /// Create a new GaMD statistics writer
    pub fn new(path: &str, title: &str) -> Result<Self, String> {
        let file = File::create(path)
            .map_err(|e| format!("Failed to create GaMD stats file {}: {}", path, e))?;

        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(writer, "# {}", title)
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# GaMD Statistics File")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "#")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Columns:")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Step      : MD step number")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_dih     : Current dihedral potential energy (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_tot     : Current total potential energy (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_max_dih : Maximum dihedral energy observed (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_min_dih : Minimum dihedral energy observed (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_avg_dih : Average dihedral energy (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# σ_dih     : Standard deviation of dihedral energy (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_max_tot : Maximum total energy observed (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_min_tot : Minimum total energy observed (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# V_avg_tot : Average total energy (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# σ_tot     : Standard deviation of total energy (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# k_dih     : GaMD boost coefficient for dihedral (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# k_tot     : GaMD boost coefficient for total (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# E_dih     : GaMD threshold energy for dihedral (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# E_tot     : GaMD threshold energy for total (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "#")
            .map_err(|e| format!("Failed to write header: {}", e))?;

        // Write column headers
        writeln!(writer, "# {:>8} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}",
            "Step", "V_dih", "V_tot", "V_max_dih", "V_min_dih", "V_avg_dih", "σ_dih",
            "V_max_tot", "V_min_tot", "V_avg_tot", "σ_tot", "k_dih", "k_tot", "E_dih", "E_tot")
            .map_err(|e| format!("Failed to write column headers: {}", e))?;

        Ok(GamdStatsWriter {
            writer,
            write_interval: 100, // Write every 100 steps by default
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
        v_dih: f64,
        v_tot: f64,
        gamd: &GamdParameters,
    ) -> Result<(), String> {
        // Only write at specified intervals
        if step % self.write_interval != 0 {
            return Ok(());
        }

        let dih_stats = &gamd.dih_stats;
        let tot_stats = &gamd.tot_stats;

        writeln!(
            self.writer,
            "  {:8} {:12.4} {:12.4} {:12.4} {:12.4} {:12.4} {:12.4} {:12.4} {:12.4} {:12.4} {:12.4} {:12.6} {:12.6} {:12.4} {:12.4}",
            step,
            v_dih,
            v_tot,
            dih_stats.v_max,
            dih_stats.v_min,
            dih_stats.v_mean,
            dih_stats.sigma_v,
            tot_stats.v_max,
            tot_stats.v_min,
            tot_stats.v_mean,
            tot_stats.sigma_v,
            gamd.k_dih,
            gamd.k_tot,
            gamd.e_dih,
            gamd.e_tot,
        )
        .map_err(|e| format!("Failed to write GaMD stats frame: {}", e))?;

        Ok(())
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<(), String> {
        self.writer.flush()
            .map_err(|e| format!("Failed to flush GaMD stats writer: {}", e))
    }

    /// Finalize and close the file
    pub fn finalize(&mut self) -> Result<(), String> {
        self.flush()?;
        writeln!(self.writer, "# END")
            .map_err(|e| format!("Failed to write END marker: {}", e))?;
        self.flush()
    }
}

/// Writer for GaMD boost potential trajectory
///
/// Writes the boost potential at each step for reweighting analysis.
/// This file is used by gamd_ana tool to compute reweighted free energies.
pub struct GamdBoostWriter {
    writer: BufWriter<File>,
}

impl GamdBoostWriter {
    /// Create a new GaMD boost writer
    pub fn new(path: &str, title: &str) -> Result<Self, String> {
        let file = File::create(path)
            .map_err(|e| format!("Failed to create GaMD boost file {}: {}", path, e))?;

        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(writer, "# {}", title)
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# GaMD Boost Potential Trajectory")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "#")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# This file contains the boost potential ΔV at each MD step")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# for reweighting analysis to recover unbiased free energies.")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "#")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Columns:")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Step  : MD step number")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# Time  : Simulation time (ps)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# ΔV    : Total boost potential (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# ΔV_dih: Dihedral boost potential (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# ΔV_tot: Total potential boost (kJ/mol)")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "#")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        writeln!(writer, "# {:>8} {:>12} {:>12} {:>12} {:>12}",
            "Step", "Time", "ΔV", "ΔV_dih", "ΔV_tot")
            .map_err(|e| format!("Failed to write column headers: {}", e))?;

        Ok(GamdBoostWriter { writer })
    }

    /// Write a boost frame
    pub fn write_frame(
        &mut self,
        step: usize,
        time: f64,
        boost_total: f64,
        boost_dih: f64,
        boost_pot: f64,
    ) -> Result<(), String> {
        writeln!(
            self.writer,
            "  {:8} {:12.6} {:12.6} {:12.6} {:12.6}",
            step, time, boost_total, boost_dih, boost_pot
        )
        .map_err(|e| format!("Failed to write boost frame: {}", e))
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<(), String> {
        self.writer.flush()
            .map_err(|e| format!("Failed to flush boost writer: {}", e))
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
    use crate::gamd::{SearchMode, BoostForm, ThresholdType};
    use std::fs;

    #[test]
    fn test_gamd_stats_writer() {
        let path = "/tmp/test_gamd_stats.dat";

        // Create writer
        let mut writer = GamdStatsWriter::new(path, "Test GaMD Stats").unwrap();
        writer.set_write_interval(1);

        // Create test parameters
        let mut gamd = GamdParameters::new(
            SearchMode::CmdSearch,
            BoostForm::TotalBoost,
            ThresholdType::LowerBound,
            6.0,
            6.0,
        );

        // Write some frames
        gamd.update_statistics(0.0, -1000.0);
        writer.write_frame(0, 0.0, -1000.0, &gamd).unwrap();

        gamd.update_statistics(0.0, -1100.0);
        writer.write_frame(1, 0.0, -1100.0, &gamd).unwrap();

        writer.finalize().unwrap();

        // Check file exists and has content
        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("GaMD Statistics"));
        assert!(content.contains("-1000.0"));
        assert!(content.contains("END"));

        // Cleanup
        fs::remove_file(path).ok();
    }

    #[test]
    fn test_gamd_boost_writer() {
        let path = "/tmp/test_gamd_boost.dat";

        let mut writer = GamdBoostWriter::new(path, "Test Boost").unwrap();

        writer.write_frame(0, 0.0, 5.0, 0.0, 5.0).unwrap();
        writer.write_frame(1, 0.002, 6.2, 0.0, 6.2).unwrap();

        writer.finalize().unwrap();

        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("Boost Potential"));
        assert!(content.contains("5.0"));
        assert!(content.contains("6.2"));

        fs::remove_file(path).ok();
    }
}
