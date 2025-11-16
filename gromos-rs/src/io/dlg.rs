//! GROMOS free energy output file writer (.dlg)
//!
//! Writes lambda derivatives (dH/dλ) for thermodynamic integration (TI) analysis.
//!
//! # Format
//! ```text
//! TITLE
//!   Free energy derivatives (dH/dλ)
//! END
//! FREEENERGYDERIVATIVES
//! # block 0 = time (ps)
//! # block 1 = lambda
//! # block 2 = dH/dλ total
//! # block 3 = dH/dλ bonds
//! # block 4 = dH/dλ angles
//! # block 5 = dH/dλ dihedrals
//! # block 6 = dH/dλ Lennard-Jones
//! # block 7 = dH/dλ Coulomb
//! # block 8 = dH/dλ restraints
//!      0.000  0.000    0.0000    0.0000    0.0000 ...
//!      0.002  0.000    1.2345    0.1234    0.2345 ...
//! END
//! ```
//!
//! # Usage in FEP/TI
//! During a free energy perturbation (FEP) or thermodynamic integration (TI) simulation:
//! 1. Run multiple simulations at different λ values (0.0 to 1.0)
//! 2. Record dH/dλ at each timestep using DlgWriter
//! 3. Analyze with GROMOS++ tools (ext_ti_ana, bar) to get ΔG
//!
//! ΔG = ∫₀¹ <dH/dλ>_λ dλ

use crate::io::IoError;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;

/// Lambda derivative data for a single timestep
///
/// Records ∂H/∂λ for different energy components during FEP/TI simulations
#[derive(Debug, Clone)]
pub struct LambdaDerivativeFrame {
    /// Simulation time (ps)
    pub time: f64,

    /// Lambda value (0.0 to 1.0)
    pub lambda: f64,

    /// Total dH/dλ (sum of all components)
    pub total: f64,

    /// dH/dλ for bonds
    pub bond: f64,

    /// dH/dλ for angles
    pub angle: f64,

    /// dH/dλ for improper dihedrals
    pub improper: f64,

    /// dH/dλ for proper dihedrals
    pub dihedral: f64,

    /// dH/dλ for Lennard-Jones interactions
    pub lj: f64,

    /// dH/dλ for Coulomb interactions
    pub coulomb: f64,

    /// dH/dλ for restraints
    pub restraint: f64,

    /// Additional lambda derivatives (optional)
    pub extra: Vec<f64>,
}

impl Default for LambdaDerivativeFrame {
    fn default() -> Self {
        Self {
            time: 0.0,
            lambda: 0.0,
            total: 0.0,
            bond: 0.0,
            angle: 0.0,
            improper: 0.0,
            dihedral: 0.0,
            lj: 0.0,
            coulomb: 0.0,
            restraint: 0.0,
            extra: Vec::new(),
        }
    }
}

impl LambdaDerivativeFrame {
    /// Create lambda derivative frame from time and lambda
    pub fn new(time: f64, lambda: f64) -> Self {
        Self {
            time,
            lambda,
            ..Default::default()
        }
    }

    /// Update total derivative (sum of all components)
    pub fn update_total(&mut self) {
        self.total = self.bond + self.angle + self.improper + self.dihedral
            + self.lj + self.coulomb + self.restraint;

        for &extra_val in &self.extra {
            self.total += extra_val;
        }
    }

    /// Set all bonded derivatives at once
    pub fn set_bonded(&mut self, bond: f64, angle: f64, improper: f64, dihedral: f64) {
        self.bond = bond;
        self.angle = angle;
        self.improper = improper;
        self.dihedral = dihedral;
    }

    /// Set all nonbonded derivatives at once
    pub fn set_nonbonded(&mut self, lj: f64, coulomb: f64) {
        self.lj = lj;
        self.coulomb = coulomb;
    }
}

/// GROMOS free energy derivative file writer (.dlg)
pub struct DlgWriter {
    writer: BufWriter<File>,
    frame_count: usize,
    write_header: bool,
}

impl DlgWriter {
    /// Create a new .dlg writer
    pub fn new<P: AsRef<Path>>(path: P, title: &str) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write TITLE block
        writeln!(writer, "TITLE")?;
        writeln!(writer, "  {}", title)?;
        writeln!(writer, "END")?;

        // Start FREEENERGYDERIVATIVES block
        writeln!(writer, "FREEENERGYDERIVATIVES")?;

        Ok(Self {
            writer,
            frame_count: 0,
            write_header: true,
        })
    }

    /// Write header comment (done once)
    fn write_header_comment(&mut self) -> Result<(), IoError> {
        writeln!(self.writer, "# Block 0: Time (ps)")?;
        writeln!(self.writer, "# Block 1: Lambda")?;
        writeln!(self.writer, "# Block 2: dH/dλ total")?;
        writeln!(self.writer, "# Block 3: dH/dλ bonds")?;
        writeln!(self.writer, "# Block 4: dH/dλ angles")?;
        writeln!(self.writer, "# Block 5: dH/dλ improper dihedrals")?;
        writeln!(self.writer, "# Block 6: dH/dλ proper dihedrals")?;
        writeln!(self.writer, "# Block 7: dH/dλ Lennard-Jones")?;
        writeln!(self.writer, "# Block 8: dH/dλ Coulomb")?;
        writeln!(self.writer, "# Block 9: dH/dλ restraints")?;
        Ok(())
    }

    /// Write a lambda derivative frame
    pub fn write_frame(&mut self, frame: &LambdaDerivativeFrame) -> Result<(), IoError> {
        // Write header on first frame
        if self.write_header {
            self.write_header_comment()?;
            self.write_header = false;
        }

        // Write frame data: time, lambda, total, bond, angle, improper, dihedral, lj, coulomb, restraint
        write!(
            self.writer,
            "{:12.6}{:12.6}{:16.6}{:16.6}{:16.6}{:16.6}{:16.6}{:16.6}{:16.6}{:16.6}",
            frame.time,
            frame.lambda,
            frame.total,
            frame.bond,
            frame.angle,
            frame.improper,
            frame.dihedral,
            frame.lj,
            frame.coulomb,
            frame.restraint,
        )?;

        // Write extra components if present
        for &extra in &frame.extra {
            write!(self.writer, "{:16.6}", extra)?;
        }

        writeln!(self.writer)?;
        self.frame_count += 1;

        Ok(())
    }

    /// Get number of frames written
    pub fn frame_count(&self) -> usize {
        self.frame_count
    }

    /// Flush and close the writer
    pub fn close(mut self) -> Result<(), IoError> {
        writeln!(self.writer, "END")?;
        self.writer.flush()?;
        Ok(())
    }
}

/// Convenience function to write lambda derivatives to a .dlg file
pub fn write_dlg<P: AsRef<Path>>(
    path: P,
    title: &str,
    frames: &[LambdaDerivativeFrame],
) -> Result<(), IoError> {
    let mut writer = DlgWriter::new(path, title)?;

    for frame in frames {
        writer.write_frame(frame)?;
    }

    writer.close()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lambda_derivative_frame() {
        let mut frame = LambdaDerivativeFrame::new(0.0, 0.5);
        frame.bond = 10.0;
        frame.angle = 5.0;
        frame.lj = -15.0;
        frame.coulomb = 20.0;
        frame.update_total();

        assert!((frame.total - 20.0).abs() < 1e-6);
        assert!((frame.lambda - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_dlg_writer() {
        let mut frame1 = LambdaDerivativeFrame::new(0.0, 0.0);
        frame1.bond = 1.0;
        frame1.angle = 2.0;
        frame1.lj = 3.0;
        frame1.update_total();

        let mut frame2 = LambdaDerivativeFrame::new(0.002, 0.1);
        frame2.bond = 1.5;
        frame2.angle = 2.5;
        frame2.lj = 3.5;
        frame2.update_total();

        let frames = vec![frame1, frame2];

        let result = write_dlg("/tmp/test.dlg", "Test FEP simulation", &frames);
        assert!(result.is_ok());

        // Read back and verify
        use std::fs;
        let content = fs::read_to_string("/tmp/test.dlg").unwrap();
        assert!(content.contains("TITLE"));
        assert!(content.contains("FREEENERGYDERIVATIVES"));
        assert!(content.contains("# Block 0: Time"));
        assert!(content.contains("# Block 2: dH/dλ total"));
        assert!(content.contains("END"));

        fs::remove_file("/tmp/test.dlg").ok();
    }

    #[test]
    fn test_set_bonded() {
        let mut frame = LambdaDerivativeFrame::new(1.0, 0.5);
        frame.set_bonded(10.0, 20.0, 30.0, 40.0);
        frame.update_total();

        assert!((frame.bond - 10.0).abs() < 1e-6);
        assert!((frame.angle - 20.0).abs() < 1e-6);
        assert!((frame.improper - 30.0).abs() < 1e-6);
        assert!((frame.dihedral - 40.0).abs() < 1e-6);
        assert!((frame.total - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_set_nonbonded() {
        let mut frame = LambdaDerivativeFrame::new(1.0, 0.5);
        frame.set_nonbonded(-50.0, 75.0);
        frame.update_total();

        assert!((frame.lj + 50.0).abs() < 1e-6);
        assert!((frame.coulomb - 75.0).abs() < 1e-6);
        assert!((frame.total - 25.0).abs() < 1e-6);
    }
}
