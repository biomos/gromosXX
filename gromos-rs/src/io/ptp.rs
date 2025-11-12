//! Perturbation Topology (.ptp) file writer
//!
//! The .ptp format specifies perturbed atoms and interactions for FEP calculations.
//! It defines dual-topology parameters (state A and state B) for atoms that change
//! during the simulation as lambda goes from 0 to 1.
//!
//! # Format
//! ```text
//! TITLE
//!   Perturbation topology for FEP calculation
//! END
//! PERTURBEDATOM
//! #  NR  IACNA  IACNB  MASNA  MASNB  CHARGA  CHARGB  ALPHLJ  ALPHCRF
//!     1      1      2   12.01   14.01   -0.100    0.200   1.510     0.500
//! END
//! PERTURBATIONPARAMETERS
//! # ALPHLJ  ALPHCRF
//!    1.510    0.500
//! END
//! ```

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use crate::topology::Topology;

/// Writer for GROMOS perturbation topology (.ptp) files
pub struct PtpWriter {
    file: File,
}

impl PtpWriter {
    /// Create a new .ptp file writer
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self { file })
    }

    /// Write perturbation topology comparing two states
    ///
    /// Compares topology_a (state A, λ=0) with topology_b (state B, λ=1)
    /// and writes the differences as a perturbation topology.
    pub fn write(
        &mut self,
        topology_a: &Topology,
        topology_b: &Topology,
        alpha_lj: f64,
        alpha_crf: f64,
        title: &str,
    ) -> io::Result<()> {
        // Check that topologies have same number of atoms
        if topology_a.num_atoms() != topology_b.num_atoms() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Topologies must have same number of atoms: {} vs {}",
                    topology_a.num_atoms(),
                    topology_b.num_atoms()
                ),
            ));
        }

        // Write title
        writeln!(self.file, "TITLE")?;
        writeln!(self.file, "{}", title)?;
        writeln!(self.file, "END")?;

        // Write perturbed atoms
        self.write_perturbed_atoms(topology_a, topology_b, alpha_lj, alpha_crf)?;

        // Write perturbation parameters
        writeln!(self.file, "PERTURBATIONPARAMETERS")?;
        writeln!(self.file, "# ALPHLJ  ALPHCRF")?;
        writeln!(self.file, "  {:8.3}  {:8.3}", alpha_lj, alpha_crf)?;
        writeln!(self.file, "END")?;

        Ok(())
    }

    /// Write perturbed atoms block
    fn write_perturbed_atoms(
        &mut self,
        topology_a: &Topology,
        topology_b: &Topology,
        alpha_lj: f64,
        alpha_crf: f64,
    ) -> io::Result<()> {
        writeln!(self.file, "PERTURBEDATOM")?;
        writeln!(
            self.file,
            "#  NR  IACNA  IACNB  MASNA  MASNB  CHARGA  CHARGB  ALPHLJ  ALPHCRF"
        )?;

        let n_atoms = topology_a.num_atoms();

        for i in 0..n_atoms {
            let iac_a = topology_a.iac[i];
            let iac_b = topology_b.iac[i];
            let mass_a = topology_a.mass[i];
            let mass_b = topology_b.mass[i];
            let charge_a = topology_a.charge[i];
            let charge_b = topology_b.charge[i];

            // Only write atoms that are actually perturbed
            let is_perturbed = iac_a != iac_b
                || (mass_a - mass_b).abs() > 1e-6
                || (charge_a - charge_b).abs() > 1e-6;

            if is_perturbed {
                writeln!(
                    self.file,
                    "{:5} {:6} {:6} {:7.2} {:7.2} {:8.3} {:8.3} {:7.3} {:8.3}",
                    i + 1,        // 1-based atom number
                    iac_a,        // IAC state A
                    iac_b,        // IAC state B
                    mass_a,       // Mass state A
                    mass_b,       // Mass state B
                    charge_a,     // Charge state A
                    charge_b,     // Charge state B
                    alpha_lj,     // LJ softness
                    alpha_crf     // CRF softness
                )?;
            }
        }

        writeln!(self.file, "END")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::Atom;
    use std::fs;

    #[test]
    fn test_ptp_writer() {
        let mut topo_a = Topology::new();
        let mut topo_b = Topology::new();

        // Add 2 atoms to both topologies
        for i in 0..2 {
            topo_a.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "MOL".to_string(),
                iac: 1,
                mass: 12.01,
                charge: -0.1,
                is_perturbed: true,
                is_polarisable: false,
                is_coarse_grained: false,
            });

            topo_b.solute.atoms.push(Atom {
                name: format!("N{}", i),
                residue_nr: 1,
                residue_name: "MOL".to_string(),
                iac: 2,
                mass: 14.01,
                charge: 0.2,
                is_perturbed: true,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }

        topo_a.iac = vec![1, 1];
        topo_a.mass = vec![12.01, 12.01];
        topo_a.charge = vec![-0.1, -0.1];

        topo_b.iac = vec![2, 2];
        topo_b.mass = vec![14.01, 14.01];
        topo_b.charge = vec![0.2, 0.2];

        let path = "/tmp/test.ptp";
        let mut writer = PtpWriter::new(path).unwrap();
        writer
            .write(&topo_a, &topo_b, 1.51, 0.5, "Test perturbation topology")
            .unwrap();

        // Read file and verify
        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("TITLE"));
        assert!(content.contains("PERTURBEDATOM"));
        assert!(content.contains("PERTURBATIONPARAMETERS"));
        assert!(content.contains("1.510"));
        assert!(content.contains("0.500"));

        fs::remove_file(path).ok();
    }
}
