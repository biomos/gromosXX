//! Validation utilities for GROMOS simulations
//!
//! Provides validation functions for topology, coordinates, energies,
//! and other simulation data to catch errors early.

use crate::math::Vec3;
use crate::topology::Topology;
use crate::configuration::Configuration;

#[derive(Debug)]
pub struct ValidationError {
    pub level: ValidationLevel,
    pub message: String,
    pub location: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationLevel {
    Warning,
    Error,
    Fatal,
}

impl ValidationError {
    pub fn warning(location: &str, message: String) -> Self {
        Self {
            level: ValidationLevel::Warning,
            message,
            location: location.to_string(),
        }
    }

    pub fn error(location: &str, message: String) -> Self {
        Self {
            level: ValidationLevel::Error,
            message,
            location: location.to_string(),
        }
    }

    pub fn fatal(location: &str, message: String) -> Self {
        Self {
            level: ValidationLevel::Fatal,
            message,
            location: location.to_string(),
        }
    }

    pub fn print(&self) {
        let level_str = match self.level {
            ValidationLevel::Warning => "WARNING",
            ValidationLevel::Error => "ERROR",
            ValidationLevel::Fatal => "FATAL",
        };
        eprintln!("[{}] {}: {}", level_str, self.location, self.message);
    }
}

pub struct ValidationReport {
    pub errors: Vec<ValidationError>,
    pub warnings: Vec<ValidationError>,
}

impl ValidationReport {
    pub fn new() -> Self {
        Self {
            errors: Vec::new(),
            warnings: Vec::new(),
        }
    }

    pub fn add(&mut self, error: ValidationError) {
        match error.level {
            ValidationLevel::Warning => self.warnings.push(error),
            ValidationLevel::Error | ValidationLevel::Fatal => self.errors.push(error),
        }
    }

    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }

    pub fn has_fatal(&self) -> bool {
        self.errors.iter().any(|e| e.level == ValidationLevel::Fatal)
    }

    pub fn print(&self) {
        for warning in &self.warnings {
            warning.print();
        }
        for error in &self.errors {
            error.print();
        }
    }

    pub fn print_summary(&self) {
        if !self.warnings.is_empty() {
            eprintln!("\n{} warning(s) found", self.warnings.len());
        }
        if !self.errors.is_empty() {
            eprintln!("{} error(s) found", self.errors.len());
        }
    }
}

impl Default for ValidationReport {
    fn default() -> Self {
        Self::new()
    }
}

/// Validate topology consistency
pub fn validate_topology(topo: &Topology) -> ValidationReport {
    let mut report = ValidationReport::new();

    // Check atom count
    let n_atoms = topo.num_atoms();
    if n_atoms == 0 {
        report.add(ValidationError::fatal("topology", "No atoms in topology".to_string()));
        return report;
    }

    // Check masses
    for (i, &mass) in topo.mass.iter().enumerate() {
        if mass <= 0.0 {
            report.add(ValidationError::error(
                "topology",
                format!("Atom {} has non-positive mass: {}", i + 1, mass),
            ));
        }
        if mass > 1000.0 {
            report.add(ValidationError::warning(
                "topology",
                format!("Atom {} has unusually large mass: {}", i + 1, mass),
            ));
        }
    }

    // Check bonds
    for (i, bond) in topo.solute.bonds.iter().enumerate() {
        if bond.i >= n_atoms || bond.j >= n_atoms {
            report.add(ValidationError::error(
                "topology",
                format!("Bond {} references out-of-range atoms: {}-{}", i, bond.i + 1, bond.j + 1),
            ));
        }
        if bond.i == bond.j {
            report.add(ValidationError::error(
                "topology",
                format!("Bond {} has same atom twice: {}", i, bond.i + 1),
            ));
        }
    }

    // Check angles
    for (i, angle) in topo.solute.angles.iter().enumerate() {
        if angle.i >= n_atoms || angle.j >= n_atoms || angle.k >= n_atoms {
            report.add(ValidationError::error(
                "topology",
                format!("Angle {} references out-of-range atoms: {}-{}-{}",
                    i, angle.i + 1, angle.j + 1, angle.k + 1),
            ));
        }
        if angle.i == angle.j || angle.j == angle.k || angle.i == angle.k {
            report.add(ValidationError::error(
                "topology",
                format!("Angle {} has duplicate atoms: {}-{}-{}",
                    i, angle.i + 1, angle.j + 1, angle.k + 1),
            ));
        }
    }

    // Check proper dihedrals
    for (i, dih) in topo.solute.proper_dihedrals.iter().enumerate() {
        if dih.i >= n_atoms || dih.j >= n_atoms || dih.k >= n_atoms || dih.l >= n_atoms {
            report.add(ValidationError::error(
                "topology",
                format!("Dihedral {} references out-of-range atoms: {}-{}-{}-{}",
                    i, dih.i + 1, dih.j + 1, dih.k + 1, dih.l + 1),
            ));
        }
    }

    report
}

/// Validate coordinates
pub fn validate_coordinates(positions: &[Vec3], box_dims: Option<Vec3>) -> ValidationReport {
    let mut report = ValidationReport::new();

    if positions.is_empty() {
        report.add(ValidationError::fatal("coordinates", "No positions found".to_string()));
        return report;
    }

    // Check for NaN or infinite values
    for (i, pos) in positions.iter().enumerate() {
        if !pos.x.is_finite() || !pos.y.is_finite() || !pos.z.is_finite() {
            report.add(ValidationError::error(
                "coordinates",
                format!("Atom {} has non-finite coordinates: ({}, {}, {})",
                    i + 1, pos.x, pos.y, pos.z),
            ));
        }
    }

    // Check distances between consecutive atoms (catch overlaps)
    let mut min_dist = f32::MAX;
    let mut min_pair = (0, 0);
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len().min(i + 100) {
            let dist = (positions[j] - positions[i]).length();
            if dist < min_dist {
                min_dist = dist;
                min_pair = (i, j);
            }
            if dist < 0.01 {  // Less than 0.01 nm = 0.1 Å
                report.add(ValidationError::warning(
                    "coordinates",
                    format!("Atoms {} and {} are very close: {:.4} nm",
                        i + 1, j + 1, dist),
                ));
            }
        }
    }

    // Check if atoms are far outside the box
    if let Some(box_dims) = box_dims {
        for (i, pos) in positions.iter().enumerate() {
            if pos.x.abs() > box_dims.x * 2.0
                || pos.y.abs() > box_dims.y * 2.0
                || pos.z.abs() > box_dims.z * 2.0 {
                report.add(ValidationError::warning(
                    "coordinates",
                    format!("Atom {} is far outside box: ({:.3}, {:.3}, {:.3})",
                        i + 1, pos.x, pos.y, pos.z),
                ));
            }
        }
    }

    report
}

/// Validate configuration (topology + coordinates)
pub fn validate_configuration(topo: &Topology, conf: &Configuration) -> ValidationReport {
    let mut report = ValidationReport::new();

    let n_topo = topo.num_atoms();
    let n_conf = conf.current().pos.len();

    if n_topo != n_conf {
        report.add(ValidationError::fatal(
            "configuration",
            format!("Topology has {} atoms but coordinates have {}", n_topo, n_conf),
        ));
        return report;
    }

    // Check for reasonable bond lengths
    for bond in &topo.solute.bonds {
        let pos_i = conf.current().pos[bond.i];
        let pos_j = conf.current().pos[bond.j];
        let r = (pos_j - pos_i).length();

        // Typical bond lengths: 0.1-0.2 nm
        if r < 0.05 || r > 0.5 {
            report.add(ValidationError::warning(
                "configuration",
                format!("Bond {}-{} has unusual length: {:.4} nm",
                    bond.i + 1, bond.j + 1, r),
            ));
        }
    }

    // Check for reasonable angles
    for angle in &topo.solute.angles {
        let pos_i = conf.current().pos[angle.i];
        let pos_j = conf.current().pos[angle.j];
        let pos_k = conf.current().pos[angle.k];

        let v1 = (pos_i - pos_j).normalize();
        let v2 = (pos_k - pos_j).normalize();
        let cos_theta = v1.dot(v2).clamp(-1.0, 1.0);
        let theta = cos_theta.acos().to_degrees();

        // Angles should be between 60° and 180°
        if theta < 60.0 || theta > 180.0 {
            report.add(ValidationError::warning(
                "configuration",
                format!("Angle {}-{}-{} has unusual value: {:.1}°",
                    angle.i + 1, angle.j + 1, angle.k + 1, theta),
            ));
        }
    }

    report
}

/// Validate energy values
pub fn validate_energy(
    kinetic: f64,
    potential: f64,
    total: f64,
    temperature: f64,
) -> ValidationReport {
    let mut report = ValidationReport::new();

    // Check for NaN or infinite
    if !kinetic.is_finite() {
        report.add(ValidationError::error("energy", "Kinetic energy is not finite".to_string()));
    }
    if !potential.is_finite() {
        report.add(ValidationError::error("energy", "Potential energy is not finite".to_string()));
    }
    if !total.is_finite() {
        report.add(ValidationError::error("energy", "Total energy is not finite".to_string()));
    }
    if !temperature.is_finite() {
        report.add(ValidationError::error("energy", "Temperature is not finite".to_string()));
    }

    // Check for unreasonable values
    if temperature < 0.0 {
        report.add(ValidationError::error(
            "energy",
            format!("Negative temperature: {:.2} K", temperature),
        ));
    }
    if temperature > 10000.0 {
        report.add(ValidationError::warning(
            "energy",
            format!("Very high temperature: {:.2} K", temperature),
        ));
    }

    // Check energy conservation (rough)
    let computed_total = kinetic + potential;
    let energy_diff = (computed_total - total).abs();
    if energy_diff > 0.01 * total.abs() {
        report.add(ValidationError::warning(
            "energy",
            format!("Energy sum mismatch: K+V = {:.4}, Total = {:.4}", computed_total, total),
        ));
    }

    report
}

/// Check for energy drift over trajectory
pub fn check_energy_drift(energies: &[(f64, f64, f64)]) -> ValidationReport {
    let mut report = ValidationReport::new();

    if energies.len() < 10 {
        return report;
    }

    let initial_total = energies[0].2;
    let final_total = energies[energies.len() - 1].2;
    let drift = (final_total - initial_total).abs();
    let relative_drift = drift / initial_total.abs();

    if relative_drift > 0.01 {  // 1% drift
        report.add(ValidationError::warning(
            "energy_drift",
            format!("Significant energy drift: {:.2}% ({:.4} kJ/mol over {} frames)",
                relative_drift * 100.0, drift, energies.len()),
        ));
    }

    // Check for sudden jumps
    for i in 1..energies.len() {
        let prev_total = energies[i - 1].2;
        let curr_total = energies[i].2;
        let jump = (curr_total - prev_total).abs();
        let relative_jump = jump / prev_total.abs();

        if relative_jump > 0.1 {  // 10% jump between frames
            report.add(ValidationError::error(
                "energy_jump",
                format!("Large energy jump at frame {}: {:.2}% ({:.4} kJ/mol)",
                    i, relative_jump * 100.0, jump),
            ));
        }
    }

    report
}
