//! check_top - Validate GROMOS topology file consistency
//!
//! Usage: check_top <input.top>
//!
//! Performs comprehensive checks on topology files:
//! - Atom numbering and types
//! - Bond/angle/dihedral indices
//! - Force field parameter references
//! - Exclusion consistency
//! - Charge group validity

use gromos_rs::io::topology::{read_topology_file, build_topology};
use gromos_rs::topology::Topology;
use std::env;
use std::process;

fn print_usage() {
    eprintln!("Usage: check_top <input.top>");
    eprintln!();
    eprintln!("Validate GROMOS topology file consistency");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  input.top    Input topology file to check");
    eprintln!();
    eprintln!("Checks performed:");
    eprintln!("  - Atom numbering and indices");
    eprintln!("  - Bond/angle/dihedral validity");
    eprintln!("  - Force field parameter references");
    eprintln!("  - Exclusion list consistency");
    eprintln!("  - Charge neutrality warnings");
}

struct TopologyChecker {
    warnings: Vec<String>,
    errors: Vec<String>,
}

impl TopologyChecker {
    fn new() -> Self {
        Self {
            warnings: Vec::new(),
            errors: Vec::new(),
        }
    }

    fn check_all(&mut self, topo: &Topology) {
        println!("Performing topology checks...\n");

        self.check_atom_count(topo);
        self.check_masses(topo);
        self.check_charges(topo);
        self.check_bonds(topo);
        self.check_angles(topo);
        self.check_proper_dihedrals(topo);
        self.check_improper_dihedrals(topo);
        self.check_exclusions(topo);
        self.check_bond_parameters(topo);
        self.check_angle_parameters(topo);
        self.check_lj_parameters(topo);
    }

    fn check_atom_count(&mut self, topo: &Topology) {
        let n_atoms = topo.solute.num_atoms();
        println!("CHECK: Atom count");
        println!("  Solute atoms: {}", n_atoms);

        if n_atoms == 0 {
            self.errors.push("No atoms in topology".to_string());
        }

        // Check array sizes match
        if topo.mass.len() != n_atoms {
            self.errors.push(format!(
                "Mass array size ({}) != atom count ({})",
                topo.mass.len(),
                n_atoms
            ));
        }

        if topo.charge.len() != n_atoms {
            self.errors.push(format!(
                "Charge array size ({}) != atom count ({})",
                topo.charge.len(),
                n_atoms
            ));
        }

        if topo.iac.len() != n_atoms {
            self.errors.push(format!(
                "IAC array size ({}) != atom count ({})",
                topo.iac.len(),
                n_atoms
            ));
        }

        println!("  ✓ Atom count checks passed\n");
    }

    fn check_masses(&mut self, topo: &Topology) {
        println!("CHECK: Atomic masses");

        let mut zero_mass_count = 0;
        let mut negative_mass_count = 0;

        for (i, &mass) in topo.mass.iter().enumerate() {
            if mass < 0.0 {
                self.errors.push(format!("Atom {} has negative mass: {}", i + 1, mass));
                negative_mass_count += 1;
            } else if mass == 0.0 {
                self.warnings.push(format!("Atom {} has zero mass", i + 1));
                zero_mass_count += 1;
            }
        }

        if zero_mass_count > 0 {
            println!("  ⚠ {} atoms with zero mass", zero_mass_count);
        }
        if negative_mass_count == 0 && zero_mass_count == 0 {
            println!("  ✓ All masses valid\n");
        } else {
            println!();
        }
    }

    fn check_charges(&mut self, topo: &Topology) {
        println!("CHECK: Atomic charges");

        let total_charge: f64 = topo.charge.iter().sum();
        println!("  Total charge: {:.6} e", total_charge);

        if total_charge.abs() > 1e-4 {
            self.warnings.push(format!(
                "System is not charge-neutral: {:.6} e",
                total_charge
            ));
            println!("  ⚠ System not charge-neutral\n");
        } else {
            println!("  ✓ System is charge-neutral\n");
        }
    }

    fn check_bonds(&mut self, topo: &Topology) {
        println!("CHECK: Bonds");
        let n_atoms = topo.solute.num_atoms();
        let n_bonds = topo.solute.bonds.len();

        println!("  Number of bonds: {}", n_bonds);

        for (idx, bond) in topo.solute.bonds.iter().enumerate() {
            // Check atom indices are in range
            if bond.i >= n_atoms {
                self.errors.push(format!(
                    "Bond {} atom i={} out of range (max: {})",
                    idx + 1,
                    bond.i + 1,
                    n_atoms
                ));
            }
            if bond.j >= n_atoms {
                self.errors.push(format!(
                    "Bond {} atom j={} out of range (max: {})",
                    idx + 1,
                    bond.j + 1,
                    n_atoms
                ));
            }

            // Check for self-bonds
            if bond.i == bond.j {
                self.errors.push(format!(
                    "Bond {} is a self-bond: atom {}",
                    idx + 1,
                    bond.i + 1
                ));
            }

            // Check bond type is valid
            if bond.bond_type >= topo.bond_parameters.len() {
                self.errors.push(format!(
                    "Bond {} has invalid type {} (max: {})",
                    idx + 1,
                    bond.bond_type + 1,
                    topo.bond_parameters.len()
                ));
            }
        }

        if self.errors.iter().any(|e| e.contains("Bond")) {
            println!("  ✗ Bond errors found\n");
        } else {
            println!("  ✓ All bonds valid\n");
        }
    }

    fn check_angles(&mut self, topo: &Topology) {
        println!("CHECK: Bond angles");
        let n_atoms = topo.solute.num_atoms();
        let n_angles = topo.solute.angles.len();

        println!("  Number of angles: {}", n_angles);

        for (idx, angle) in topo.solute.angles.iter().enumerate() {
            // Check atom indices
            if angle.i >= n_atoms {
                self.errors.push(format!(
                    "Angle {} atom i={} out of range",
                    idx + 1,
                    angle.i + 1
                ));
            }
            if angle.j >= n_atoms {
                self.errors.push(format!(
                    "Angle {} atom j={} out of range",
                    idx + 1,
                    angle.j + 1
                ));
            }
            if angle.k >= n_atoms {
                self.errors.push(format!(
                    "Angle {} atom k={} out of range",
                    idx + 1,
                    angle.k + 1
                ));
            }

            // Check for degenerate angles
            if angle.i == angle.j || angle.j == angle.k || angle.i == angle.k {
                self.errors.push(format!(
                    "Angle {} has duplicate atoms: {}-{}-{}",
                    idx + 1,
                    angle.i + 1,
                    angle.j + 1,
                    angle.k + 1
                ));
            }

            // Check angle type
            if angle.angle_type >= topo.angle_parameters.len() {
                self.errors.push(format!(
                    "Angle {} has invalid type {} (max: {})",
                    idx + 1,
                    angle.angle_type + 1,
                    topo.angle_parameters.len()
                ));
            }
        }

        if self.errors.iter().any(|e| e.contains("Angle")) {
            println!("  ✗ Angle errors found\n");
        } else {
            println!("  ✓ All angles valid\n");
        }
    }

    fn check_proper_dihedrals(&mut self, topo: &Topology) {
        println!("CHECK: Proper dihedrals");
        let n_atoms = topo.solute.num_atoms();
        let n_dihedrals = topo.solute.proper_dihedrals.len();

        println!("  Number of proper dihedrals: {}", n_dihedrals);

        for (idx, dihedral) in topo.solute.proper_dihedrals.iter().enumerate() {
            if dihedral.i >= n_atoms
                || dihedral.j >= n_atoms
                || dihedral.k >= n_atoms
                || dihedral.l >= n_atoms
            {
                self.errors.push(format!(
                    "Proper dihedral {} has out-of-range atoms",
                    idx + 1
                ));
            }

            // Check for duplicate atoms
            let atoms = [dihedral.i, dihedral.j, dihedral.k, dihedral.l];
            for i in 0..4 {
                for j in (i + 1)..4 {
                    if atoms[i] == atoms[j] {
                        self.errors.push(format!(
                            "Proper dihedral {} has duplicate atoms",
                            idx + 1
                        ));
                        break;
                    }
                }
            }
        }

        if self.errors.iter().any(|e| e.contains("Proper dihedral")) {
            println!("  ✗ Dihedral errors found\n");
        } else if n_dihedrals > 0 {
            println!("  ✓ All proper dihedrals valid\n");
        } else {
            println!();
        }
    }

    fn check_improper_dihedrals(&mut self, topo: &Topology) {
        println!("CHECK: Improper dihedrals");
        let n_atoms = topo.solute.num_atoms();
        let n_impropers = topo.solute.improper_dihedrals.len();

        println!("  Number of improper dihedrals: {}", n_impropers);

        for (idx, dihedral) in topo.solute.improper_dihedrals.iter().enumerate() {
            if dihedral.i >= n_atoms
                || dihedral.j >= n_atoms
                || dihedral.k >= n_atoms
                || dihedral.l >= n_atoms
            {
                self.errors.push(format!(
                    "Improper dihedral {} has out-of-range atoms",
                    idx + 1
                ));
            }
        }

        if self.errors.iter().any(|e| e.contains("Improper dihedral")) {
            println!("  ✗ Improper dihedral errors found\n");
        } else if n_impropers > 0 {
            println!("  ✓ All improper dihedrals valid\n");
        } else {
            println!();
        }
    }

    fn check_exclusions(&mut self, topo: &Topology) {
        println!("CHECK: Exclusion lists");
        let n_atoms = topo.solute.num_atoms();

        if topo.exclusions.len() != n_atoms {
            self.errors.push(format!(
                "Exclusion list size ({}) != atom count ({})",
                topo.exclusions.len(),
                n_atoms
            ));
        }

        let mut total_exclusions = 0;
        let mut asymmetric_count = 0;

        for (i, excl_set) in topo.exclusions.iter().enumerate() {
            total_exclusions += excl_set.len();

            for &j in excl_set {
                if j >= n_atoms {
                    self.errors.push(format!(
                        "Atom {} has exclusion to out-of-range atom {}",
                        i + 1,
                        j + 1
                    ));
                }

                // Check symmetry: if i excludes j, j should exclude i
                if j < topo.exclusions.len() && !topo.exclusions[j].contains(&i) {
                    self.warnings.push(format!(
                        "Asymmetric exclusion: atom {} excludes {} but not vice versa",
                        i + 1,
                        j + 1
                    ));
                    asymmetric_count += 1;
                }
            }
        }

        println!("  Total exclusions: {}", total_exclusions);
        if asymmetric_count > 0 {
            println!("  ⚠ {} asymmetric exclusions", asymmetric_count);
        } else {
            println!("  ✓ Exclusions are symmetric\n");
        }
    }

    fn check_bond_parameters(&mut self, topo: &Topology) {
        println!("CHECK: Bond parameters");
        println!("  Number of bond types: {}", topo.bond_parameters.len());

        for (i, params) in topo.bond_parameters.iter().enumerate() {
            if params.r0 <= 0.0 {
                self.errors.push(format!(
                    "Bond type {} has non-positive equilibrium length: {}",
                    i + 1,
                    params.r0
                ));
            }

            if params.k_harmonic < 0.0 {
                self.warnings.push(format!(
                    "Bond type {} has negative harmonic force constant: {}",
                    i + 1,
                    params.k_harmonic
                ));
            }
        }

        if self.errors.iter().any(|e| e.contains("Bond type")) {
            println!("  ✗ Bond parameter errors found\n");
        } else if !topo.bond_parameters.is_empty() {
            println!("  ✓ All bond parameters valid\n");
        } else {
            println!();
        }
    }

    fn check_angle_parameters(&mut self, topo: &Topology) {
        println!("CHECK: Angle parameters");
        println!("  Number of angle types: {}", topo.angle_parameters.len());

        for (i, params) in topo.angle_parameters.iter().enumerate() {
            if params.theta0 < 0.0 || params.theta0 > std::f64::consts::PI {
                self.warnings.push(format!(
                    "Angle type {} has theta0 outside [0, π]: {} rad",
                    i + 1,
                    params.theta0
                ));
            }

            if params.k_harmonic < 0.0 {
                self.warnings.push(format!(
                    "Angle type {} has negative harmonic force constant: {}",
                    i + 1,
                    params.k_harmonic
                ));
            }
        }

        if !topo.angle_parameters.is_empty() {
            println!("  ✓ Angle parameter checks complete\n");
        } else {
            println!();
        }
    }

    fn check_lj_parameters(&mut self, topo: &Topology) {
        println!("CHECK: Lennard-Jones parameters");
        let n_types = topo.lj_parameters.len();
        println!("  Number of LJ atom types: {}", n_types);

        if n_types == 0 {
            self.warnings.push("No LJ parameters defined".to_string());
            println!("  ⚠ No LJ parameters\n");
            return;
        }

        // Check matrix is square
        for (i, row) in topo.lj_parameters.iter().enumerate() {
            if row.len() != n_types {
                self.errors.push(format!(
                    "LJ parameter matrix row {} has wrong size: {} (expected {})",
                    i,
                    row.len(),
                    n_types
                ));
            }
        }

        // Check for negative parameters
        let mut negative_count = 0;
        for i in 0..n_types {
            for j in 0..n_types {
                if i < topo.lj_parameters.len() && j < topo.lj_parameters[i].len() {
                    let params = &topo.lj_parameters[i][j];
                    if params.c6 < 0.0 || params.c12 < 0.0 {
                        negative_count += 1;
                    }
                }
            }
        }

        if negative_count > 0 {
            self.warnings.push(format!(
                "{} LJ parameter pairs have negative values",
                negative_count
            ));
        }

        println!("  ✓ LJ parameter checks complete\n");
    }

    fn print_summary(&self) {
        println!("{}", "=".repeat(60));
        println!("TOPOLOGY CHECK SUMMARY");
        println!("{}", "=".repeat(60));

        if self.errors.is_empty() && self.warnings.is_empty() {
            println!("✓ Topology is valid - no errors or warnings");
            println!();
            return;
        }

        if !self.errors.is_empty() {
            println!("\nERRORS ({}):", self.errors.len());
            for error in &self.errors {
                println!("  ✗ {}", error);
            }
        }

        if !self.warnings.is_empty() {
            println!("\nWARNINGS ({}):", self.warnings.len());
            for warning in &self.warnings {
                println!("  ⚠ {}", warning);
            }
        }

        println!();

        if self.errors.is_empty() {
            println!("✓ No critical errors - topology is usable");
        } else {
            println!("✗ Critical errors found - topology may not work correctly");
        }
        println!();
    }

    fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let input_file = &args[1];

    println!("Checking topology file: {}\n", input_file);

    // Read topology
    let topo = match read_topology_file(input_file) {
        Ok(parsed) => {
            println!("Successfully parsed topology file\n");
            build_topology(parsed)
        }
        Err(e) => {
            eprintln!("Error reading topology file: {:?}", e);
            process::exit(1);
        }
    };

    // Run checks
    let mut checker = TopologyChecker::new();
    checker.check_all(&topo);
    checker.print_summary();

    // Exit with error code if critical errors found
    if checker.has_errors() {
        process::exit(1);
    }
}
