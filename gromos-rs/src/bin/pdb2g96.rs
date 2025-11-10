//! pdb2g96 - Convert PDB files to GROMOS96 format
//!
//! Usage: pdb2g96 <input.pdb> [output.g96]
//!
//! Converts Protein Data Bank (PDB) format coordinate files to GROMOS96 format.
//! - Automatically converts Angstrom to nanometers
//! - Preserves atom and residue information
//! - Outputs to stdout if no output file specified

use gromos_rs::io::pdb::PDBStructure;
use gromos_rs::io::g96::G96Writer;
use std::env;
use std::process;

fn print_usage() {
    eprintln!("Usage: pdb2g96 <input.pdb> [output.g96]");
    eprintln!();
    eprintln!("Convert PDB format to GROMOS96 format");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  input.pdb    Input PDB file");
    eprintln!("  output.g96   Output GROMOS96 file (optional, default: stdout)");
    eprintln!();
    eprintln!("Features:");
    eprintln!("  - Automatic Angstrom to nanometer conversion");
    eprintln!("  - Preserves residue and atom information");
    eprintln!("  - Supports ATOM and HETATM records");
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let input_pdb = &args[1];
    let output_g96 = if args.len() > 2 {
        Some(&args[2])
    } else {
        None
    };

    // Read PDB file
    println!("Reading PDB file: {}", input_pdb);
    let pdb = match PDBStructure::read_pdb(input_pdb) {
        Ok(pdb) => pdb,
        Err(e) => {
            eprintln!("Error reading PDB file: {}", e);
            process::exit(1);
        }
    };

    let num_atoms = pdb.num_atoms();
    let num_residues = pdb.residues.len();
    println!("  Found {} atoms in {} residues", num_atoms, num_residues);

    // Collect all coordinates
    let mut positions = Vec::with_capacity(num_atoms);
    for residue in &pdb.residues {
        for atom in &residue.atoms {
            positions.push(atom.coord);
        }
    }

    // Create title
    let title = if pdb.title.is_empty() {
        format!("pdb2g96: Converted from {}", input_pdb)
    } else {
        format!("pdb2g96: {}", pdb.title.join(" "))
    };

    // Write GROMOS96 file
    let default_output = "output.g96".to_string();
    let output_path = output_g96.unwrap_or(&default_output);
    println!("Writing GROMOS96 file: {}", output_path);

    match write_g96_from_pdb(output_path, &title, &pdb) {
        Ok(_) => {
            println!("Successfully converted {} to {}", input_pdb, output_path);
            println!("  {} atoms written", positions.len());
        }
        Err(e) => {
            eprintln!("Error writing GROMOS96 file: {}", e);
            process::exit(1);
        }
    }
}

/// Write a GROMOS96 file from PDB structure
fn write_g96_from_pdb(
    path: &str,
    title: &str,
    pdb: &PDBStructure,
) -> Result<(), String> {
    let mut writer = G96Writer::new(path)?;

    writer.write_title(title)?;

    // Write POSITION block
    use std::io::Write;
    writeln!(writer.writer.get_mut(), "POSITION")
        .map_err(|e| format!("Write error: {}", e))?;

    let mut atom_counter = 1;
    for residue in &pdb.residues {
        for atom in &residue.atoms {
            // GROMOS96 format: resnum resname atomname atomnum x y z
            writeln!(
                writer.writer.get_mut(),
                "{:5} {:5} {:5}{:7}{:15.9}{:15.9}{:15.9}",
                residue.number,
                residue.name,
                atom.name,
                atom_counter,
                atom.coord.x as f64,
                atom.coord.y as f64,
                atom.coord.z as f64
            )
            .map_err(|e| format!("Write error: {}", e))?;
            atom_counter += 1;
        }
    }

    writeln!(writer.writer.get_mut(), "END")
        .map_err(|e| format!("Write error: {}", e))?;

    writer.close()
}
