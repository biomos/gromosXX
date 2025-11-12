//! sim_box - Solvate a solute in a box of pre-equilibrated solvent
//!
//! Usage: sim_box @topo <topology> @pbc <r|t> @pos <solute.cnf> @solvent <solvent.g96>
//!               [@minwall <distance>] [@thresh <distance>] [@boxsize]
//!
//! Solvates a molecular system by:
//! 1. Reading solute coordinates
//! 2. Replicating a pre-equilibrated solvent box
//! 3. Removing solvent molecules that clash with the solute
//! 4. Writing the solvated system to stdout

use gromos_rs::io::g96::write_g96;
use gromos_rs::math::Vec3;
use std::env;
use std::process;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn print_usage() {
    eprintln!("sim_box - Solvate a solute in a box of pre-equilibrated solvent");
    eprintln!();
    eprintln!("Usage: sim_box @topo <topology> @pbc <r|t> @pos <solute.cnf> @solvent <solvent.g96> \\");
    eprintln!("               [@minwall <distance>] [@thresh <distance>] [@boxsize]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo      Molecular topology file");
    eprintln!("  @pbc       Periodic boundary conditions:");
    eprintln!("             r = rectangular box");
    eprintln!("             t = truncated octahedron (not yet implemented)");
    eprintln!("  @pos       Input coordinate file for the solute");
    eprintln!("  @solvent   Input coordinate file for pre-equilibrated solvent");
    eprintln!("  @minwall   Minimum solute-to-wall distance (nm)");
    eprintln!("             - 1 value: cubic box");
    eprintln!("             - 3 values: rectangular box (x y z)");
    eprintln!("  @thresh    Minimum solvent-to-solute distance (default: 0.23 nm)");
    eprintln!("  @boxsize   Use box dimensions from solute coordinate file");
    eprintln!();
    eprintln!("Description:");
    eprintln!("  Solvates a solute molecule in a pre-equilibrated solvent box.");
    eprintln!("  The solvent configuration should contain a BOX block.");
    eprintln!();
    eprintln!("  Box sizing options:");
    eprintln!("  1. @boxsize: Use box from solute file");
    eprintln!("  2. @minwall: Calculate box from solute extent + minwall distance");
    eprintln!();
    eprintln!("  Solvent removal criteria:");
    eprintln!("  - Solvent molecule center of geometry must be:");
    eprintln!("    a) Inside the target box");
    eprintln!("    b) At least @thresh distance from any solute atom");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  # Cubic box with 1.4 nm minimum wall distance");
    eprintln!("  sim_box @topo system.top @pbc r @pos solute.cnf @solvent h2o.g96 @minwall 1.4");
    eprintln!();
    eprintln!("  # Use box size from solute file");
    eprintln!("  sim_box @topo system.top @pbc r @pos solute.cnf @solvent h2o.g96 @boxsize");
    eprintln!();
    eprintln!("  # Rectangular box with different wall distances");
    eprintln!("  sim_box @topo system.top @pbc r @pos solute.cnf @solvent h2o.g96 @minwall 1.0 1.5 2.0");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  GROMOS96 format to stdout containing solvated system");
}

#[derive(Debug)]
struct SimBoxArgs {
    topo: String,
    pbc: char,
    pos: String,
    solvent: String,
    minwall: Vec<f32>,
    thresh: f32,
    use_boxsize: bool,
}

impl Default for SimBoxArgs {
    fn default() -> Self {
        Self {
            topo: String::new(),
            pbc: 'r',
            pos: String::new(),
            solvent: String::new(),
            minwall: Vec::new(),
            thresh: 0.23, // Default from GROMOS++
            use_boxsize: false,
        }
    }
}

fn parse_args(args: Vec<String>) -> Result<SimBoxArgs, String> {
    let mut sb_args = SimBoxArgs::default();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing topology file for @topo".to_string());
                }
                sb_args.topo = args[i].clone();
            }
            "@pbc" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing boundary type for @pbc".to_string());
                }
                let pbc_str = args[i].clone();
                if pbc_str.len() != 1 {
                    return Err(format!("Invalid @pbc: {}", pbc_str));
                }
                sb_args.pbc = pbc_str.chars().next().unwrap();
                if sb_args.pbc != 'r' && sb_args.pbc != 't' {
                    return Err(format!("Invalid @pbc: {} (must be 'r' or 't')", sb_args.pbc));
                }
            }
            "@pos" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing solute file for @pos".to_string());
                }
                sb_args.pos = args[i].clone();
            }
            "@solvent" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing solvent file for @solvent".to_string());
                }
                sb_args.solvent = args[i].clone();
            }
            "@minwall" => {
                // Read 1 or 3 values
                loop {
                    i += 1;
                    if i >= args.len() || args[i].starts_with('@') {
                        i -= 1;
                        break;
                    }
                    let val: f32 = args[i]
                        .parse()
                        .map_err(|_| format!("Invalid @minwall value: {}", args[i]))?;
                    sb_args.minwall.push(val);
                }
                if sb_args.minwall.len() != 1 && sb_args.minwall.len() != 3 {
                    return Err(format!(
                        "@minwall requires 1 or 3 values, got {}",
                        sb_args.minwall.len()
                    ));
                }
            }
            "@thresh" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing threshold for @thresh".to_string());
                }
                sb_args.thresh = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid @thresh: {}", args[i]))?;
            }
            "@boxsize" => {
                sb_args.use_boxsize = true;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    // Validate
    if sb_args.topo.is_empty() {
        return Err("Missing required argument: @topo".to_string());
    }
    if sb_args.pos.is_empty() {
        return Err("Missing required argument: @pos".to_string());
    }
    if sb_args.solvent.is_empty() {
        return Err("Missing required argument: @solvent".to_string());
    }
    if !sb_args.use_boxsize && sb_args.minwall.is_empty() {
        return Err("Either @boxsize or @minwall must be specified".to_string());
    }
    if sb_args.use_boxsize && !sb_args.minwall.is_empty() {
        return Err("Cannot specify both @boxsize and @minwall".to_string());
    }

    Ok(sb_args)
}

/// Simple coordinate file reader that also parses box
fn read_coords_with_box(path: &str) -> Result<(Vec<Vec3>, Option<Vec3>), String> {
    let file = File::open(path)
        .map_err(|e| format!("Cannot open {}: {}", path, e))?;
    let reader = BufReader::new(file);

    let mut positions = Vec::new();
    let mut box_dims: Option<Vec3> = None;
    let mut in_position = false;
    let mut in_box = false;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Read error: {}", e))?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed == "TITLE" {
            continue;
        }

        if trimmed == "POSITION" || trimmed == "POSITIONRED" {
            in_position = true;
            in_box = false;
            continue;
        } else if trimmed == "BOX" {
            in_position = false;
            in_box = true;
            continue;
        } else if trimmed == "END" {
            in_position = false;
            in_box = false;
            continue;
        }

        if in_position && line.len() >= 24 {
            let coords = &line[24..];
            let parts: Vec<&str> = coords.split_whitespace().collect();
            if parts.len() >= 3 {
                let x: f32 = parts[0].parse()
                    .map_err(|_| format!("Invalid x coordinate: {}", parts[0]))?;
                let y: f32 = parts[1].parse()
                    .map_err(|_| format!("Invalid y coordinate: {}", parts[1]))?;
                let z: f32 = parts[2].parse()
                    .map_err(|_| format!("Invalid z coordinate: {}", parts[2]))?;
                positions.push(Vec3::new(x, y, z));
            }
        } else if in_box {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 3 {
                let lx: f32 = parts[0].parse()
                    .map_err(|_| format!("Invalid box x: {}", parts[0]))?;
                let ly: f32 = parts[1].parse()
                    .map_err(|_| format!("Invalid box y: {}", parts[1]))?;
                let lz: f32 = parts[2].parse()
                    .map_err(|_| format!("Invalid box z: {}", parts[2]))?;
                box_dims = Some(Vec3::new(lx, ly, lz));
            }
        }
    }

    if positions.is_empty() {
        return Err(format!("No positions found in {}", path));
    }

    Ok((positions, box_dims))
}

/// Calculate center of geometry
fn calc_cog(positions: &[Vec3]) -> Vec3 {
    let mut cog = Vec3::ZERO;
    for pos in positions {
        cog = cog + *pos;
    }
    cog / positions.len() as f32
}

/// Shift all positions by a vector
fn shift_positions(positions: &mut [Vec3], shift: Vec3) {
    for pos in positions {
        *pos = *pos + shift;
    }
}

/// Calculate maximum extent of solute in each dimension
fn calc_max_extent(positions: &[Vec3]) -> Vec3 {
    if positions.is_empty() {
        return Vec3::ZERO;
    }

    let mut min = positions[0];
    let mut max = positions[0];

    for pos in positions {
        min.x = min.x.min(pos.x);
        min.y = min.y.min(pos.y);
        min.z = min.z.min(pos.z);
        max.x = max.x.max(pos.x);
        max.y = max.y.max(pos.y);
        max.z = max.z.max(pos.z);
    }

    Vec3::new(max.x - min.x, max.y - min.y, max.z - min.z)
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args[1] == "-h" || args[1] == "--help" {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let sb_args = match parse_args(args) {
        Ok(args) => args,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    if sb_args.pbc == 't' {
        eprintln!("Error: Truncated octahedron (@pbc t) not yet implemented");
        eprintln!("       Use rectangular boxes (@pbc r) for now");
        process::exit(1);
    }

    eprintln!("sim_box - Solvating system");
    eprintln!("  Topology:   {}", sb_args.topo);
    eprintln!("  PBC:        {}", sb_args.pbc);
    eprintln!("  Solute:     {}", sb_args.pos);
    eprintln!("  Solvent:    {}", sb_args.solvent);
    eprintln!("  Threshold:  {} nm", sb_args.thresh);
    eprintln!();

    // Read solute coordinates
    eprintln!("Reading solute coordinates...");
    let (mut solute_pos, solute_box) = match read_coords_with_box(&sb_args.pos) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("Error reading solute: {}", e);
            process::exit(1);
        }
    };
    eprintln!("  Solute: {} atoms", solute_pos.len());

    // Read solvent coordinates
    eprintln!("Reading solvent coordinates...");
    let (solvent_pos_orig, solvent_box_opt) = match read_coords_with_box(&sb_args.solvent) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("Error reading solvent: {}", e);
            process::exit(1);
        }
    };

    let solvent_box = match solvent_box_opt {
        Some(b) => b,
        None => {
            eprintln!("Error: Solvent file must contain a BOX block");
            process::exit(1);
        }
    };
    eprintln!("  Solvent: {} atoms in box ({:.3}, {:.3}, {:.3}) nm",
              solvent_pos_orig.len(), solvent_box.x, solvent_box.y, solvent_box.z);

    // Determine target box size
    let target_box: Vec3 = if sb_args.use_boxsize {
        // Use box from solute file
        match solute_box {
            Some(b) => {
                eprintln!("Using box size from solute file: ({:.3}, {:.3}, {:.3}) nm", b.x, b.y, b.z);
                b
            }
            None => {
                eprintln!("Error: @boxsize specified but solute file has no BOX block");
                process::exit(1);
            }
        }
    } else {
        // Calculate from solute extent + minwall
        let extent = calc_max_extent(&solute_pos);
        eprintln!("  Solute extent: ({:.3}, {:.3}, {:.3}) nm", extent.x, extent.y, extent.z);

        let box_size = if sb_args.minwall.len() == 1 {
            // Cubic box
            let max_extent = extent.x.max(extent.y).max(extent.z);
            let size = max_extent + 2.0 * sb_args.minwall[0];
            eprintln!("  Creating cubic box: {:.3} nm (minwall: {:.3} nm)", size, sb_args.minwall[0]);
            Vec3::new(size, size, size)
        } else {
            // Rectangular box
            let size_x = extent.x + 2.0 * sb_args.minwall[0];
            let size_y = extent.y + 2.0 * sb_args.minwall[1];
            let size_z = extent.z + 2.0 * sb_args.minwall[2];
            eprintln!("  Creating rectangular box: ({:.3}, {:.3}, {:.3}) nm", size_x, size_y, size_z);
            eprintln!("  Minwall distances: ({:.3}, {:.3}, {:.3}) nm",
                      sb_args.minwall[0], sb_args.minwall[1], sb_args.minwall[2]);
            Vec3::new(size_x, size_y, size_z)
        };

        box_size
    };

    // Move solute to center
    let solute_cog = calc_cog(&solute_pos);
    let shift_to_origin = Vec3::ZERO - solute_cog;
    shift_positions(&mut solute_pos, shift_to_origin);
    eprintln!("  Shifted solute to origin");

    // Calculate how many solvent boxes needed in each dimension
    let nx = (target_box.x / solvent_box.x).ceil() as i32 + 1;
    let ny = (target_box.y / solvent_box.y).ceil() as i32 + 1;
    let nz = (target_box.z / solvent_box.z).ceil() as i32 + 1;

    eprintln!();
    eprintln!("Replicating solvent box...");
    eprintln!("  Need {} x {} x {} = {} boxes", nx, ny, nz, nx * ny * nz);

    // Replicate solvent boxes
    let mut all_solvent = Vec::new();
    let half_box = target_box * 0.5;

    // Center the replicated grid
    let start_x = -((nx as f32 - 1.0) * 0.5 * solvent_box.x);
    let start_y = -((ny as f32 - 1.0) * 0.5 * solvent_box.y);
    let start_z = -((nz as f32 - 1.0) * 0.5 * solvent_box.z);

    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let shift = Vec3::new(
                    start_x + ix as f32 * solvent_box.x,
                    start_y + iy as f32 * solvent_box.y,
                    start_z + iz as f32 * solvent_box.z,
                );
                for pos in &solvent_pos_orig {
                    all_solvent.push(*pos + shift);
                }
            }
        }
    }

    eprintln!("  Total solvent atoms before filtering: {}", all_solvent.len());

    // Assume solvent is water (3 atoms per molecule)
    // TODO: Get this from topology
    let atoms_per_molecule = 3;
    let num_molecules = all_solvent.len() / atoms_per_molecule;

    // Filter solvent molecules
    eprintln!("Filtering solvent molecules...");
    let mut kept_solvent = Vec::new();
    let thresh2 = sb_args.thresh * sb_args.thresh;

    for mol_idx in 0..num_molecules {
        let start_atom = mol_idx * atoms_per_molecule;
        let end_atom = start_atom + atoms_per_molecule;
        let mol_positions = &all_solvent[start_atom..end_atom];

        // Calculate center of geometry of this solvent molecule
        let mol_cog = calc_cog(mol_positions);

        // Check if inside target box (rectangular PBC)
        if mol_cog.x.abs() > half_box.x || mol_cog.y.abs() > half_box.y || mol_cog.z.abs() > half_box.z {
            continue; // Outside box
        }

        // Check minimum distance to solute atoms
        let mut min_dist2 = f32::MAX;
        for solute_atom in &solute_pos {
            let dx = mol_cog.x - solute_atom.x;
            let dy = mol_cog.y - solute_atom.y;
            let dz = mol_cog.z - solute_atom.z;
            let dist2 = dx * dx + dy * dy + dz * dz;
            if dist2 < min_dist2 {
                min_dist2 = dist2;
            }
        }

        // Keep if far enough from solute
        if min_dist2 > thresh2 {
            for atom_pos in mol_positions {
                kept_solvent.push(*atom_pos);
            }
        }
    }

    let kept_molecules = kept_solvent.len() / atoms_per_molecule;
    eprintln!("  Kept {} solvent molecules ({} atoms)", kept_molecules, kept_solvent.len());
    eprintln!("  Removed {} molecules due to clashes or being outside box",
              num_molecules - kept_molecules);

    // Combine solute and solvent
    let mut final_positions = solute_pos.clone();
    final_positions.extend(kept_solvent);

    eprintln!();
    eprintln!("Final system:");
    eprintln!("  Total atoms: {}", final_positions.len());
    eprintln!("  Solute atoms: {}", solute_pos.len());
    eprintln!("  Solvent atoms: {}", final_positions.len() - solute_pos.len());
    eprintln!("  Box: ({:.3}, {:.3}, {:.3}) nm", target_box.x, target_box.y, target_box.z);
    eprintln!();

    // Write to stdout
    let title = format!("Solvated system: {} in {}\nBox: {:.3} x {:.3} x {:.3} nm",
                        sb_args.pos, sb_args.solvent, target_box.x, target_box.y, target_box.z);

    if let Err(e) = write_g96(
        "/dev/stdout",
        &title,
        &final_positions,
        None,
        Some(target_box),
        None,
    ) {
        eprintln!("Error writing output: {}", e);
        process::exit(1);
    }

    eprintln!("Done!");
}
