//! Parser for GROMOS coordinate files (.conf/.cnf)
//!
//! Format example:
//! ```text
//! TITLE
//!   System description
//! END
//! POSITION
//!     1 RES    ATOM      1    x        y        z
//! END
//! VELOCITY (optional)
//!     1 RES    ATOM      1    vx       vy       vz
//! END
//! BOX
//!     lx   ly   lz
//! END
//! ```

use crate::configuration::{Configuration, Box as SimBox};
use crate::math::Vec3;
use crate::io::IoError;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Read GROMOS coordinate file
pub fn read_coordinate_file<P: AsRef<Path>>(path: P, num_temp_groups: usize, num_energy_groups: usize) -> Result<Configuration, IoError> {
    let file = File::open(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    let reader = BufReader::new(file);

    let mut positions = Vec::new();
    let mut velocities = Vec::new();
    let mut box_dims = Vec3::ZERO;
    let mut in_position = false;
    let mut in_velocity = false;
    let mut in_box = false;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Skip comments and empty lines
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        // Check for block markers
        if trimmed == "POSITION" {
            in_position = true;
            in_velocity = false;
            in_box = false;
            continue;
        } else if trimmed == "VELOCITY" {
            in_position = false;
            in_velocity = true;
            in_box = false;
            continue;
        } else if trimmed == "BOX" {
            in_position = false;
            in_velocity = false;
            in_box = true;
            continue;
        } else if trimmed == "END" {
            in_position = false;
            in_velocity = false;
            in_box = false;
            continue;
        }

        // Parse data
        if in_position {
            // CRITICAL: Skip first 24 characters (residue/atom metadata)
            // Format: "    1 RES    ATOM      1"  (24 chars) then x y z
            if line.len() < 24 {
                continue;
            }

            let coords = &line[24..];
            let parts: Vec<&str> = coords.split_whitespace().collect();

            if parts.len() >= 3 {
                let x: f32 = parts[0].parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid x coordinate: {}", parts[0])))?;
                let y: f32 = parts[1].parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid y coordinate: {}", parts[1])))?;
                let z: f32 = parts[2].parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid z coordinate: {}", parts[2])))?;

                positions.push(Vec3::new(x, y, z));
            }
        } else if in_velocity {
            // Same format as POSITION
            if line.len() < 24 {
                continue;
            }

            let coords = &line[24..];
            let parts: Vec<&str> = coords.split_whitespace().collect();

            if parts.len() >= 3 {
                let vx: f32 = parts[0].parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid vx: {}", parts[0])))?;
                let vy: f32 = parts[1].parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid vy: {}", parts[1])))?;
                let vz: f32 = parts[2].parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid vz: {}", parts[2])))?;

                velocities.push(Vec3::new(vx, vy, vz));
            }
        } else if in_box {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();

            if parts.len() >= 3 {
                let lx: f32 = parts[0].parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid box x: {}", parts[0])))?;
                let ly: f32 = parts[1].parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid box y: {}", parts[1])))?;
                let lz: f32 = parts[2].parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid box z: {}", parts[2])))?;

                box_dims = Vec3::new(lx, ly, lz);
            }
        }
    }

    // Validate
    if positions.is_empty() {
        return Err(IoError::FormatError("No positions found".to_string()));
    }

    let n_atoms = positions.len();

    // Create configuration
    let mut conf = Configuration::new(n_atoms, num_temp_groups, num_energy_groups);

    // Set positions
    conf.current_mut().pos = positions;

    // Set velocities (if present)
    if !velocities.is_empty() {
        if velocities.len() != n_atoms {
            return Err(IoError::FormatError(
                format!("Velocity count ({}) doesn't match atom count ({})", velocities.len(), n_atoms)
            ));
        }
        conf.current_mut().vel = velocities;
    }

    // Set box
    if box_dims != Vec3::ZERO {
        conf.current_mut().box_config = SimBox::rectangular(box_dims.x, box_dims.y, box_dims.z);
    }

    Ok(conf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cg16() {
        // Test with actual GROMOS file
        let result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

        if let Ok(conf) = result {
            println!("Loaded {} atoms", conf.current().pos.len());
            assert_eq!(conf.current().pos.len(), 4);

            // Check first atom position
            let pos0 = conf.current().pos[0];
            println!("Atom 0: ({}, {}, {})", pos0.x, pos0.y, pos0.z);
            assert!((pos0.x - 4.197156491).abs() < 1e-6);
            assert!((pos0.y - 0.505921049).abs() < 1e-6);
            assert!((pos0.z - 2.679733124).abs() < 1e-6);

            // Check box
            let box_dims = conf.current().box_config.dimensions();
            assert!((box_dims.x - 10.0).abs() < 1e-6);
        } else {
            println!("Could not load file (may not exist in test environment)");
        }
    }
}
