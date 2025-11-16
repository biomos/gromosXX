/// REPLICA block parser for input files
///
/// Parses REPLICA blocks from .imd input files following GROMOS format:
///
/// ```
/// REPLICA
/// # NRET  LRESCALE  RET(1:NRET)
///   4     1         300.0 310.0 320.0 330.0
/// # TRIALS  EQUILIBRATE
///   1000    100
/// END
/// ```

use std::io::{BufRead, BufReader};
use std::fs::File;

#[derive(Debug, Clone)]
pub struct ReplicaBlock {
    /// Number of replicas
    pub num_replicas: usize,

    /// Whether to rescale velocities after temperature exchange (0 or 1)
    pub rescale_velocities: bool,

    /// Temperature values for each replica (K)
    pub temperatures: Vec<f64>,

    /// Number of MD steps between exchange attempts
    pub exchange_interval: usize,

    /// Number of equilibration steps before exchanges start
    pub equilibration_steps: usize,
}

impl Default for ReplicaBlock {
    fn default() -> Self {
        ReplicaBlock {
            num_replicas: 2,
            rescale_velocities: true,
            temperatures: vec![300.0, 310.0],
            exchange_interval: 1000,
            equilibration_steps: 100,
        }
    }
}

impl ReplicaBlock {
    /// Parse REPLICA block from input file
    pub fn parse(reader: &mut dyn BufRead) -> Result<Option<Self>, String> {
        let mut replica_block = None;
        let mut in_replica_block = false;
        let mut line_num = 0;
        let mut parsed_header = false;
        let mut parsed_trials = false;

        for line in reader.lines() {
            line_num += 1;
            let line = line.map_err(|e| format!("Error reading line {}: {}", line_num, e))?;
            let trimmed = line.trim();

            // Skip comments and empty lines
            if trimmed.starts_with('#') || trimmed.is_empty() {
                continue;
            }

            // Check for REPLICA block start
            if trimmed == "REPLICA" {
                in_replica_block = true;
                replica_block = Some(ReplicaBlock::default());
                continue;
            }

            // Check for block end
            if trimmed == "END" && in_replica_block {
                break;
            }

            // Parse REPLICA parameters
            if in_replica_block {
                let block = replica_block.as_mut().unwrap();
                let parts: Vec<&str> = trimmed.split_whitespace().collect();

                if parts.is_empty() {
                    continue;
                }

                // Line 1: NRET  LRESCALE  RET(1:NRET)
                if !parsed_header {
                    if parts.len() < 3 {
                        return Err(format!("Invalid REPLICA header at line {}: expected NRET LRESCALE TEMPS", line_num));
                    }

                    block.num_replicas = parts[0].parse::<usize>()
                        .map_err(|_| format!("Invalid NRET at line {}: {}", line_num, parts[0]))?;

                    let rescale = parts[1].parse::<i32>()
                        .map_err(|_| format!("Invalid LRESCALE at line {}: {}", line_num, parts[1]))?;

                    block.rescale_velocities = rescale == 1;

                    // Parse temperature values (rest of the line)
                    block.temperatures = parts[2..].iter()
                        .map(|s| s.parse::<f64>())
                        .collect::<Result<Vec<f64>, _>>()
                        .map_err(|_| format!("Invalid temperature values at line {}", line_num))?;

                    // Validate temperature count
                    if block.temperatures.len() != block.num_replicas {
                        return Err(format!(
                            "Expected {} temperatures for {} replicas, got {} at line {}",
                            block.num_replicas, block.num_replicas, block.temperatures.len(), line_num
                        ));
                    }

                    // Validate temperatures are positive and increasing
                    for (i, temp) in block.temperatures.iter().enumerate() {
                        if *temp <= 0.0 {
                            return Err(format!("Temperature {} must be positive, got {} at line {}", i, temp, line_num));
                        }
                        if i > 0 && block.temperatures[i] <= block.temperatures[i - 1] {
                            return Err(format!(
                                "Temperatures must be strictly increasing, but T[{}]={} <= T[{}]={} at line {}",
                                i, block.temperatures[i], i - 1, block.temperatures[i - 1], line_num
                            ));
                        }
                    }

                    parsed_header = true;
                }
                // Line 2: TRIALS  EQUILIBRATE
                else if !parsed_trials {
                    if parts.len() != 2 {
                        return Err(format!("Invalid TRIALS line at {}: expected TRIALS EQUILIBRATE", line_num));
                    }

                    block.exchange_interval = parts[0].parse::<usize>()
                        .map_err(|_| format!("Invalid TRIALS at line {}: {}", line_num, parts[0]))?;

                    block.equilibration_steps = parts[1].parse::<usize>()
                        .map_err(|_| format!("Invalid EQUILIBRATE at line {}: {}", line_num, parts[1]))?;

                    if block.exchange_interval == 0 {
                        return Err(format!("TRIALS must be greater than 0 at line {}", line_num));
                    }

                    parsed_trials = true;
                }
            }
        }

        Ok(replica_block)
    }

    /// Parse from file path
    pub fn parse_file(path: &str) -> Result<Option<Self>, String> {
        let file = File::open(path)
            .map_err(|e| format!("Could not open file {}: {}", path, e))?;
        let mut reader = BufReader::new(file);
        Self::parse(&mut reader)
    }

    /// Write REPLICA block to string
    pub fn to_string(&self) -> String {
        let mut output = String::new();
        output.push_str("REPLICA\n");
        output.push_str("# NRET  LRESCALE  RET(1:NRET)\n");

        let rescale_val = if self.rescale_velocities { 1 } else { 0 };

        output.push_str(&format!("  {}     {}         ", self.num_replicas, rescale_val));
        for (i, temp) in self.temperatures.iter().enumerate() {
            if i > 0 {
                output.push_str(" ");
            }
            output.push_str(&format!("{:.1}", temp));
        }
        output.push('\n');

        output.push_str("# TRIALS  EQUILIBRATE\n");
        output.push_str(&format!("  {}    {}\n", self.exchange_interval, self.equilibration_steps));
        output.push_str("END\n");

        output
    }

    /// Validate the replica block configuration
    pub fn validate(&self) -> Result<(), String> {
        if self.num_replicas < 2 {
            return Err("At least 2 replicas are required for REMD".to_string());
        }

        if self.temperatures.len() != self.num_replicas {
            return Err(format!(
                "Number of temperatures ({}) does not match number of replicas ({})",
                self.temperatures.len(), self.num_replicas
            ));
        }

        if self.exchange_interval == 0 {
            return Err("Exchange interval (TRIALS) must be greater than 0".to_string());
        }

        // Check that temperatures are strictly increasing
        for i in 1..self.temperatures.len() {
            if self.temperatures[i] <= self.temperatures[i - 1] {
                return Err(format!(
                    "Temperatures must be strictly increasing, but T[{}]={} <= T[{}]={}",
                    i, self.temperatures[i], i - 1, self.temperatures[i - 1]
                ));
            }
        }

        Ok(())
    }

    /// Get temperature spacing statistics (useful for optimal setup)
    pub fn temperature_spacing_stats(&self) -> (f64, f64, f64) {
        if self.temperatures.len() < 2 {
            return (0.0, 0.0, 0.0);
        }

        let mut spacings = Vec::new();
        for i in 1..self.temperatures.len() {
            spacings.push(self.temperatures[i] - self.temperatures[i - 1]);
        }

        let min_spacing = spacings.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_spacing = spacings.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let avg_spacing = spacings.iter().sum::<f64>() / spacings.len() as f64;

        (min_spacing, max_spacing, avg_spacing)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_replica_block() {
        let input = r"
REPLICA
# NRET  LRESCALE  RET(1:NRET)
  4     1         300.0 310.0 320.0 330.0
# TRIALS  EQUILIBRATE
  1000    100
END
";
        let mut cursor = Cursor::new(input);
        let block = ReplicaBlock::parse(&mut cursor).unwrap().unwrap();

        assert_eq!(block.num_replicas, 4);
        assert!(block.rescale_velocities);
        assert_eq!(block.temperatures, vec![300.0, 310.0, 320.0, 330.0]);
        assert_eq!(block.exchange_interval, 1000);
        assert_eq!(block.equilibration_steps, 100);
    }

    #[test]
    fn test_parse_replica_block_no_rescale() {
        let input = r"
REPLICA
# NRET  LRESCALE  RET(1:NRET)
  3     0         298.15 308.15 318.15
# TRIALS  EQUILIBRATE
  500    50
END
";
        let mut cursor = Cursor::new(input);
        let block = ReplicaBlock::parse(&mut cursor).unwrap().unwrap();

        assert_eq!(block.num_replicas, 3);
        assert!(!block.rescale_velocities);
        assert_eq!(block.temperatures.len(), 3);
        assert_eq!(block.exchange_interval, 500);
        assert_eq!(block.equilibration_steps, 50);
    }

    #[test]
    fn test_invalid_temperature_count() {
        let input = r"
REPLICA
# NRET  LRESCALE  RET(1:NRET)
  4     1         300.0 310.0 320.0
# TRIALS  EQUILIBRATE
  1000    100
END
";
        let mut cursor = Cursor::new(input);
        let result = ReplicaBlock::parse(&mut cursor);

        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Expected 4 temperatures"));
    }

    #[test]
    fn test_invalid_temperature_order() {
        let input = r"
REPLICA
# NRET  LRESCALE  RET(1:NRET)
  3     1         300.0 320.0 310.0
# TRIALS  EQUILIBRATE
  1000    100
END
";
        let mut cursor = Cursor::new(input);
        let result = ReplicaBlock::parse(&mut cursor);

        assert!(result.is_err());
        assert!(result.unwrap_err().contains("strictly increasing"));
    }

    #[test]
    fn test_negative_temperature() {
        let input = r"
REPLICA
# NRET  LRESCALE  RET(1:NRET)
  2     1         -300.0 310.0
# TRIALS  EQUILIBRATE
  1000    100
END
";
        let mut cursor = Cursor::new(input);
        let result = ReplicaBlock::parse(&mut cursor);

        assert!(result.is_err());
        assert!(result.unwrap_err().contains("must be positive"));
    }

    #[test]
    fn test_zero_exchange_interval() {
        let input = r"
REPLICA
# NRET  LRESCALE  RET(1:NRET)
  2     1         300.0 310.0
# TRIALS  EQUILIBRATE
  0    100
END
";
        let mut cursor = Cursor::new(input);
        let result = ReplicaBlock::parse(&mut cursor);

        assert!(result.is_err());
        assert!(result.unwrap_err().contains("TRIALS must be greater than 0"));
    }

    #[test]
    fn test_roundtrip() {
        let block = ReplicaBlock {
            num_replicas: 4,
            rescale_velocities: true,
            temperatures: vec![300.0, 310.0, 320.0, 330.0],
            exchange_interval: 1000,
            equilibration_steps: 100,
        };

        let output = block.to_string();
        let mut cursor = Cursor::new(output.as_bytes());
        let parsed = ReplicaBlock::parse(&mut cursor).unwrap().unwrap();

        assert_eq!(block.num_replicas, parsed.num_replicas);
        assert_eq!(block.rescale_velocities, parsed.rescale_velocities);
        assert_eq!(block.temperatures, parsed.temperatures);
        assert_eq!(block.exchange_interval, parsed.exchange_interval);
        assert_eq!(block.equilibration_steps, parsed.equilibration_steps);
    }

    #[test]
    fn test_validate() {
        let mut block = ReplicaBlock {
            num_replicas: 4,
            rescale_velocities: true,
            temperatures: vec![300.0, 310.0, 320.0, 330.0],
            exchange_interval: 1000,
            equilibration_steps: 100,
        };

        assert!(block.validate().is_ok());

        // Test with only 1 replica
        block.num_replicas = 1;
        block.temperatures = vec![300.0];
        assert!(block.validate().is_err());

        // Test with mismatched temperature count
        block.num_replicas = 4;
        block.temperatures = vec![300.0, 310.0];
        assert!(block.validate().is_err());

        // Test with zero exchange interval
        block.temperatures = vec![300.0, 310.0, 320.0, 330.0];
        block.exchange_interval = 0;
        assert!(block.validate().is_err());
    }

    #[test]
    fn test_temperature_spacing_stats() {
        let block = ReplicaBlock {
            num_replicas: 4,
            rescale_velocities: true,
            temperatures: vec![300.0, 310.0, 320.0, 330.0],
            exchange_interval: 1000,
            equilibration_steps: 100,
        };

        let (min, max, avg) = block.temperature_spacing_stats();
        assert_eq!(min, 10.0);
        assert_eq!(max, 10.0);
        assert_eq!(avg, 10.0);
    }

    #[test]
    fn test_temperature_spacing_stats_variable() {
        let block = ReplicaBlock {
            num_replicas: 4,
            rescale_velocities: true,
            temperatures: vec![300.0, 305.0, 315.0, 330.0],
            exchange_interval: 1000,
            equilibration_steps: 100,
        };

        let (min, max, avg) = block.temperature_spacing_stats();
        assert_eq!(min, 5.0);
        assert_eq!(max, 15.0);
        assert_eq!(avg, 10.0); // (5 + 10 + 15) / 3 = 10
    }
}
