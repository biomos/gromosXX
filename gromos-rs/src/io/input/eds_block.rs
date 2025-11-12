/// EDS block parser for input files
///
/// Parses EDS blocks from .imd input files following GROMOS format:
///
/// ```
/// EDS
/// # NUMSTATES  FORM  S (or multiple S values for multi-s)
///   3          1     0.5
/// # E_OFFSETS (E_i^R for each state)
///   0.0  10.0  20.0
/// # TEMP
///   300.0
/// # SEARCH  E_MAX  E_MIN  (optional, for AEDS)
///   1       10.0   -50.0
/// END
/// ```

use crate::eds::{EDSForm, EDSParameters, AEDSParameters};
use std::io::{BufRead, BufReader};
use std::fs::File;

#[derive(Debug, Clone)]
pub struct EdsBlock {
    pub num_states: usize,
    pub form: EDSForm,
    pub s_values: Vec<f64>,
    pub e_offsets: Vec<f64>,
    pub temperature: f64,
    pub search_enabled: bool,
    pub e_max: f64,
    pub e_min: f64,
}

impl Default for EdsBlock {
    fn default() -> Self {
        EdsBlock {
            num_states: 2,
            form: EDSForm::SingleS,
            s_values: vec![0.5],
            e_offsets: vec![0.0, 10.0],
            temperature: 300.0,
            search_enabled: false,
            e_max: 10.0,
            e_min: -50.0,
        }
    }
}

impl EdsBlock {
    /// Parse EDS block from input file
    pub fn parse(reader: &mut dyn BufRead) -> Result<Option<Self>, String> {
        let mut eds_block = None;
        let mut in_eds_block = false;
        let mut line_num = 0;
        let mut parsed_header = false;
        let mut parsed_offsets = false;
        let mut parsed_temp = false;

        for line in reader.lines() {
            line_num += 1;
            let line = line.map_err(|e| format!("Error reading line {}: {}", line_num, e))?;
            let trimmed = line.trim();

            // Skip comments and empty lines
            if trimmed.starts_with('#') || trimmed.is_empty() {
                continue;
            }

            // Check for EDS block start
            if trimmed == "EDS" {
                in_eds_block = true;
                eds_block = Some(EdsBlock::default());
                continue;
            }

            // Check for block end
            if trimmed == "END" && in_eds_block {
                break;
            }

            // Parse EDS parameters
            if in_eds_block {
                let block = eds_block.as_mut().unwrap();
                let parts: Vec<&str> = trimmed.split_whitespace().collect();

                if parts.is_empty() {
                    continue;
                }

                // Line 1: NUMSTATES  FORM  S (and possibly more S values)
                if !parsed_header {
                    if parts.len() < 3 {
                        return Err(format!("Invalid EDS header at line {}: expected NUMSTATES FORM S", line_num));
                    }

                    block.num_states = parts[0].parse::<usize>()
                        .map_err(|_| format!("Invalid NUMSTATES at line {}: {}", line_num, parts[0]))?;

                    let form_val = parts[1].parse::<i32>()
                        .map_err(|_| format!("Invalid FORM at line {}: {}", line_num, parts[1]))?;

                    block.form = match form_val {
                        1 => EDSForm::SingleS,
                        2 => EDSForm::MultiS,
                        3 => EDSForm::PairS,
                        _ => return Err(format!("Invalid FORM value at line {}: {}", line_num, form_val)),
                    };

                    // Parse S values (rest of the line)
                    block.s_values = parts[2..].iter()
                        .map(|s| s.parse::<f64>())
                        .collect::<Result<Vec<f64>, _>>()
                        .map_err(|_| format!("Invalid S values at line {}", line_num))?;

                    // Validate S values count based on form
                    match block.form {
                        EDSForm::SingleS => {
                            if block.s_values.len() != 1 {
                                return Err(format!("SingleS form requires exactly 1 s value, got {} at line {}", block.s_values.len(), line_num));
                            }
                        }
                        EDSForm::MultiS => {
                            let expected = (block.num_states * (block.num_states - 1)) / 2;
                            if block.s_values.len() != expected {
                                return Err(format!("MultiS form requires {} s values for {} states, got {} at line {}",
                                    expected, block.num_states, block.s_values.len(), line_num));
                            }
                        }
                        EDSForm::PairS => {
                            // Variable number allowed for PairS
                        }
                    }

                    parsed_header = true;
                }
                // Line 2: E_OFFSETS (all on one line)
                else if !parsed_offsets {
                    block.e_offsets = parts.iter()
                        .map(|s| s.parse::<f64>())
                        .collect::<Result<Vec<f64>, _>>()
                        .map_err(|_| format!("Invalid E_OFFSETS at line {}", line_num))?;

                    if block.e_offsets.len() != block.num_states {
                        return Err(format!("Expected {} energy offsets for {} states, got {} at line {}",
                            block.num_states, block.num_states, block.e_offsets.len(), line_num));
                    }

                    parsed_offsets = true;
                }
                // Line 3: TEMP
                else if !parsed_temp {
                    if parts.len() != 1 {
                        return Err(format!("Invalid TEMP line at {}: expected single value", line_num));
                    }

                    block.temperature = parts[0].parse::<f64>()
                        .map_err(|_| format!("Invalid TEMP at line {}: {}", line_num, parts[0]))?;

                    parsed_temp = true;
                }
                // Line 4 (optional): SEARCH  E_MAX  E_MIN
                else if parts.len() == 3 {
                    let search = parts[0].parse::<i32>()
                        .map_err(|_| format!("Invalid SEARCH at line {}: {}", line_num, parts[0]))?;

                    block.search_enabled = search == 1;

                    block.e_max = parts[1].parse::<f64>()
                        .map_err(|_| format!("Invalid E_MAX at line {}: {}", line_num, parts[1]))?;

                    block.e_min = parts[2].parse::<f64>()
                        .map_err(|_| format!("Invalid E_MIN at line {}: {}", line_num, parts[2]))?;
                }
            }
        }

        Ok(eds_block)
    }

    /// Parse from file path
    pub fn parse_file(path: &str) -> Result<Option<Self>, String> {
        let file = File::open(path)
            .map_err(|e| format!("Could not open file {}: {}", path, e))?;
        let mut reader = BufReader::new(file);
        Self::parse(&mut reader)
    }

    /// Convert to EDSParameters
    pub fn to_parameters(&self, num_atoms: usize) -> Result<EDSParameters, String> {
        EDSParameters::new(
            self.form,
            self.s_values.clone(),
            self.e_offsets.clone(),
            self.temperature,
            num_atoms,
        )
    }

    /// Convert to AEDSParameters
    pub fn to_aeds_parameters(&self, num_atoms: usize) -> Result<AEDSParameters, String> {
        let eds = self.to_parameters(num_atoms)?;
        Ok(AEDSParameters::new(
            eds,
            self.e_max,
            self.e_min,
            self.search_enabled,
        ))
    }

    /// Write EDS block to string
    pub fn to_string(&self) -> String {
        let mut output = String::new();
        output.push_str("EDS\n");

        // Header line
        output.push_str("# NUMSTATES  FORM  S");
        if self.form != EDSForm::SingleS {
            output.push_str(" (values)");
        }
        output.push('\n');

        let form_val = match self.form {
            EDSForm::SingleS => 1,
            EDSForm::MultiS => 2,
            EDSForm::PairS => 3,
        };

        output.push_str(&format!("  {}          {}     ", self.num_states, form_val));
        for (i, s) in self.s_values.iter().enumerate() {
            if i > 0 {
                output.push_str("  ");
            }
            output.push_str(&format!("{:.4}", s));
        }
        output.push('\n');

        // Energy offsets
        output.push_str("# E_OFFSETS (E_i^R for each state)\n");
        output.push_str("  ");
        for (i, offset) in self.e_offsets.iter().enumerate() {
            if i > 0 {
                output.push_str("  ");
            }
            output.push_str(&format!("{:.2}", offset));
        }
        output.push('\n');

        // Temperature
        output.push_str("# TEMP\n");
        output.push_str(&format!("  {:.1}\n", self.temperature));

        // Optional AEDS parameters
        if self.search_enabled {
            output.push_str("# SEARCH  E_MAX  E_MIN\n");
            output.push_str(&format!("  1       {:.2}   {:.2}\n", self.e_max, self.e_min));
        }

        output.push_str("END\n");

        output
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_eds_block_single_s() {
        let input = r"
EDS
# NUMSTATES  FORM  S
  3          1     0.5
# E_OFFSETS
  0.0  10.0  20.0
# TEMP
  300.0
END
";
        let mut cursor = Cursor::new(input);
        let block = EdsBlock::parse(&mut cursor).unwrap().unwrap();

        assert_eq!(block.num_states, 3);
        assert_eq!(block.form, EDSForm::SingleS);
        assert_eq!(block.s_values, vec![0.5]);
        assert_eq!(block.e_offsets, vec![0.0, 10.0, 20.0]);
        assert_eq!(block.temperature, 300.0);
        assert!(!block.search_enabled);
    }

    #[test]
    fn test_parse_eds_block_with_aeds() {
        let input = r"
EDS
# NUMSTATES  FORM  S
  2          1     0.5
# E_OFFSETS
  0.0  10.0
# TEMP
  300.0
# SEARCH  E_MAX  E_MIN
  1       10.0   -50.0
END
";
        let mut cursor = Cursor::new(input);
        let block = EdsBlock::parse(&mut cursor).unwrap().unwrap();

        assert_eq!(block.num_states, 2);
        assert!(block.search_enabled);
        assert_eq!(block.e_max, 10.0);
        assert_eq!(block.e_min, -50.0);
    }

    #[test]
    fn test_parse_eds_block_multi_s() {
        let input = r"
EDS
# NUMSTATES  FORM  S (3 values for 3 states = 3 pairs)
  3          2     0.4  0.5  0.6
# E_OFFSETS
  0.0  10.0  20.0
# TEMP
  300.0
END
";
        let mut cursor = Cursor::new(input);
        let block = EdsBlock::parse(&mut cursor).unwrap().unwrap();

        assert_eq!(block.num_states, 3);
        assert_eq!(block.form, EDSForm::MultiS);
        assert_eq!(block.s_values.len(), 3); // N(N-1)/2 = 3
    }

    #[test]
    fn test_invalid_offset_count() {
        let input = r"
EDS
# NUMSTATES  FORM  S
  3          1     0.5
# E_OFFSETS
  0.0  10.0
# TEMP
  300.0
END
";
        let mut cursor = Cursor::new(input);
        let result = EdsBlock::parse(&mut cursor);

        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Expected 3 energy offsets"));
    }

    #[test]
    fn test_roundtrip() {
        let block = EdsBlock {
            num_states: 3,
            form: EDSForm::SingleS,
            s_values: vec![0.5],
            e_offsets: vec![0.0, 10.0, 20.0],
            temperature: 300.0,
            search_enabled: true,
            e_max: 10.0,
            e_min: -50.0,
        };

        let output = block.to_string();
        let mut cursor = Cursor::new(output.as_bytes());
        let parsed = EdsBlock::parse(&mut cursor).unwrap().unwrap();

        assert_eq!(block.num_states, parsed.num_states);
        assert_eq!(block.form, parsed.form);
        assert_eq!(block.s_values, parsed.s_values);
        assert_eq!(block.e_offsets, parsed.e_offsets);
        assert_eq!(block.temperature, parsed.temperature);
        assert_eq!(block.search_enabled, parsed.search_enabled);
    }
}
