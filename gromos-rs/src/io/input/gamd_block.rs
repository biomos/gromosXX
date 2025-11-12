/// GAMD block parser for input files
///
/// Parses GAMD blocks from .imd input files following GROMOS format:
///
/// ```
/// GAMD
/// # SEARCH  FORM  THRESH
///   1       2     1          # CMD search, total boost, lower bound
/// # SIGMA0_DIH  SIGMA0_TOT
///   6.0         6.0
/// # K_DIH    K_TOT      E_DIH    E_TOT (for production mode)
///   0.0      0.001234   0.0      -1234.56
/// # EQUIL  WINDOW
///   1000   0
/// END
/// ```

use crate::gamd::{GamdParameters, SearchMode, BoostForm, ThresholdType};
use std::io::{BufRead, BufReader};
use std::fs::File;

#[derive(Debug, Clone)]
pub struct GamdBlock {
    pub search_mode: SearchMode,
    pub boost_form: BoostForm,
    pub threshold_type: ThresholdType,
    pub sigma0_dih: f64,
    pub sigma0_tot: f64,
    pub k_dih: Option<f64>,
    pub k_tot: Option<f64>,
    pub e_dih: Option<f64>,
    pub e_tot: Option<f64>,
    pub equilibration_steps: usize,
    pub window_size: usize,
}

impl Default for GamdBlock {
    fn default() -> Self {
        GamdBlock {
            search_mode: SearchMode::CmdSearch,
            boost_form: BoostForm::TotalBoost,
            threshold_type: ThresholdType::LowerBound,
            sigma0_dih: 6.0,
            sigma0_tot: 6.0,
            k_dih: None,
            k_tot: None,
            e_dih: None,
            e_tot: None,
            equilibration_steps: 0,
            window_size: 0,
        }
    }
}

impl GamdBlock {
    /// Parse GAMD block from input file
    pub fn parse(reader: &mut dyn BufRead) -> Result<Option<Self>, String> {
        let mut gamd_block = None;
        let mut in_gamd_block = false;
        let mut line_num = 0;

        for line in reader.lines() {
            line_num += 1;
            let line = line.map_err(|e| format!("Error reading line {}: {}", line_num, e))?;
            let trimmed = line.trim();

            // Skip comments and empty lines
            if trimmed.starts_with('#') || trimmed.is_empty() {
                continue;
            }

            // Check for GAMD block start
            if trimmed == "GAMD" {
                in_gamd_block = true;
                gamd_block = Some(GamdBlock::default());
                continue;
            }

            // Check for block end
            if trimmed == "END" && in_gamd_block {
                break;
            }

            // Parse GAMD parameters
            if in_gamd_block {
                let block = gamd_block.as_mut().unwrap();
                let parts: Vec<&str> = trimmed.split_whitespace().collect();

                if parts.is_empty() {
                    continue;
                }

                // Try to parse as different parameter types
                // Line 1: SEARCH  FORM  THRESH
                if parts.len() == 3 && block.search_mode == SearchMode::CmdSearch {
                    let search = parts[0].parse::<i32>()
                        .map_err(|_| format!("Invalid SEARCH value at line {}: {}", line_num, parts[0]))?;
                    let form = parts[1].parse::<i32>()
                        .map_err(|_| format!("Invalid FORM value at line {}: {}", line_num, parts[1]))?;
                    let thresh = parts[2].parse::<i32>()
                        .map_err(|_| format!("Invalid THRESH value at line {}: {}", line_num, parts[2]))?;

                    block.search_mode = match search {
                        0 => SearchMode::NoSearch,
                        1 => SearchMode::CmdSearch,
                        2 => SearchMode::GamdSearch,
                        _ => return Err(format!("Invalid SEARCH value at line {}: {}", line_num, search)),
                    };

                    block.boost_form = match form {
                        1 => BoostForm::DihedralBoost,
                        2 => BoostForm::TotalBoost,
                        3 => BoostForm::DualBoost,
                        _ => return Err(format!("Invalid FORM value at line {}: {}", line_num, form)),
                    };

                    block.threshold_type = match thresh {
                        1 => ThresholdType::LowerBound,
                        2 => ThresholdType::UpperBound,
                        _ => return Err(format!("Invalid THRESH value at line {}: {}", line_num, thresh)),
                    };
                }
                // Line 2: SIGMA0_DIH  SIGMA0_TOT
                else if parts.len() == 2 && block.sigma0_dih == 6.0 {
                    block.sigma0_dih = parts[0].parse::<f64>()
                        .map_err(|_| format!("Invalid SIGMA0_DIH at line {}: {}", line_num, parts[0]))?;
                    block.sigma0_tot = parts[1].parse::<f64>()
                        .map_err(|_| format!("Invalid SIGMA0_TOT at line {}: {}", line_num, parts[1]))?;
                }
                // Line 3 (optional): K_DIH  K_TOT  E_DIH  E_TOT
                else if parts.len() == 4 && block.k_dih.is_none() {
                    block.k_dih = Some(parts[0].parse::<f64>()
                        .map_err(|_| format!("Invalid K_DIH at line {}: {}", line_num, parts[0]))?);
                    block.k_tot = Some(parts[1].parse::<f64>()
                        .map_err(|_| format!("Invalid K_TOT at line {}: {}", line_num, parts[1]))?);
                    block.e_dih = Some(parts[2].parse::<f64>()
                        .map_err(|_| format!("Invalid E_DIH at line {}: {}", line_num, parts[2]))?);
                    block.e_tot = Some(parts[3].parse::<f64>()
                        .map_err(|_| format!("Invalid E_TOT at line {}: {}", line_num, parts[3]))?);
                }
                // Line 4 or 3: EQUIL  WINDOW
                else if parts.len() == 2 {
                    block.equilibration_steps = parts[0].parse::<usize>()
                        .map_err(|_| format!("Invalid EQUIL at line {}: {}", line_num, parts[0]))?;
                    block.window_size = parts[1].parse::<usize>()
                        .map_err(|_| format!("Invalid WINDOW at line {}: {}", line_num, parts[1]))?;
                }
            }
        }

        Ok(gamd_block)
    }

    /// Parse from file path
    pub fn parse_file(path: &str) -> Result<Option<Self>, String> {
        let file = File::open(path)
            .map_err(|e| format!("Could not open file {}: {}", path, e))?;
        let mut reader = BufReader::new(file);
        Self::parse(&mut reader)
    }

    /// Convert to GamdParameters
    pub fn to_parameters(&self) -> GamdParameters {
        let mut params = GamdParameters::new(
            self.search_mode,
            self.boost_form,
            self.threshold_type,
            self.sigma0_dih,
            self.sigma0_tot,
        );

        // Set fixed parameters if in production mode
        if self.search_mode == SearchMode::NoSearch {
            if let Some(k) = self.k_dih {
                params.k_dih = k;
            }
            if let Some(k) = self.k_tot {
                params.k_tot = k;
            }
            if let Some(e) = self.e_dih {
                params.e_dih = e;
            }
            if let Some(e) = self.e_tot {
                params.e_tot = e;
            }
        }

        params.equilibration_steps = self.equilibration_steps;
        params.window_size = self.window_size;

        params
    }

    /// Write GAMD block to string
    pub fn to_string(&self) -> String {
        let mut output = String::new();
        output.push_str("GAMD\n");
        output.push_str("# SEARCH  FORM  THRESH\n");

        let search_val = match self.search_mode {
            SearchMode::NoSearch => 0,
            SearchMode::CmdSearch => 1,
            SearchMode::GamdSearch => 2,
        };

        let form_val = match self.boost_form {
            BoostForm::DihedralBoost => 1,
            BoostForm::TotalBoost => 2,
            BoostForm::DualBoost => 3,
        };

        let thresh_val = match self.threshold_type {
            ThresholdType::LowerBound => 1,
            ThresholdType::UpperBound => 2,
        };

        output.push_str(&format!("  {}       {}     {}\n", search_val, form_val, thresh_val));
        output.push_str("# SIGMA0_DIH  SIGMA0_TOT\n");
        output.push_str(&format!("  {:.1}         {:.1}\n", self.sigma0_dih, self.sigma0_tot));

        if let (Some(k_dih), Some(k_tot), Some(e_dih), Some(e_tot)) =
            (self.k_dih, self.k_tot, self.e_dih, self.e_tot)
        {
            output.push_str("# K_DIH    K_TOT      E_DIH    E_TOT\n");
            output.push_str(&format!(
                "  {:.6}  {:.6}   {:.2}     {:.2}\n",
                k_dih, k_tot, e_dih, e_tot
            ));
        }

        output.push_str("# EQUIL  WINDOW\n");
        output.push_str(&format!("  {}   {}\n", self.equilibration_steps, self.window_size));
        output.push_str("END\n");

        output
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_gamd_block_cmd_search() {
        let input = r"
GAMD
# SEARCH  FORM  THRESH
  1       2     1
# SIGMA0_DIH  SIGMA0_TOT
  6.0         6.0
# EQUIL  WINDOW
  1000   0
END
";
        let mut cursor = Cursor::new(input);
        let block = GamdBlock::parse(&mut cursor).unwrap().unwrap();

        assert_eq!(block.search_mode, SearchMode::CmdSearch);
        assert_eq!(block.boost_form, BoostForm::TotalBoost);
        assert_eq!(block.threshold_type, ThresholdType::LowerBound);
        assert_eq!(block.sigma0_dih, 6.0);
        assert_eq!(block.sigma0_tot, 6.0);
        assert_eq!(block.equilibration_steps, 1000);
        assert_eq!(block.window_size, 0);
    }

    #[test]
    fn test_parse_gamd_block_production() {
        let input = r"
GAMD
# SEARCH  FORM  THRESH
  0       2     1
# SIGMA0_DIH  SIGMA0_TOT
  6.0         6.0
# K_DIH    K_TOT      E_DIH    E_TOT
  0.0      0.001234   0.0      -1234.56
# EQUIL  WINDOW
  0      0
END
";
        let mut cursor = Cursor::new(input);
        let block = GamdBlock::parse(&mut cursor).unwrap().unwrap();

        assert_eq!(block.search_mode, SearchMode::NoSearch);
        assert_eq!(block.k_tot, Some(0.001234));
        assert_eq!(block.e_tot, Some(-1234.56));
    }

    #[test]
    fn test_roundtrip() {
        let block = GamdBlock {
            search_mode: SearchMode::GamdSearch,
            boost_form: BoostForm::DualBoost,
            threshold_type: ThresholdType::UpperBound,
            sigma0_dih: 6.0,
            sigma0_tot: 6.0,
            k_dih: None,
            k_tot: None,
            e_dih: None,
            e_tot: None,
            equilibration_steps: 500,
            window_size: 1000,
        };

        let output = block.to_string();
        let mut cursor = Cursor::new(output.as_bytes());
        let parsed = GamdBlock::parse(&mut cursor).unwrap().unwrap();

        assert_eq!(block.search_mode, parsed.search_mode);
        assert_eq!(block.boost_form, parsed.boost_form);
        assert_eq!(block.threshold_type, parsed.threshold_type);
    }
}
