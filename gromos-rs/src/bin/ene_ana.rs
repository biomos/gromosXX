//! ener_ana - Analyze energy trajectory files
//!
//! Usage: ener_ana @traj <energy_file> [@stats] [@plot]
//!
//! Analyzes GROMOS energy trajectory files (.tre) and provides statistics
//! including averages, standard deviations, minima, and maxima for all energy components.

use gromos_rs::io::energy::EnergyReader;
use std::env;
use std::process;

fn print_usage() {
    eprintln!("ener_ana - Energy trajectory analysis");
    eprintln!();
    eprintln!("Usage: ener_ana @traj <energy_file> [@stats] [@components <list>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @traj         Input energy trajectory file (.tre)");
    eprintln!("  @stats        Show detailed statistics (default: true)");
    eprintln!("  @components   Energy components to analyze (default: all)");
    eprintln!("                Available: kinetic, potential, total, temp, volume, pressure,");
    eprintln!("                          bond, angle, improper, dihedral, lj, coul, shake, restraint");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  ener_ana @traj output.tre");
    eprintln!("  ener_ana @traj output.tre @components kinetic,potential,total");
}

#[derive(Debug)]
struct EnergyStats {
    count: usize,
    mean: f64,
    std_dev: f64,
    min: f64,
    max: f64,
}

impl EnergyStats {
    fn from_values(values: &[f64]) -> Self {
        if values.is_empty() {
            return Self {
                count: 0,
                mean: 0.0,
                std_dev: 0.0,
                min: 0.0,
                max: 0.0,
            };
        }

        let count = values.len();
        let mean = values.iter().sum::<f64>() / count as f64;

        let variance = values.iter()
            .map(|v| (v - mean).powi(2))
            .sum::<f64>() / count as f64;
        let std_dev = variance.sqrt();

        let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        Self {
            count,
            mean,
            std_dev,
            min,
            max,
        }
    }

    fn print(&self, label: &str, unit: &str) {
        println!("  {:20} {:12.4} ± {:10.4} {} (min: {:10.4}, max: {:10.4})",
            label, self.mean, self.std_dev, unit, self.min, self.max);
    }
}

#[derive(Debug)]
struct AnalysisArgs {
    traj_file: String,
    show_stats: bool,
    components: Vec<String>,
}

fn parse_args(args: Vec<String>) -> Result<AnalysisArgs, String> {
    let mut traj_file = None;
    let mut show_stats = true;
    let mut components = Vec::new();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing value for @traj".to_string());
                }
                traj_file = Some(args[i].clone());
            }
            "@stats" => {
                show_stats = true;
            }
            "@components" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing value for @components".to_string());
                }
                components = args[i].split(',').map(|s| s.trim().to_string()).collect();
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    let traj_file = traj_file.ok_or("Missing required argument @traj")?;

    // If no components specified, use all
    if components.is_empty() {
        components = vec![
            "kinetic".to_string(),
            "potential".to_string(),
            "total".to_string(),
            "temperature".to_string(),
            "volume".to_string(),
            "pressure".to_string(),
            "bond".to_string(),
            "angle".to_string(),
            "improper".to_string(),
            "dihedral".to_string(),
            "lj".to_string(),
            "coul".to_string(),
            "shake".to_string(),
            "restraint".to_string(),
        ];
    }

    Ok(AnalysisArgs {
        traj_file,
        show_stats,
        components,
    })
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    // Parse arguments
    let parsed_args = match parse_args(args) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    println!("ener_ana - Energy trajectory analysis");
    println!("  Energy file: {}", parsed_args.traj_file);
    println!();

    // Open energy reader
    let mut reader = match EnergyReader::new(&parsed_args.traj_file) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error opening energy file: {}", e);
            process::exit(1);
        }
    };

    println!("  Title: {}", reader.title());
    println!();

    // Read all frames
    println!("  Reading energy frames...");
    let frames = match reader.read_all_frames() {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error reading energy frames: {}", e);
            process::exit(1);
        }
    };

    println!("  Total frames: {}", frames.len());
    println!();

    if frames.is_empty() {
        eprintln!("Warning: No energy frames found in file");
        process::exit(0);
    }

    // Collect values for each component
    let times: Vec<f64> = frames.iter().map(|f| f.time).collect();
    let kinetics: Vec<f64> = frames.iter().map(|f| f.kinetic).collect();
    let potentials: Vec<f64> = frames.iter().map(|f| f.potential).collect();
    let totals: Vec<f64> = frames.iter().map(|f| f.total).collect();
    let temperatures: Vec<f64> = frames.iter().map(|f| f.temperature).collect();
    let volumes: Vec<f64> = frames.iter().map(|f| f.volume).collect();
    let pressures: Vec<f64> = frames.iter().map(|f| f.pressure).collect();
    let bonds: Vec<f64> = frames.iter().map(|f| f.bond).collect();
    let angles: Vec<f64> = frames.iter().map(|f| f.angle).collect();
    let impropers: Vec<f64> = frames.iter().map(|f| f.improper).collect();
    let dihedrals: Vec<f64> = frames.iter().map(|f| f.dihedral).collect();
    let ljs: Vec<f64> = frames.iter().map(|f| f.lj).collect();
    let couls: Vec<f64> = frames.iter().map(|f| f.coul_real + f.coul_recip + f.coul_self).collect();
    let shakes: Vec<f64> = frames.iter().map(|f| f.shake).collect();
    let restraints: Vec<f64> = frames.iter().map(|f| f.restraint).collect();

    // Time range
    let time_stats = EnergyStats::from_values(&times);
    println!("Time range:");
    println!("  Start:  {:10.4} ps", time_stats.min);
    println!("  End:    {:10.4} ps", time_stats.max);
    println!("  Frames: {}", time_stats.count);
    println!();

    // Statistics for requested components
    if parsed_args.show_stats {
        println!("Energy statistics:");
        println!();

        for component in &parsed_args.components {
            match component.as_str() {
                "kinetic" => {
                    let stats = EnergyStats::from_values(&kinetics);
                    stats.print("Kinetic", "kJ/mol");
                }
                "potential" => {
                    let stats = EnergyStats::from_values(&potentials);
                    stats.print("Potential", "kJ/mol");
                }
                "total" => {
                    let stats = EnergyStats::from_values(&totals);
                    stats.print("Total", "kJ/mol");
                }
                "temperature" | "temp" => {
                    let stats = EnergyStats::from_values(&temperatures);
                    stats.print("Temperature", "K     ");
                }
                "volume" => {
                    let stats = EnergyStats::from_values(&volumes);
                    stats.print("Volume", "nm³   ");
                }
                "pressure" => {
                    let stats = EnergyStats::from_values(&pressures);
                    stats.print("Pressure", "bar   ");
                }
                "bond" => {
                    let stats = EnergyStats::from_values(&bonds);
                    stats.print("Bond", "kJ/mol");
                }
                "angle" => {
                    let stats = EnergyStats::from_values(&angles);
                    stats.print("Angle", "kJ/mol");
                }
                "improper" => {
                    let stats = EnergyStats::from_values(&impropers);
                    stats.print("Improper", "kJ/mol");
                }
                "dihedral" => {
                    let stats = EnergyStats::from_values(&dihedrals);
                    stats.print("Dihedral", "kJ/mol");
                }
                "lj" => {
                    let stats = EnergyStats::from_values(&ljs);
                    stats.print("Lennard-Jones", "kJ/mol");
                }
                "coul" | "coulomb" => {
                    let stats = EnergyStats::from_values(&couls);
                    stats.print("Coulomb (total)", "kJ/mol");
                }
                "shake" => {
                    let stats = EnergyStats::from_values(&shakes);
                    stats.print("SHAKE", "kJ/mol");
                }
                "restraint" => {
                    let stats = EnergyStats::from_values(&restraints);
                    stats.print("Restraint", "kJ/mol");
                }
                _ => {
                    eprintln!("Warning: Unknown component '{}'", component);
                }
            }
        }

        println!();
    }

    // Energy conservation check
    let total_drift = totals.last().unwrap() - totals.first().unwrap();
    let total_avg = totals.iter().sum::<f64>() / totals.len() as f64;
    let drift_percent = (total_drift / total_avg.abs()) * 100.0;

    println!("Energy conservation:");
    println!("  Initial total: {:12.4} kJ/mol", totals.first().unwrap());
    println!("  Final total:   {:12.4} kJ/mol", totals.last().unwrap());
    println!("  Drift:         {:12.4} kJ/mol ({:.6}%)", total_drift, drift_percent);
    println!();

    if drift_percent.abs() > 0.1 {
        println!("  WARNING: Energy drift exceeds 0.1%");
    } else {
        println!("  Energy is well conserved");
    }

    println!();
    println!("Analysis complete!");
}
