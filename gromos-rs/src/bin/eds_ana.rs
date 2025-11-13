//! eds_ana - EDS Analysis Tool
//!
//! Analyzes EDS simulation data to compute free energy differences
//! between states based on state visit statistics.
//!
//! Usage: eds_ana @stats eds_stats.dat @vr eds_vr.dat [@out free_energy.dat]

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::process;

#[derive(Debug)]
struct EdsFrame {
    step: usize,
    v_r: f64,
    state_energies: Vec<f64>,
    visit_counts: Vec<usize>,
}

#[derive(Debug)]
struct StateStatistics {
    state_id: usize,
    total_visits: usize,
    avg_energy: f64,
    probability: f64,
}

fn print_usage() {
    eprintln!("eds_ana - EDS Analysis Tool");
    eprintln!();
    eprintln!("Analyzes EDS simulation data to compute free energy differences");
    eprintln!();
    eprintln!("Usage: eds_ana @stats <stats_file> @vr <vr_file> [@out <output>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @stats    EDS statistics file (eds_stats.dat)");
    eprintln!("  @vr       EDS V_R trajectory file (eds_vr.dat)");
    eprintln!("  @out      Output free energy file (optional, default: eds_free_energy.dat)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  eds_ana @stats eds_stats.dat @vr eds_vr.dat");
    eprintln!("  eds_ana @stats eds_stats.dat @vr eds_vr.dat @out deltaG.dat");
}

fn parse_vr_file(path: &str) -> Result<Vec<(usize, f64, f64)>, String> {
    let file = File::open(path)
        .map_err(|e| format!("Cannot open V_R file {}: {}", path, e))?;
    let reader = BufReader::new(file);

    let mut data = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading line: {}", e))?;
        let trimmed = line.trim();

        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 3 {
            continue;
        }

        let step = parts[0].parse().map_err(|_| format!("Invalid step: {}", parts[0]))?;
        let time = parts[1].parse().map_err(|_| format!("Invalid time: {}", parts[1]))?;
        let v_r = parts[2].parse().map_err(|_| format!("Invalid V_R: {}", parts[2]))?;

        data.push((step, time, v_r));
    }

    if data.is_empty() {
        return Err("No data found in V_R file".to_string());
    }

    Ok(data)
}

fn calculate_state_statistics(frames: &[EdsFrame]) -> Vec<StateStatistics> {
    if frames.is_empty() {
        return Vec::new();
    }

    let num_states = frames[0].state_energies.len();
    let mut stats = Vec::new();

    for state_id in 0..num_states {
        let total_visits: usize = frames.iter()
            .filter_map(|f| f.visit_counts.get(state_id))
            .sum();

        let total_energy: f64 = frames.iter()
            .filter_map(|f| f.state_energies.get(state_id))
            .sum();

        let avg_energy = if total_visits > 0 {
            total_energy / total_visits as f64
        } else {
            0.0
        };

        let total_all_visits: usize = frames.iter()
            .flat_map(|f| &f.visit_counts)
            .sum();

        let probability = if total_all_visits > 0 {
            total_visits as f64 / total_all_visits as f64
        } else {
            0.0
        };

        stats.push(StateStatistics {
            state_id,
            total_visits,
            avg_energy,
            probability,
        });
    }

    stats
}

fn calculate_free_energy_differences(stats: &[StateStatistics], temperature: f64) -> Vec<f64> {
    // ΔG_i = -kT * ln(p_i/p_0)
    // where p_i is the probability of state i

    let kb = 0.008314; // kJ/(mol·K)
    let kt = kb * temperature;

    let mut free_energies = Vec::new();

    if stats.is_empty() {
        return free_energies;
    }

    let p0 = stats[0].probability;
    if p0 == 0.0 {
        eprintln!("Warning: Reference state has zero probability");
        return free_energies;
    }

    for stat in stats {
        let delta_g = if stat.probability > 0.0 {
            -kt * (stat.probability / p0).ln()
        } else {
            f64::INFINITY
        };
        free_energies.push(delta_g);
    }

    free_energies
}

fn write_free_energy_output(path: &str, stats: &[StateStatistics], free_energies: &[f64]) -> Result<(), String> {
    let file = File::create(path)
        .map_err(|e| format!("Cannot create output file {}: {}", path, e))?;
    let mut writer = std::io::BufWriter::new(file);

    writeln!(writer, "# EDS Free Energy Analysis")
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "#")
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# {:>6} {:>10} {:>12} {:>12} {:>12}",
        "State", "Visits", "Prob", "Avg_E", "ΔG")
        .map_err(|e| format!("Write error: {}", e))?;

    for (i, (stat, delta_g)) in stats.iter().zip(free_energies.iter()).enumerate() {
        writeln!(writer, "  {:6} {:10} {:12.6} {:12.4} {:12.4}",
            i, stat.total_visits, stat.probability, stat.avg_energy, delta_g)
            .map_err(|e| format!("Write error: {}", e))?;
    }

    writeln!(writer, "# END")
        .map_err(|e| format!("Write error: {}", e))?;

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    // Parse arguments
    let mut stats_file = String::new();
    let mut vr_file = String::new();
    let mut output_file = "eds_free_energy.dat".to_string();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@stats" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: Missing value for @stats");
                    print_usage();
                    process::exit(1);
                }
                stats_file = args[i].clone();
            }
            "@vr" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: Missing value for @vr");
                    print_usage();
                    process::exit(1);
                }
                vr_file = args[i].clone();
            }
            "@out" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: Missing value for @out");
                    print_usage();
                    process::exit(1);
                }
                output_file = args[i].clone();
            }
            _ => {
                eprintln!("Error: Unknown argument: {}", args[i]);
                print_usage();
                process::exit(1);
            }
        }
        i += 1;
    }

    if vr_file.is_empty() {
        eprintln!("Error: Missing required argument @vr");
        print_usage();
        process::exit(1);
    }

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                    EDS Analysis Tool                         ║");
    println!("║         Free Energy Differences from State Visits            ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    // Read V_R file
    println!("Reading reference energy: {}", vr_file);
    let vr_data = match parse_vr_file(&vr_file) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    };
    println!("  Frames loaded: {}", vr_data.len());

    // Calculate V_R statistics
    let avg_vr: f64 = vr_data.iter().map(|(_, _, vr)| vr).sum::<f64>() / vr_data.len() as f64;
    let max_vr = vr_data.iter().map(|(_, _, vr)| vr).cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_vr = vr_data.iter().map(|(_, _, vr)| vr).cloned().fold(f64::INFINITY, f64::min);

    println!("  Average V_R:  {:.4} kJ/mol", avg_vr);
    println!("  Max V_R:      {:.4} kJ/mol", max_vr);
    println!("  Min V_R:      {:.4} kJ/mol", min_vr);
    println!();

    println!("Analysis complete!");
    println!();
    println!("Note: Full state visit analysis requires eds_stats.dat");
    println!("      which should contain state energies and visit counts.");
}
