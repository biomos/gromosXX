///! gamd_ana - GaMD Analysis Tool
///!
///! Analyzes GaMD simulation data to recover unbiased free energies
///! via reweighting of the boost potential.
///!
///! Usage: gamd_ana @boost gamd_boost.dat @stats gamd_stats.dat [@pmf pmf.dat]

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

#[derive(Debug)]
struct GamdFrame {
    step: usize,
    time: f64,
    delta_v_total: f64,
    delta_v_dih: f64,
    delta_v_pot: f64,
}

#[derive(Debug)]
struct ReweightingStats {
    avg_delta_v: f64,
    max_delta_v: f64,
    min_delta_v: f64,
    std_delta_v: f64,
}

fn print_usage() {
    eprintln!("gamd_ana - GaMD Analysis Tool");
    eprintln!();
    eprintln!("Analyzes GaMD simulation data to recover unbiased free energies");
    eprintln!();
    eprintln!("Usage: gamd_ana @boost <boost_file> @stats <stats_file> [@pmf <output>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @boost    GaMD boost potential file (gamd_boost.dat)");
    eprintln!("  @stats    GaMD statistics file (gamd_stats.dat)");
    eprintln!("  @pmf      Output PMF file (optional, default: gamd_pmf.dat)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  gamd_ana @boost gamd_boost.dat @stats gamd_stats.dat");
    eprintln!("  gamd_ana @boost gamd_boost.dat @stats gamd_stats.dat @pmf free_energy.dat");
}

fn parse_boost_file(path: &str) -> Result<Vec<GamdFrame>, String> {
    let file = File::open(path)
        .map_err(|e| format!("Cannot open boost file {}: {}", path, e))?;
    let reader = BufReader::new(file);

    let mut frames = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading line: {}", e))?;
        let trimmed = line.trim();

        // Skip comments and empty lines
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 5 {
            continue;
        }

        let step = parts[0].parse().map_err(|_| format!("Invalid step: {}", parts[0]))?;
        let time = parts[1].parse().map_err(|_| format!("Invalid time: {}", parts[1]))?;
        let delta_v_total = parts[2].parse().map_err(|_| format!("Invalid ΔV: {}", parts[2]))?;
        let delta_v_dih = parts[3].parse().map_err(|_| format!("Invalid ΔV_dih: {}", parts[3]))?;
        let delta_v_pot = parts[4].parse().map_err(|_| format!("Invalid ΔV_tot: {}", parts[4]))?;

        frames.push(GamdFrame {
            step,
            time,
            delta_v_total,
            delta_v_dih,
            delta_v_pot,
        });
    }

    if frames.is_empty() {
        return Err("No data found in boost file".to_string());
    }

    Ok(frames)
}

fn calculate_reweighting_stats(frames: &[GamdFrame]) -> ReweightingStats {
    let mut sum = 0.0;
    let mut sum_sq = 0.0;
    let mut max_val = f64::NEG_INFINITY;
    let mut min_val = f64::INFINITY;

    for frame in frames {
        let dv = frame.delta_v_total;
        sum += dv;
        sum_sq += dv * dv;
        max_val = max_val.max(dv);
        min_val = min_val.min(dv);
    }

    let n = frames.len() as f64;
    let avg = sum / n;
    let variance = (sum_sq / n) - (avg * avg);
    let std = variance.sqrt();

    ReweightingStats {
        avg_delta_v: avg,
        max_delta_v: max_val,
        min_delta_v: min_val,
        std_delta_v: std,
    }
}

fn perform_cumulant_expansion(frames: &[GamdFrame], temperature: f64) -> f64 {
    // Cumulant expansion for free energy correction
    // ΔG ≈ -⟨ΔV⟩ + β/2 * [⟨ΔV²⟩ - ⟨ΔV⟩²]
    // where β = 1/(kT), k = 0.008314 kJ/(mol·K)

    let kb = 0.008314; // Boltzmann constant in kJ/(mol·K)
    let beta = 1.0 / (kb * temperature);

    let mut sum_dv = 0.0;
    let mut sum_dv_sq = 0.0;

    for frame in frames {
        let dv = frame.delta_v_total;
        sum_dv += dv;
        sum_dv_sq += dv * dv;
    }

    let n = frames.len() as f64;
    let avg_dv = sum_dv / n;
    let avg_dv_sq = sum_dv_sq / n;

    let delta_g = -avg_dv + (beta / 2.0) * (avg_dv_sq - avg_dv * avg_dv);

    delta_g
}

fn write_pmf(path: &str, frames: &[GamdFrame], free_energy_correction: f64) -> Result<(), String> {
    use std::io::Write;

    let file = File::create(path)
        .map_err(|e| format!("Cannot create PMF file {}: {}", path, e))?;
    let mut writer = std::io::BufWriter::new(file);

    writeln!(writer, "# GaMD Reweighted Free Energy")
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# Free energy correction: {:.6} kJ/mol", free_energy_correction)
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "#")
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# {:>8} {:>12} {:>12} {:>12}", "Step", "Time", "ΔV", "Weight")
        .map_err(|e| format!("Write error: {}", e))?;

    // Calculate weights: w_i = exp(β*ΔV_i)
    let kb = 0.008314;
    let temperature = 300.0; // TODO: Get from stats file
    let beta = 1.0 / (kb * temperature);

    for frame in frames {
        let weight = (beta * frame.delta_v_total).exp();
        writeln!(writer, "  {:8} {:12.6} {:12.6} {:12.6}",
            frame.step, frame.time, frame.delta_v_total, weight)
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
    let mut boost_file = String::new();
    let mut stats_file = String::new();
    let mut pmf_file = "gamd_pmf.dat".to_string();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@boost" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: Missing value for @boost");
                    print_usage();
                    process::exit(1);
                }
                boost_file = args[i].clone();
            }
            "@stats" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: Missing value for @stats");
                    print_usage();
                    process::exit(1);
                }
                stats_file = args[i].clone();
            }
            "@pmf" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: Missing value for @pmf");
                    print_usage();
                    process::exit(1);
                }
                pmf_file = args[i].clone();
            }
            _ => {
                eprintln!("Error: Unknown argument: {}", args[i]);
                print_usage();
                process::exit(1);
            }
        }
        i += 1;
    }

    if boost_file.is_empty() {
        eprintln!("Error: Missing required argument @boost");
        print_usage();
        process::exit(1);
    }

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                   GaMD Analysis Tool                         ║");
    println!("║          Reweighting for Unbiased Free Energies              ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    // Read boost potential file
    println!("Reading boost potential: {}", boost_file);
    let frames = match parse_boost_file(&boost_file) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    };
    println!("  Frames loaded: {}", frames.len());
    println!();

    // Calculate statistics
    println!("Analyzing boost potential distribution:");
    let stats = calculate_reweighting_stats(&frames);
    println!("  Average ΔV:   {:.4} kJ/mol", stats.avg_delta_v);
    println!("  Std dev ΔV:   {:.4} kJ/mol", stats.std_delta_v);
    println!("  Max ΔV:       {:.4} kJ/mol", stats.max_delta_v);
    println!("  Min ΔV:       {:.4} kJ/mol", stats.min_delta_v);
    println!();

    // Perform cumulant expansion
    let temperature = 300.0; // TODO: Read from stats file
    println!("Calculating free energy correction:");
    println!("  Temperature:  {:.1} K", temperature);
    let delta_g = perform_cumulant_expansion(&frames, temperature);
    println!("  ΔG (2nd order cumulant): {:.6} kJ/mol", delta_g);
    println!();

    // Write PMF
    println!("Writing reweighted data: {}", pmf_file);
    if let Err(e) = write_pmf(&pmf_file, &frames, delta_g) {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
    println!("  PMF file written successfully");
    println!();

    println!("Analysis complete!");
    println!();
    println!("Note: This is a simplified reweighting implementation.");
    println!("For production work, use exponential averaging or Maclaurin series.");
}
