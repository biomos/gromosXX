//! rep_ana - Replica Exchange Analysis Tool
//!
//! Analyzes REMD simulation data to compute exchange statistics,
//! acceptance rates, and replica diffusion.
//!
//! Usage: rep_ana @repdat replica.dat [@out exchange_stats.dat]

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::process;

#[derive(Debug, Clone)]
struct ExchangeEvent {
    step: usize,
    time: f64,
    replica_i: usize,
    replica_j: usize,
    temp_i: f64,
    temp_j: f64,
    probability: f64,
    accepted: bool,
}

#[derive(Debug)]
struct ExchangeStatistics {
    total_attempts: usize,
    total_accepted: usize,
    acceptance_rate: f64,
    pair_stats: Vec<Vec<PairStat>>,
}

#[derive(Debug, Clone)]
struct PairStat {
    attempts: usize,
    accepted: usize,
    acceptance_rate: f64,
}

fn print_usage() {
    eprintln!("rep_ana - Replica Exchange Analysis Tool");
    eprintln!();
    eprintln!("Analyzes REMD simulation data for exchange statistics");
    eprintln!();
    eprintln!("Usage: rep_ana @repdat <replica_file> [@out <output>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @repdat   Replica exchange data file (replica.dat)");
    eprintln!("  @out      Output statistics file (optional, default: exchange_stats.dat)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  rep_ana @repdat replica.dat");
    eprintln!("  rep_ana @repdat replica.dat @out stats.dat");
}

fn parse_replica_file(path: &str) -> Result<Vec<ExchangeEvent>, String> {
    let file = File::open(path)
        .map_err(|e| format!("Cannot open replica file {}: {}", path, e))?;
    let reader = BufReader::new(file);

    let mut events = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading line: {}", e))?;
        let trimmed = line.trim();

        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }

        // Expected format: step time rep_i rep_j T_i T_j prob accepted
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 8 {
            continue;
        }

        let step = parts[0].parse().map_err(|_| format!("Invalid step: {}", parts[0]))?;
        let time = parts[1].parse().map_err(|_| format!("Invalid time: {}", parts[1]))?;
        let replica_i = parts[2].parse().map_err(|_| format!("Invalid replica_i: {}", parts[2]))?;
        let replica_j = parts[3].parse().map_err(|_| format!("Invalid replica_j: {}", parts[3]))?;
        let temp_i = parts[4].parse().map_err(|_| format!("Invalid temp_i: {}", parts[4]))?;
        let temp_j = parts[5].parse().map_err(|_| format!("Invalid temp_j: {}", parts[5]))?;
        let probability = parts[6].parse().map_err(|_| format!("Invalid prob: {}", parts[6]))?;
        let accepted = parts[7].parse::<i32>().map_err(|_| format!("Invalid accepted: {}", parts[7]))? == 1;

        events.push(ExchangeEvent {
            step,
            time,
            replica_i,
            replica_j,
            temp_i,
            temp_j,
            probability,
            accepted,
        });
    }

    if events.is_empty() {
        return Err("No exchange events found in replica file".to_string());
    }

    Ok(events)
}

fn calculate_statistics(events: &[ExchangeEvent], num_replicas: usize) -> ExchangeStatistics {
    let total_attempts = events.len();
    let total_accepted = events.iter().filter(|e| e.accepted).count();
    let acceptance_rate = if total_attempts > 0 {
        total_accepted as f64 / total_attempts as f64
    } else {
        0.0
    };

    // Initialize pair statistics
    let mut pair_stats = vec![vec![PairStat { attempts: 0, accepted: 0, acceptance_rate: 0.0 }; num_replicas]; num_replicas];

    for event in events {
        let i = event.replica_i;
        let j = event.replica_j;

        pair_stats[i][j].attempts += 1;
        pair_stats[j][i].attempts += 1;

        if event.accepted {
            pair_stats[i][j].accepted += 1;
            pair_stats[j][i].accepted += 1;
        }
    }

    // Calculate acceptance rates for each pair
    for i in 0..num_replicas {
        for j in 0..num_replicas {
            if pair_stats[i][j].attempts > 0 {
                pair_stats[i][j].acceptance_rate =
                    pair_stats[i][j].accepted as f64 / pair_stats[i][j].attempts as f64;
            }
        }
    }

    ExchangeStatistics {
        total_attempts,
        total_accepted,
        acceptance_rate,
        pair_stats,
    }
}

fn write_statistics(path: &str, stats: &ExchangeStatistics, num_replicas: usize) -> Result<(), String> {
    let file = File::create(path)
        .map_err(|e| format!("Cannot create output file {}: {}", path, e))?;
    let mut writer = std::io::BufWriter::new(file);

    writeln!(writer, "# Replica Exchange Statistics")
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "#")
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# Overall Statistics:")
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "#   Total attempts: {}", stats.total_attempts)
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "#   Total accepted: {}", stats.total_accepted)
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "#   Acceptance rate: {:.2}%", stats.acceptance_rate * 100.0)
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "#")
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# Pair-wise Exchange Statistics:")
        .map_err(|e| format!("Write error: {}", e))?;
    writeln!(writer, "# {:>6} {:>6} {:>10} {:>10} {:>12}",
        "Rep_i", "Rep_j", "Attempts", "Accepted", "Rate")
        .map_err(|e| format!("Write error: {}", e))?;

    for i in 0..num_replicas {
        for j in (i+1)..num_replicas {
            if stats.pair_stats[i][j].attempts > 0 {
                writeln!(writer, "  {:6} {:6} {:10} {:10} {:12.2}%",
                    i, j,
                    stats.pair_stats[i][j].attempts,
                    stats.pair_stats[i][j].accepted,
                    stats.pair_stats[i][j].acceptance_rate * 100.0)
                    .map_err(|e| format!("Write error: {}", e))?;
            }
        }
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
    let mut replica_file = String::new();
    let mut output_file = "exchange_stats.dat".to_string();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@repdat" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: Missing value for @repdat");
                    print_usage();
                    process::exit(1);
                }
                replica_file = args[i].clone();
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

    if replica_file.is_empty() {
        eprintln!("Error: Missing required argument @repdat");
        print_usage();
        process::exit(1);
    }

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║              Replica Exchange Analysis Tool                  ║");
    println!("║          Exchange Statistics and Acceptance Rates            ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    // Read replica exchange file
    println!("Reading replica exchange data: {}", replica_file);
    let events = match parse_replica_file(&replica_file) {
        Ok(e) => e,
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    };
    println!("  Exchange events loaded: {}", events.len());
    println!();

    // Determine number of replicas
    let num_replicas = events.iter()
        .flat_map(|e| vec![e.replica_i, e.replica_j])
        .max()
        .unwrap_or(0) + 1;

    println!("Detected {} replicas", num_replicas);
    println!();

    // Calculate statistics
    println!("Calculating exchange statistics...");
    let stats = calculate_statistics(&events, num_replicas);

    println!("Overall Statistics:");
    println!("  Total attempts:   {}", stats.total_attempts);
    println!("  Total accepted:   {}", stats.total_accepted);
    println!("  Acceptance rate:  {:.2}%", stats.acceptance_rate * 100.0);
    println!();

    // Write output
    println!("Writing statistics: {}", output_file);
    if let Err(e) = write_statistics(&output_file, &stats, num_replicas) {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
    println!("  Statistics file written successfully");
    println!();

    println!("Analysis complete!");
}
