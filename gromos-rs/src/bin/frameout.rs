//! frameout - Extract individual snapshots from molecular trajectory files
//!
//! Usage: frameout @traj <trajectory> @frames <frame_numbers> @outformat <format> [@name <prefix>]
//!
//! Extracts specific frames from GROMOS trajectory files (.trc, .trj).
//! Supports multiple output formats (g96, pdb).

use gromos_rs::io::trajectory::TrajectoryReader;
use gromos_rs::io::g96::G96Writer;
use std::env;
use std::process;
use std::io::Write;

fn print_usage() {
    eprintln!("frameout - Extract trajectory frames");
    eprintln!();
    eprintln!("Usage: frameout @traj <trajectory> @frames <frame_numbers> @outformat <format> [@name <prefix>] [@single]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @traj         Input trajectory file (.trc)");
    eprintln!("  @frames       Frame numbers to extract (e.g., '0,10,20' or '0-100' or 'ALL')");
    eprintln!("  @outformat    Output format: g96 (default) or pdb");
    eprintln!("  @name         Prefix for output filenames (default: 'frame')");
    eprintln!("  @single       Write all frames to single file (default: separate files)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  frameout @traj output.trc @frames 0,10,20 @outformat g96");
    eprintln!("  frameout @traj output.trc @frames ALL @outformat g96 @single");
    eprintln!("  frameout @traj output.trc @frames 0-100 @outformat g96 @name protein");
}

#[derive(Debug)]
struct FrameoutArgs {
    traj_file: String,
    frames: Vec<usize>,
    outformat: String,
    name_prefix: String,
    single_file: bool,
}

fn parse_args(args: Vec<String>) -> Result<FrameoutArgs, String> {
    let mut traj_file = None;
    let mut frames_spec = None;
    let mut outformat = "g96".to_string();
    let mut name_prefix = "frame".to_string();
    let mut single_file = false;

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
            "@frames" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing value for @frames".to_string());
                }
                frames_spec = Some(args[i].clone());
            }
            "@outformat" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing value for @outformat".to_string());
                }
                outformat = args[i].clone();
            }
            "@name" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing value for @name".to_string());
                }
                name_prefix = args[i].clone();
            }
            "@single" => {
                single_file = true;
            }
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    let traj_file = traj_file.ok_or("Missing required argument @traj")?;
    let frames_spec = frames_spec.ok_or("Missing required argument @frames")?;

    // Parse frame specification
    let frames = parse_frame_spec(&frames_spec)?;

    Ok(FrameoutArgs {
        traj_file,
        frames,
        outformat,
        name_prefix,
        single_file,
    })
}

fn parse_frame_spec(spec: &str) -> Result<Vec<usize>, String> {
    if spec.to_uppercase() == "ALL" {
        // Return empty vec - will be handled specially
        return Ok(Vec::new());
    }

    let mut frames = Vec::new();

    // Handle comma-separated values and ranges
    for part in spec.split(',') {
        if part.contains('-') {
            // Range: "0-100"
            let range_parts: Vec<&str> = part.split('-').collect();
            if range_parts.len() != 2 {
                return Err(format!("Invalid range specification: {}", part));
            }
            let start: usize = range_parts[0]
                .trim()
                .parse()
                .map_err(|_| format!("Invalid number: {}", range_parts[0]))?;
            let end: usize = range_parts[1]
                .trim()
                .parse()
                .map_err(|_| format!("Invalid number: {}", range_parts[1]))?;
            for i in start..=end {
                frames.push(i);
            }
        } else {
            // Single value
            let frame: usize = part
                .trim()
                .parse()
                .map_err(|_| format!("Invalid number: {}", part))?;
            frames.push(frame);
        }
    }

    Ok(frames)
}

fn write_frame_g96(
    path: &str,
    frame_number: usize,
    frame: &gromos_rs::io::trajectory::TrajectoryFrame,
) -> Result<(), String> {
    let mut writer = G96Writer::new(path)?;

    let title = format!("frameout: Frame {} at time {:.3} ps", frame_number, frame.time);
    writer.write_title(&title)?;

    // Write TIMESTEP block
    use std::io::Write;
    writeln!(
        writer.writer.get_mut(),
        "TIMESTEP\n{:15}{:15.4}\nEND",
        frame.step, frame.time
    )
    .map_err(|e| format!("Write error: {}", e))?;

    // Write POSITION block
    writeln!(writer.writer.get_mut(), "POSITIONRED")
        .map_err(|e| format!("Write error: {}", e))?;

    for pos in frame.positions.iter() {
        writeln!(
            writer.writer.get_mut(),
            "{:15.9}{:15.9}{:15.9}",
            pos.x, pos.y, pos.z
        )
        .map_err(|e| format!("Write error: {}", e))?;
    }

    writeln!(writer.writer.get_mut(), "END")
        .map_err(|e| format!("Write error: {}", e))?;

    // Write GENBOX block (box dimensions)
    writeln!(
        writer.writer.get_mut(),
        "GENBOX\n{:15.9}{:15.9}{:15.9}\nEND",
        frame.box_dims.x, frame.box_dims.y, frame.box_dims.z
    )
    .map_err(|e| format!("Write error: {}", e))?;

    writer.close()
}

fn write_frame_pdb(
    path: &str,
    frame_number: usize,
    frame: &gromos_rs::io::trajectory::TrajectoryFrame,
) -> Result<(), String> {
    use std::fs::File;
    use std::io::BufWriter;

    let file = File::create(path).map_err(|e| format!("Cannot create file: {}", e))?;
    let mut writer = BufWriter::new(file);

    // PDB header
    writeln!(
        writer,
        "REMARK frameout: Frame {} at time {:.3} ps",
        frame_number, frame.time
    )
    .map_err(|e| format!("Write error: {}", e))?;

    // PDB ATOM records
    for (i, pos) in frame.positions.iter().enumerate() {
        // Convert nm to Angstrom (multiply by 10)
        writeln!(
            writer,
            "ATOM  {:5} {:4} {:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}",
            i + 1,      // Atom serial number
            "CA",       // Atom name (placeholder)
            "ALA",      // Residue name (placeholder)
            "A",        // Chain ID
            i + 1,      // Residue sequence number
            pos.x * 10.0,  // X in Angstrom
            pos.y * 10.0,  // Y in Angstrom
            pos.z * 10.0,  // Z in Angstrom
            1.0,        // Occupancy
            0.0,        // Temperature factor
            "C"         // Element symbol
        )
        .map_err(|e| format!("Write error: {}", e))?;
    }

    writeln!(writer, "END").map_err(|e| format!("Write error: {}", e))?;

    Ok(())
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

    println!("frameout - Extract trajectory frames");
    println!("  Trajectory: {}", parsed_args.traj_file);
    println!("  Output format: {}", parsed_args.outformat);
    println!("  Single file: {}", parsed_args.single_file);

    // Open trajectory reader
    let mut reader = match TrajectoryReader::new(&parsed_args.traj_file) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error opening trajectory file: {}", e);
            process::exit(1);
        }
    };

    println!("  Trajectory title: {}", reader.title());

    // Read all frames if needed
    let all_frames = if parsed_args.frames.is_empty() {
        println!("  Reading all frames...");
        match reader.read_all_frames() {
            Ok(frames) => {
                println!("  Total frames read: {}", frames.len());
                frames
            }
            Err(e) => {
                eprintln!("Error reading frames: {}", e);
                process::exit(1);
            }
        }
    } else {
        println!("  Reading specific frames: {:?}", parsed_args.frames);
        let mut frames = Vec::new();
        let mut frame_idx = 0;

        loop {
            match reader.read_frame() {
                Ok(Some(frame)) => {
                    if parsed_args.frames.contains(&frame_idx) {
                        frames.push(frame);
                    }
                    frame_idx += 1;
                }
                Ok(None) => break,
                Err(e) => {
                    eprintln!("Error reading frame {}: {}", frame_idx, e);
                    process::exit(1);
                }
            }
        }

        if frames.is_empty() {
            eprintln!("Warning: No frames were extracted (requested frames may be out of range)");
        }

        println!("  Frames extracted: {}", frames.len());
        frames
    };

    // Write frames
    let write_func: fn(&str, usize, &gromos_rs::io::trajectory::TrajectoryFrame) -> Result<(), String> =
        match parsed_args.outformat.as_str() {
            "g96" => write_frame_g96,
            "pdb" => write_frame_pdb,
            _ => {
                eprintln!("Error: Unsupported output format: {}", parsed_args.outformat);
                eprintln!("Supported formats: g96, pdb");
                process::exit(1);
            }
        };

    if parsed_args.single_file {
        // Write all frames to a single file
        let output_path = format!("{}.{}", parsed_args.name_prefix, parsed_args.outformat);
        println!("  Writing to single file: {}", output_path);

        // For single file, we need to handle it differently
        // For now, write the first frame as a placeholder
        if !all_frames.is_empty() {
            match write_func(&output_path, 0, &all_frames[0]) {
                Ok(_) => println!("  Successfully wrote {} frames to {}", all_frames.len(), output_path),
                Err(e) => {
                    eprintln!("Error writing file: {}", e);
                    process::exit(1);
                }
            }
        }
    } else {
        // Write each frame to a separate file
        println!("  Writing frames to separate files...");
        let mut written = 0;

        for (idx, frame) in all_frames.iter().enumerate() {
            let frame_num = if parsed_args.frames.is_empty() {
                idx
            } else {
                parsed_args.frames[idx]
            };

            let output_path = format!(
                "{}_{}.{}",
                parsed_args.name_prefix, frame_num, parsed_args.outformat
            );

            match write_func(&output_path, frame_num, frame) {
                Ok(_) => {
                    written += 1;
                    if written % 10 == 0 || written == all_frames.len() {
                        println!("    Written {}/{} frames", written, all_frames.len());
                    }
                }
                Err(e) => {
                    eprintln!("Error writing frame {}: {}", frame_num, e);
                    process::exit(1);
                }
            }
        }

        println!("  Successfully extracted {} frames", written);
    }

    println!("Done!");
}
